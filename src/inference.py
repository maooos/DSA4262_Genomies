"""Score sites with saved modelling-notebook artifacts, without refitting."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import joblib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from threadpoolctl import threadpool_limits

from .data_parser import iter_sites
from .feature_engineering import (
    RAW_FEATURE_NAMES, REQUIRED_COLUMNS, SITE_KEY_COLUMNS, build_features,
)

KEYS = SITE_KEY_COLUMNS


def load_model_bundle(path):
    """Accept the dictionary saved by data0_modelling_pipeline.ipynb."""
    bundle = joblib.load(path)
    if not isinstance(bundle, dict) or not {"model", "feature_columns"} <= bundle.keys():
        raise ValueError(
            "Use a selected_model.joblib saved by data0_modelling_pipeline.ipynb "
            "(a dictionary containing model and feature_columns). "
            "newadjustmentfolder bundles use a different preprocessing pipeline."
        )
    if not hasattr(bundle["model"], "predict_proba"):
        raise ValueError("The saved model must support predict_proba.")
    if bundle.get("kmer_residuals", False):
        raise ValueError("This loader supports the original features without k-mer residuals.")
    columns = bundle["feature_columns"]
    if not columns or len(set(columns)) != len(columns):
        raise ValueError("Saved feature_columns must be nonempty and unique.")
    return bundle


def _validate_keys(frame):
    missing = set(KEYS) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing site identifiers: {sorted(missing)}")
    if frame[KEYS].isna().any().any() or frame[KEYS[0]].astype(str).str.strip().eq("").any():
        raise ValueError("Site identifiers must not be empty.")
    positions = pd.to_numeric(frame[KEYS[1]], errors="raise")
    if not np.isfinite(positions).all() or not (positions == np.floor(positions)).all():
        raise ValueError("Site positions must be finite integers.")


def _read_batches(path, batch_size):
    """Stream reads, carrying the last site across Parquet batch boundaries."""
    if path.suffix.lower() == ".parquet":
        pending = pd.DataFrame()
        for batch in pq.ParquetFile(path).iter_batches(
            batch_size=batch_size, columns=REQUIRED_COLUMNS,
        ):
            frame = pd.concat([pending, batch.to_pandas()], ignore_index=True)
            _validate_keys(frame)
            last_site = frame[KEYS].eq(frame[KEYS].iloc[-1]).all(axis=1)
            pending = frame.loc[last_site].copy()
            complete = frame.loc[~last_site]
            if not complete.empty:
                yield complete
        if not pending.empty:
            yield pending
        return
    if not path.name.lower().endswith((".json", ".jsonl", ".json.gz", ".jsonl.gz")):
        raise ValueError("Read input must be project JSON/JSONL (optionally gzip) or reads Parquet.")
    pieces, count = [], 0
    for transcript, position, sequence, values in iter_sites(path):
        frame = pd.DataFrame(values, columns=RAW_FEATURE_NAMES)
        frame[KEYS[0]], frame[KEYS[1]] = transcript, position
        frame["sequence"], frame["n_reads"] = sequence, len(values)
        pieces.append(frame)
        count += len(frame)
        if count >= batch_size:
            yield pd.concat(pieces, ignore_index=True)
            pieces, count = [], 0
    if pieces:
        yield pd.concat(pieces, ignore_index=True)


def features_from_reads(reads):
    """Reproduce the training notebook's engineered features and raw means."""
    _validate_keys(reads)
    missing = set(REQUIRED_COLUMNS) - set(reads.columns)
    if missing:
        raise ValueError(f"Missing read columns: {sorted(missing)}")
    if reads[REQUIRED_COLUMNS].isna().any().any():
        raise ValueError("Read data contains missing values.")
    if not np.isfinite(reads[RAW_FEATURE_NAMES].to_numpy(dtype=float)).all():
        raise ValueError("Read measurements must be finite.")
    if (reads[[name for name in RAW_FEATURE_NAMES if name.endswith('_dwell')]] <= 0).any().any():
        raise ValueError("Dwell times must be positive for the trained log transformation.")
    groups = reads.groupby(KEYS, sort=False)
    if not groups["n_reads"].nunique().eq(1).all() or not groups["sequence"].nunique().eq(1).all():
        raise ValueError("Each site must have one consistent sequence and read count.")
    if not groups.size().eq(groups["n_reads"].first()).all():
        raise ValueError("Incomplete or duplicate site reads: n_reads does not match the data.")
    engineered = build_features(reads).set_index(KEYS)
    raw_means = groups[RAW_FEATURE_NAMES].mean().add_prefix("raw_mean__")
    return engineered.join(raw_means).reset_index()


def _feature_batches(path, input_kind, batch_size):
    if input_kind == "reads":
        for reads in _read_batches(path, batch_size):
            yield features_from_reads(reads)
    elif input_kind == "features":
        if path.suffix.lower() == ".parquet":
            frame = pd.read_parquet(path)
            if set(KEYS) <= set(frame.index.names):
                frame = frame.reset_index()
            for start in range(0, len(frame), batch_size):
                yield frame.iloc[start:start + batch_size]
        elif path.suffix.lower() == ".csv":
            yield from pd.read_csv(path, chunksize=batch_size, dtype={KEYS[0]: str})
        else:
            raise ValueError("Feature input must be a site-level Parquet or CSV table.")
    else:
        raise ValueError("input_kind must be 'features' or 'reads'.")


def predict_to_csv(model_path, input_path, output_path, *, input_kind="features",
                   batch_size=200_000, threads=4):
    """Write one score per site observation and a JSON model/input record.

    Pooled feature tables retain dataset in their keys and output CSV, so an
    overlapping human site in data0/data1 remains two distinct observations.

    Existing outputs are preserved. A failed run removes its temporary files.
    Only load model artifacts that you trust, as joblib loads Python objects.
    """
    if batch_size < 1 or threads < 1:
        raise ValueError("batch_size and threads must be positive.")
    model_path, input_path, output_path = (
        Path(p).expanduser().resolve() for p in (model_path, input_path, output_path)
    )
    if output_path.suffix.lower() != ".csv":
        raise ValueError("Choose an output filename ending in .csv.")
    metadata_path = output_path.with_suffix(".metadata.json")
    for path in (output_path, metadata_path):
        if path.exists():
            raise FileExistsError(f"Choose a new output filename; already exists: {path}")
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    bundle = load_model_bundle(model_path)
    model, columns = bundle["model"], bundle["feature_columns"]
    classes = np.asarray(model.classes_)
    if classes.shape != (2,) or set(classes.tolist()) != {0, 1}:
        raise ValueError("Expected binary model classes 0 and 1.")
    positive_column = int(np.flatnonzero(classes == 1)[0])
    print(f"Model: {bundle.get('model_name', type(model).__name__)}", flush=True)
    print(f"Input: {input_path} ({input_kind})", flush=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.{uuid4().hex}.partial")
    temporary_metadata = temporary.with_suffix(".json")
    seen, count = set(), 0
    output_keys = None
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            for frame in _feature_batches(input_path, input_kind, batch_size):
                if frame.empty:
                    continue
                _validate_keys(frame)
                current_keys = ["dataset", *KEYS] if "dataset" in frame else KEYS
                if "dataset" in frame and (frame["dataset"].isna().any() or frame["dataset"].astype(str).str.strip().eq("").any()):
                    raise ValueError("Dataset identifiers must not be empty.")
                if output_keys is not None and current_keys != output_keys:
                    raise ValueError("Inconsistent dataset identifiers across batches.")
                output_keys = current_keys
                keys = list(frame[output_keys].itertuples(index=False, name=None))
                if len(set(keys)) != len(keys) or seen.intersection(keys):
                    raise ValueError("Duplicate sites found; expected one complete feature row per site.")
                seen.update(keys)
                missing = set(columns) - set(frame.columns)
                if missing:
                    raise ValueError(f"Dataset is missing model features: {sorted(missing)}")
                X = frame.loc[:, columns]
                if not np.isfinite(X.to_numpy(dtype=float)).all():
                    raise ValueError("Model features must be finite.")
                with threadpool_limits(limits=threads):
                    scores = np.asarray(model.predict_proba(X))[:, positive_column]
                if scores.shape != (len(frame),) or not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
                    raise ValueError("Model returned invalid scores; expected one value in [0, 1] per site.")
                result = frame[output_keys].copy()
                result["score"] = scores
                result.to_csv(handle, index=False, header=count == 0)
                count += len(result)
                print(f"Scored {count:,} sites", flush=True)
        if not count:
            raise ValueError("No sites were found in the input.")
        with model_path.open("rb") as handle:
            model_hash = hashlib.file_digest(handle, "sha256").hexdigest()
        stat = input_path.stat()
        metadata = {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "model_path": str(model_path), "model_sha256": model_hash,
            "model_name": bundle.get("model_name", type(model).__name__),
            "training_run_id": bundle.get("run_id"),
            "validation_threshold": bundle.get("threshold"),
            "input_path": str(input_path), "input_kind": input_kind,
            "input_size_bytes": stat.st_size, "input_mtime_ns": stat.st_mtime_ns,
            "output_path": str(output_path), "site_count": count,
            "feature_columns": columns,
            "identifier_columns": output_keys,
        }
        temporary_metadata.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        temporary_metadata.rename(metadata_path)
        temporary.rename(output_path)
    finally:
        temporary.unlink(missing_ok=True)
        temporary_metadata.unlink(missing_ok=True)
    print(f"Saved {count:,} site scores to {output_path}", flush=True)
    return metadata
