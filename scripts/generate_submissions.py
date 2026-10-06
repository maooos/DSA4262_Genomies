"""Generate all three competition CSVs with one saved pooled H04 model."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from src.data_parser import file_sha256
from src.inference import KEYS, load_model_bundle
from src.modelling_data import DATASETS, OBSERVATION_KEYS, SPLIT_NAMES


def generate_submissions(run_dir, *, raw_dir=Path("data/raw"),
                         output_dir=Path("data/processed/submissions"),
                         team_name="Genomies", threads=4):
    """Score the union of train/validation/test features; never refit the model.

    Original metadata supplies only site identifiers and output order, not
    labels. Every dataset must have exactly the same sites as its metadata.
    All three CSVs are staged and validated before publishing any of them.
    """
    if not re.fullmatch(r"[A-Za-z0-9_-]+", team_name) or threads < 1:
        raise ValueError("Use a filename-safe team name and a positive thread count.")
    run_dir, raw_dir, output_dir = (
        Path(p).expanduser().resolve() for p in (run_dir, raw_dir, output_dir)
    )
    filenames = [f"{team_name}_dataset{i}.csv" for i in range(3)]
    manifest_name = f"{team_name}_submission_metadata.json"
    for name in [*filenames, manifest_name]:
        if (output_dir / name).exists():
            raise FileExistsError(f"Existing submission preserved: {output_dir / name}")

    model_path = run_dir / "selected_model.joblib"
    model_hash = file_sha256(model_path)
    bundle = load_model_bundle(model_path)
    if bundle.get("experiment", {}).get("experiment_id") != "H04":
        raise ValueError("The selected model must be H04.")
    if set(bundle.get("datasets", [])) != set(DATASETS):
        raise ValueError("Use the pooled H04 model developed on data0, data1 and data2.")
    model, columns = bundle["model"], bundle["feature_columns"]
    classes = np.asarray(model.classes_)
    if classes.shape != (2,) or set(classes.tolist()) != {0, 1}:
        raise ValueError("Expected binary model classes 0 and 1.")
    positive_column = int(np.flatnonzero(classes == 1)[0])

    feature_paths = [run_dir / f"X_{split}.parquet" for split in SPLIT_NAMES]
    parts = [pd.read_parquet(path) for path in feature_paths]
    if any(list(part.index.names) != OBSERVATION_KEYS for part in parts):
        raise ValueError("Expected pooled features indexed by dataset and site identifiers.")
    features = pd.concat(parts, verify_integrity=True)
    if set(features.index.get_level_values("dataset")) != set(DATASETS):
        raise ValueError("The feature tables must contain exactly data0, data1 and data2.")
    if not np.isfinite(features.loc[:, columns].to_numpy(dtype=float)).all():
        raise ValueError("All model features must be finite.")

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "team_name": team_name,
        "model_path": str(model_path), "model_sha256": model_hash,
        "model_name": bundle.get("model_name"), "run_id": bundle.get("run_id"),
        "training_partition": bundle.get("training_partition"),
        "prediction_method": "predict_proba for class 1; no thresholding or refitting",
        "coverage": "all sites across train, val and test; not a held-out evaluation",
        "feature_sources": [
            {"path": str(path), "sha256": file_sha256(path)} for path in feature_paths
        ],
        "files": [],
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".submission-", dir=output_dir.parent) as temporary:
        stage = Path(temporary)
        for dataset, filename in zip(DATASETS, filenames):
            metadata_path = raw_dir / dataset / ("data.info.labelled" if dataset == "data0" else "data.info")
            identifiers = pd.read_csv(metadata_path, usecols=KEYS, dtype={KEYS[0]: str})[KEYS]
            positions = pd.to_numeric(identifiers[KEYS[1]], errors="raise")
            if (identifiers.empty or identifiers.isna().any().any()
                    or identifiers[KEYS[0]].str.strip().eq("").any()
                    or not np.isfinite(positions).all()
                    or not positions.eq(np.floor(positions)).all()):
                raise ValueError(f"{dataset}: invalid site identifiers.")
            identifiers[KEYS[1]] = positions.astype("int64")
            expected = pd.MultiIndex.from_frame(identifiers)
            table = features.xs(dataset, level="dataset")
            if (not expected.is_unique or len(table) != len(expected)
                    or len(expected.difference(table.index))
                    or len(table.index.difference(expected))):
                raise ValueError(f"{dataset}: feature sites do not match all original metadata sites.")
            with threadpool_limits(limits=threads):
                scores = np.asarray(model.predict_proba(table.loc[expected, columns]))[:, positive_column]
            if (scores.shape != (len(expected),) or not np.isfinite(scores).all()
                    or not ((scores >= 0) & (scores <= 1)).all()):
                raise ValueError(f"{dataset}: invalid probability scores.")
            result = identifiers.assign(score=scores)
            path = stage / filename
            result.to_csv(path, index=False)
            saved = pd.read_csv(path, float_precision="round_trip", dtype={KEYS[0]: str})
            if list(saved.columns) != [*KEYS, "score"]:
                raise ValueError(f"{dataset}: incorrect submission columns.")
            pd.testing.assert_frame_equal(saved, result)
            manifest["files"].append({
                "dataset": dataset, "filename": filename, "site_count": len(result),
                "score_min": float(scores.min()), "score_max": float(scores.max()),
                "sha256": file_sha256(path),
                "identifier_source": str(metadata_path),
                "identifier_source_sha256": file_sha256(metadata_path),
            })
            print(f"Validated {filename}: {len(result):,} sites", flush=True)
        if file_sha256(model_path) != model_hash:
            raise ValueError("The model artifact changed during prediction.")
        (stage / manifest_name).write_text(json.dumps(manifest, indent=2) + "\n")
        output_dir.mkdir(parents=True, exist_ok=True)
        for name in [*filenames, manifest_name]:
            # Exclusive creation preserves any output created by another run.
            with (output_dir / name).open("xb") as target:
                target.write((stage / name).read_bytes())
    print(f"Saved three submission CSVs to {output_dir}", flush=True)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True, help="Completed pooled H04 run directory")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/submissions"))
    parser.add_argument("--team-name", default="Genomies")
    parser.add_argument("--threads", type=int, default=4)
    generate_submissions(**vars(parser.parse_args(argv)))


if __name__ == "__main__":
    main()
