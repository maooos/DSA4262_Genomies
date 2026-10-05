import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# Order defined by the input-data specification.
FEATURE_NAMES = [
    "minus1_dwell",
    "minus1_signal_sd",
    "minus1_signal_mean",
    "central_dwell",
    "central_signal_sd",
    "central_signal_mean",
    "plus1_dwell",
    "plus1_signal_sd",
    "plus1_signal_mean",
]

SITE_KEYS = ["transcript_id", "transcript_position"]


def file_sha256(path):
    """Identify exact inputs and outputs for resumable data preparation."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_labels(label_path):
    """Load labels and validate that each site has one metadata row."""
    labels = pd.read_csv(
        label_path,
        dtype={
            "transcript_id": "string",
            "transcript_position": "Int64",
            "label": "Float64",
        },
    )

    required = [*SITE_KEYS, "label"]
    missing_columns = set(required) - set(labels.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}. "
            f"Found: {list(labels.columns)}"
        )

    # Required fields must be populated.
    if labels[required].isna().any().any():
        raise ValueError("Missing transcript ID, position, or label.")

    if labels["transcript_id"].str.strip().eq("").any():
        raise ValueError("Empty transcript ID.")

    if labels["label"].isna().any():
        raise ValueError("Labels contain missing values.")

    label_values = labels["label"].to_numpy(dtype=float)

    if not np.isfinite(label_values).all():
        raise ValueError("Labels contain non-finite values.")

    if not labels["label"].between(0, 1).all():
        raise ValueError("Expected labels in the range [0, 1].")

    duplicates = labels.duplicated(SITE_KEYS, keep=False)

    if duplicates.any():
        examples = labels.loc[duplicates, SITE_KEYS].head()
        raise ValueError(
            "Duplicate site keys in the label file:\n"
            + examples.to_string(index=False)
        )

    # Preserve gene IDs when supplied; otherwise use missing values.
    if "gene_id" in labels.columns:
        labels["gene_id"] = (
            labels["gene_id"]
            .astype("string")
            .replace(r"^\s*$", pd.NA, regex=True)
        )
    else:
        labels["gene_id"] = pd.Series(
            pd.NA,
            index=labels.index,
            dtype="string",
        )

    return labels[["gene_id", *SITE_KEYS, "label"]].copy()


def iter_sites(signal_path):
    """
    Yield one validated site at a time.

    Returns:
        transcript_id, transcript_position, sequence, reads_array

    reads_array has shape (number_of_reads, 9).
    """
    signal_path = Path(signal_path)
    opener = gzip.open if signal_path.suffix.lower() == ".gz" else open

    with opener(signal_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            # Ignore blank lines.
            if not line.strip():
                continue

            try:
                record = json.loads(line)

                if not isinstance(record, dict) or len(record) != 1:
                    raise ValueError(
                        "Expected exactly one transcript per JSON line."
                    )

                transcript_id, positions = next(iter(record.items()))

                if not transcript_id.strip():
                    raise ValueError("Empty transcript ID.")

                if not isinstance(positions, dict) or len(positions) != 1:
                    raise ValueError(
                        "Expected exactly one position per JSON line."
                    )

                position_text, sequences = next(iter(positions.items()))
                position = int(position_text)

                if not isinstance(sequences, dict) or len(sequences) != 1:
                    raise ValueError(
                        "Expected exactly one sequence per site."
                    )

                sequence, reads = next(iter(sequences.items()))

                if re.fullmatch(r"[ACGT]{7}", sequence) is None:
                    raise ValueError(
                        f"Expected a 7-mer using A/C/G/T; got {sequence!r}."
                    )

                features = np.asarray(reads, dtype=np.float64)

                if (
                    features.ndim != 2
                    or features.shape[0] == 0
                    or features.shape[1] != len(FEATURE_NAMES)
                ):
                    raise ValueError(
                        "Expected a nonempty reads × 9 array; "
                        f"got {features.shape}."
                    )

                if not np.isfinite(features).all():
                    raise ValueError(
                        "Features contain missing or non-finite values."
                    )

            except (ValueError, TypeError, OverflowError) as exc:
                raise ValueError(
                    f"{signal_path}, line {line_number}: {exc}"
                ) from exc

            yield transcript_id, position, sequence, features


def make_schema(include_labels):
    """Use consistent Parquet column types across all batches."""
    fields = [
        pa.field("transcript_id", pa.string()),
        pa.field("transcript_position", pa.int64()),
        pa.field("sequence", pa.string()),
        pa.field("central_5mer", pa.string()),
        pa.field("n_reads", pa.int64()),
        pa.field("read_idx", pa.int64()),
    ]

    fields.extend(
        pa.field(name, pa.float64())
        for name in FEATURE_NAMES
    )

    if include_labels:
        fields.extend([
            pa.field("gene_id", pa.string()),
            pa.field("label", pa.float64()),
        ])

    return pa.schema(fields)


def parse_to_parquet(
    signal_path,
    output_path,
    label_path=None,
    batch_size=100_000,
    max_sites=None,
):
    """
    Write a flat per-read Parquet table.

    Labels are optional so the same parser can process prediction data.
    Existing output files are not overwritten.
    """
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")

    if max_sites is not None and max_sites < 1:
        raise ValueError("max_sites must be positive.")

    output_path = Path(output_path)
    audit_path = output_path.with_suffix(".audit.json")
    temporary_path = output_path.with_name(output_path.name + ".partial")

    for path in [output_path, audit_path, temporary_path]:
        if path.exists():
            raise FileExistsError(
                f"{path} already exists. Choose a new output filename."
            )

    labels = load_labels(label_path) if label_path is not None else None
    label_keys = (
        set(labels[SITE_KEYS].itertuples(index=False, name=None))
        if labels is not None
        else None
    )

    schema = make_schema(include_labels=labels is not None)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seen_sites = set()
    batch = []
    site_count = 0
    read_count = 0
    unmatched_site_count = 0
    unexpected_motif_count = 0

    def write_batch(writer):
        """Attach optional metadata and write the current batch."""
        frame = pd.DataFrame(
            batch,
            columns=[
                *SITE_KEYS,
                "sequence",
                "central_5mer",
                "n_reads",
                "read_idx",
                *FEATURE_NAMES,
            ],
        )

        if labels is not None:
            frame = frame.merge(
                labels,
                on=SITE_KEYS,
                how="left",
                validate="many_to_one",
                sort=False,
            )

        table = pa.Table.from_pandas(
            frame,
            schema=schema,
            preserve_index=False,
            safe=True,
        )

        writer.write_table(table)
        batch.clear()

    # Write to a temporary file so a failed parse does not leave a complete final Parquet file.
    try:
        with pq.ParquetWriter(
            temporary_path,
            schema=schema,
            compression="snappy",
        ) as writer:
            for transcript_id, position, sequence, features in iter_sites(
                signal_path
            ):
                site_key = (transcript_id, position)

                if site_key in seen_sites:
                    raise ValueError(
                        f"Duplicate site in signal data: {site_key}"
                    )

                seen_sites.add(site_key)
                site_count += 1

                if label_keys is not None and site_key not in label_keys:
                    unmatched_site_count += 1

                # Report unexpected motifs, but do not remove these sites.
                if re.fullmatch(
                    r"[AGT][AG]AC[ACT]", sequence[1:6]
                ) is None:
                    unexpected_motif_count += 1

                central_5mer = sequence[1:6]
                n_reads = features.shape[0]

                for read_idx, values in enumerate(features):
                    batch.append((
                        transcript_id,
                        position,
                        sequence,
                        central_5mer,
                        n_reads,
                        read_idx,
                        *values.tolist(),
                    ))
                    read_count += 1

                    if len(batch) >= batch_size:
                        write_batch(writer)

                if max_sites is not None and site_count >= max_sites:
                    break

            if site_count == 0:
                raise ValueError("No sites were found in the input.")

            if batch:
                write_batch(writer)

        temporary_path.rename(output_path)

    except Exception:
        if temporary_path.exists():
            temporary_path.unlink()
        raise

    audit = {
        "signal_path": str(signal_path),
        "signal_sha256": file_sha256(signal_path),
        "label_path": str(label_path) if label_path is not None else None,
        "label_sha256": file_sha256(label_path) if label_path is not None else None,
        "output_path": str(output_path),
        "output_sha256": file_sha256(output_path),
        "max_sites_requested": max_sites,
        "parsed_sites": site_count,
        "parsed_read_rows": read_count,
        "sites_outside_expected_DRACH_motif": unexpected_motif_count,
        "label_rows": len(labels) if labels is not None else None,
        "sites_without_labels": (
            unmatched_site_count if labels is not None else None
        ),
        "labels_without_parsed_sites": (
            len(label_keys - seen_sites)
            if label_keys is not None
            else None
        ),
    }

    audit_path.write_text(
        json.dumps(audit, indent=2),
        encoding="utf-8",
    )

    return audit
