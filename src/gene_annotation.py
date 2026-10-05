"""Add recovered site annotations to complete per-read Parquets in batches."""

import json
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .data_parser import FEATURE_NAMES, SITE_KEYS, file_sha256


ANNOTATION_TYPES = {
    "gene_id_status": pa.string(),
    "gene_name": pa.string(),
    "transcript_version": pa.int64(),
    "gene_version": pa.int64(),
    "transcript_biotype": pa.string(),
    "seqname": pa.string(),
    "annotation_source": pa.string(),
    "sequence_matches_reference": pa.bool_(),
    "fasta_gene_matches_gtf": pa.bool_(),
    "geo_reference_id": pa.string(),
    "geo_reference_center_0based": pa.int64(),
    "epinano_reference_id": pa.string(),
    "epinano_reference_center_0based": pa.int64(),
    "sequence_only_candidate_count": pa.int64(),
    "reference_mapping_status": pa.string(),
    "reference_source": pa.string(),
}


def annotate_parquet(input_path, mapping_path, output_path, *, dataset, batch_size=100_000):
    """Fill only missing gene IDs; preserve all reads, labels and signal values.

    A complete site mapping is required. Unknown genes in human data cause an
    error. Synthetic constructs retain null gene_id plus their explicit status
    and reference assignments, including evidence strength.
    """
    if dataset not in {"data0", "data1", "data2"}:
        raise ValueError("dataset must be data0, data1 or data2")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    input_path, mapping_path, output_path = map(Path, (input_path, mapping_path, output_path))
    audit_path = output_path.with_suffix(".audit.json")
    temporary = output_path.with_name(output_path.name + ".partial")
    temporary_audit = audit_path.with_name(audit_path.name + ".partial")
    for path in (output_path, audit_path, temporary, temporary_audit):
        if path.exists():
            raise FileExistsError(f"{path} already exists; choose a new output path")

    mapping = pd.read_csv(mapping_path, dtype={"gene_id": "string", "seqname": "string"})
    required = [*SITE_KEYS, "gene_id", "gene_id_status", "sequence", "raw_n_reads", "label"]
    if missing := set(required) - set(mapping):
        raise ValueError(f"Missing annotation columns: {sorted(missing)}")
    if mapping.empty or mapping[SITE_KEYS].isna().any().any() or mapping.duplicated(SITE_KEYS).any():
        raise ValueError("Annotation must have one nonempty row per site")
    synthetic = dataset == "data2"
    expected_status = "not_applicable_synthetic_construct" if synthetic else "confirmed_ensembl91_mapping"
    if not mapping.gene_id_status.eq(expected_status).all():
        raise ValueError("Unexpected gene_id_status in annotation")
    if synthetic:
        if mapping.gene_id.notna().any():
            raise ValueError("Synthetic constructs must not be assigned biological gene IDs")
        for column in ("geo_reference_id", "reference_mapping_status", "geo_reference_center_0based"):
            if column not in mapping or mapping[column].isna().any():
                raise ValueError(f"Incomplete synthetic reference annotation: {column}")
    else:
        if not mapping.gene_id.str.fullmatch(r"ENSG\d{11}", na=False).all():
            raise ValueError("Missing or invalid human gene IDs")
        for column in ("sequence_matches_reference", "fasta_gene_matches_gtf"):
            if column not in mapping or not mapping[column].eq(True).all():
                raise ValueError(f"Reference verification failed: {column}")

    index = pd.MultiIndex.from_frame(mapping[SITE_KEYS])
    source = pq.ParquetFile(input_path)
    required_reads = {*SITE_KEYS, "gene_id", "sequence", "n_reads", "read_idx", "label", *FEATURE_NAMES}
    if missing := required_reads - set(source.schema_arrow.names):
        raise ValueError(f"Missing parsed-read columns: {sorted(missing)}")
    extra = [name for name in ANNOTATION_TYPES if name in mapping]
    if set(extra) & set(source.schema_arrow.names):
        raise ValueError("Input already has annotation columns; supply the original parsed Parquet")
    schema = pa.schema([*source.schema_arrow, *(pa.field(c, ANNOTATION_TYPES[c]) for c in extra)])
    counts = np.zeros(len(mapping), dtype=np.int64)
    rows = filled = missing_before = missing_after = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_digest, mapping_digest = file_sha256(input_path), file_sha256(mapping_path)
    try:
        with pq.ParquetWriter(temporary, schema, compression="snappy") as writer:
            for batch in source.iter_batches(batch_size=batch_size):
                table = pa.Table.from_batches([batch]).replace_schema_metadata(None)
                keys = table.select(SITE_KEYS).to_pandas()
                positions = index.get_indexer(pd.MultiIndex.from_frame(keys))
                if (positions < 0).any():
                    raise ValueError("Parsed data contains a site absent from the annotation")
                expected = mapping.iloc[positions].reset_index(drop=True)
                for actual_column, expected_column in (
                    ("sequence", "sequence"), ("n_reads", "raw_n_reads"), ("label", "label")
                ):
                    actual = table[actual_column].to_numpy()
                    if not np.array_equal(actual, expected[expected_column].to_numpy()):
                        raise ValueError(f"Parsed data and annotation disagree on {actual_column}")
                position_series = pd.Series(positions)
                expected_read_idx = (
                    position_series.groupby(position_series, sort=False).cumcount().to_numpy() + counts[positions]
                )
                if not np.array_equal(table["read_idx"].to_numpy(), expected_read_idx):
                    raise ValueError("Missing, duplicated or reordered read_idx within a site")
                counts += np.bincount(positions, minlength=len(mapping))
                old_gene = table["gene_id"].to_pandas().astype("string")
                new_gene = expected.gene_id
                conflict = old_gene.notna() & (new_gene.isna() | old_gene.ne(new_gene).fillna(False))
                if conflict.any():
                    raise ValueError("Recovered gene ID conflicts with an existing gene ID")
                replacement = old_gene.fillna(new_gene)
                missing_before += int(old_gene.isna().sum())
                missing_after += int(replacement.isna().sum())
                filled += int((old_gene.isna() & replacement.notna()).sum())
                table = table.set_column(table.schema.get_field_index("gene_id"), "gene_id",
                                         pa.array(replacement, type=pa.string(), from_pandas=True))
                for column in extra:
                    values = expected[column]
                    if pa.types.is_string(ANNOTATION_TYPES[column]):
                        values = values.astype("string")
                    table = table.append_column(column, pa.array(values, type=ANNOTATION_TYPES[column], from_pandas=True))
                writer.write_table(table)
                rows += len(batch)
        if not np.array_equal(counts, mapping.raw_n_reads.to_numpy()):
            raise ValueError("Parsed data does not contain every annotated read/site")

        # Read back the complete file; compare every original column except the
        # intentionally supplemented gene_id, preserving row order and types.
        immutable = [column for column in source.schema_arrow.names if column != "gene_id"]
        restored = pq.ParquetFile(temporary)
        for old, new in zip_longest(source.iter_batches(batch_size=batch_size),
                                    restored.iter_batches(batch_size=batch_size)):
            if old is None or new is None:
                raise ValueError("Output row count changed")
            before = pa.Table.from_batches([old]).select(immutable)
            after = pa.Table.from_batches([new]).select(immutable)
            if not before.equals(after):
                raise ValueError("Original signal, label, key or row order changed during annotation")
        if file_sha256(input_path) != input_digest or file_sha256(mapping_path) != mapping_digest:
            raise ValueError("Inputs changed during annotation")
        audit = {
            "dataset": dataset, "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "input_path": str(input_path), "input_sha256": input_digest,
            "mapping_path": str(mapping_path), "mapping_sha256": mapping_digest,
            "output_path": str(output_path), "output_sha256": file_sha256(temporary),
            "read_rows": rows, "sites": len(mapping), "transcripts": mapping.transcript_id.nunique(),
            "genes": mapping.gene_id.nunique(), "gene_id_missing_before": missing_before,
            "gene_id_filled_reads": filled, "gene_id_missing_after": missing_after,
            "gene_id_not_applicable_reads": rows if synthetic else 0,
            "unresolved_applicable_gene_ids": 0, "all_original_non_gene_columns_verified": True,
            "all_site_read_counts_verified": True, "added_columns": extra,
        }
        temporary_audit.write_text(json.dumps(audit, indent=2) + "\n")
        temporary.rename(output_path)
        temporary_audit.rename(audit_path)
        return audit
    except Exception:
        temporary.unlink(missing_ok=True)
        temporary_audit.unlink(missing_ok=True)
        raise
