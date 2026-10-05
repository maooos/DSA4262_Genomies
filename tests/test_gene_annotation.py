"""Verify the full-read annotation stage and resumable raw-to-dataset workflow."""

import gzip
import json

import numpy as np
import pandas as pd
import pytest

from scripts.prepare_datasets import main as prepare_datasets
from src.data_parser import FEATURE_NAMES, file_sha256, parse_to_parquet
from src.feature_engineering import build_features
from src.gene_annotation import annotate_parquet


def make_inputs(tmp_path, *, dataset="data1"):
    raw_dir = tmp_path / "raw" / dataset
    raw_dir.mkdir(parents=True)
    signals = raw_dir / f"dataset{dataset[-1]}.json.gz"
    labels = raw_dir / "data.info"
    rows = []
    # Non-lexical transcript order and sites spanning multiple read batches.
    with gzip.open(signals, "wt") as handle:
        for i, (tx, count) in enumerate((("T9", 5), ("T1", 3))):
            values = np.arange(count * 9).reshape(count, 9) / 100 + 0.01
            handle.write(json.dumps({tx: {str(i + 10): {"AAGACCA": values.tolist()}}}) + "\n")
            rows.append({"transcript_id": tx, "transcript_position": i + 10,
                         "label": i, "sequence": "AAGACCA", "raw_n_reads": count,
                         "gene_id": f"ENSG{i + 1:011d}", "gene_name": f"GENE{i}",
                         "gene_id_status": "confirmed_ensembl91_mapping",
                         "sequence_matches_reference": True, "fasta_gene_matches_gtf": True})
    mapping = pd.DataFrame(rows)
    mapping[["transcript_id", "transcript_position", "label"]].to_csv(labels, index=False)
    map_path = tmp_path / "mapping.csv"
    # Deliberately reverse mapping order to catch positional joins.
    mapping.iloc[::-1].to_csv(map_path, index=False)
    parsed = tmp_path / "parsed.parquet"
    parse_to_parquet(signals, parsed, labels, batch_size=2)
    return signals, labels, parsed, map_path, mapping


def test_full_dataset_keeps_every_read_feature_label_and_order(tmp_path):
    _, _, parsed, mapping, _ = make_inputs(tmp_path)
    before = pd.read_parquet(parsed)
    output = tmp_path / "annotated.parquet"
    audit = annotate_parquet(parsed, mapping, output, dataset="data1", batch_size=3)
    after = pd.read_parquet(output)
    immutable = before.columns.drop("gene_id")
    pd.testing.assert_frame_equal(before[immutable], after[immutable])
    assert after.gene_id.tolist() == ["ENSG00000000001"] * 5 + ["ENSG00000000002"] * 3
    assert audit["gene_id_filled_reads"] == 8
    assert audit["gene_id_missing_after"] == 0
    assert audit["sites"] == 2
    pd.testing.assert_frame_equal(build_features(before), build_features(after))


@pytest.mark.parametrize("problem", ["duplicate", "missing_site", "wrong_sequence", "wrong_count", "conflict"])
def test_inconsistent_mapping_fails_without_final_output(tmp_path, problem):
    _, _, parsed, map_path, mapping = make_inputs(tmp_path)
    if problem == "duplicate":
        mapping = pd.concat([mapping, mapping.iloc[:1]])
    elif problem == "missing_site":
        mapping = mapping.iloc[:1]
    elif problem == "wrong_sequence":
        mapping.loc[0, "sequence"] = "TGGACTT"
    elif problem == "wrong_count":
        mapping.loc[0, "raw_n_reads"] += 1
    else:
        reads = pd.read_parquet(parsed)
        reads["gene_id"] = "ENSG99999999999"
        reads.to_parquet(parsed, index=False)
    mapping.to_csv(map_path, index=False)
    output = tmp_path / "annotated.parquet"
    with pytest.raises(ValueError):
        annotate_parquet(parsed, map_path, output, dataset="data1", batch_size=3)
    assert not output.exists()
    assert not output.with_suffix(".audit.json").exists()
    assert not list(tmp_path.glob("*.partial"))


def test_synthetic_gene_id_remains_not_applicable_with_reference_on_every_read(tmp_path):
    _, _, parsed, map_path, mapping = make_inputs(tmp_path, dataset="data2")
    mapping["gene_id"] = pd.NA
    mapping["gene_id_status"] = "not_applicable_synthetic_construct"
    mapping["geo_reference_id"] = "cc6m_2244_t7_ecorv"
    mapping["geo_reference_center_0based"] = [28, 35]
    mapping["reference_mapping_status"] = ["unique_7mer_in_curlcake_reference", "inferred_unique_ordered_alignment"]
    mapping.to_csv(map_path, index=False)
    output = tmp_path / "annotated.parquet"
    audit = annotate_parquet(parsed, map_path, output, dataset="data2", batch_size=3)
    result = pd.read_parquet(output)
    assert result.gene_id.isna().all()
    assert result.geo_reference_id.notna().all()
    assert result.reference_mapping_status.str.startswith("inferred").sum() == 3
    assert audit["gene_id_not_applicable_reads"] == 8
    assert audit["unresolved_applicable_gene_ids"] == 0


def test_pipeline_invokes_parse_then_annotation_and_resumes_by_checksum(tmp_path):
    signals, labels, _, _, mapping = make_inputs(tmp_path)
    processed = tmp_path / "processed"
    recovery = processed / "gene_id_recovery"
    recovery.mkdir(parents=True)
    map_path = recovery / "data1.gene_ids.csv"
    mapping.to_csv(map_path, index=False)
    recovery.joinpath("recovery.audit.json").write_text(json.dumps({
        "datasets": {"data1": {"signal_sha256": file_sha256(signals),
                                "metadata_sha256": file_sha256(labels)}},
        "outputs": [{"path": str(map_path), "sha256": file_sha256(map_path)}],
    }))
    args = ["--raw-dir", str(tmp_path / "raw"), "--processed-dir", str(processed), "--datasets", "data1"]
    first = prepare_datasets(args)
    output = processed / "annotated/data1_reads.parquet"
    fingerprint = file_sha256(output)
    second = prepare_datasets(args)
    assert first == second and fingerprint == file_sha256(output)
    # Changed raw inputs must invalidate the parse cache, never silently reuse it.
    with gzip.open(signals, "at") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="Stale/partial"):
        prepare_datasets(args)
