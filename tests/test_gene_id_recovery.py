"""Guard against fabricated gene IDs and ambiguous short-sequence matches."""

import gzip
import json

import pandas as pd
import pytest

from scripts.recover_gene_ids import (
    load_sites,
    map_human_sites,
    map_synthetic_sites,
    unique_ordered_embedding,
    verify_mixture_signals,
)


def test_raw_metadata_join_preserves_unsorted_input_order(tmp_path):
    directory = tmp_path / "data1"
    directory.mkdir()
    info = pd.DataFrame({"transcript_id": ["T9", "T1"], "transcript_position": [10, 2],
                         "n_reads": [1, 1], "label": [0, 1]})
    info.to_csv(directory / "data.info", index=False)
    with gzip.open(directory / "dataset1.json.gz", "wt") as handle:
        for tx, pos in (("T1", 2), ("T9", 10)):
            handle.write(json.dumps({tx: {str(pos): {"AAGACCA": [[1.0] * 9]}}}) + "\n")
    result, _ = load_sites(tmp_path, "data1")
    pd.testing.assert_frame_equal(result[info.columns], info)


def test_repeated_contexts_need_a_unique_complete_alignment():
    assert unique_ordered_embedding(["A", "B", "A"], ["X", "A", "B", "A"]) == [1, 2, 3]
    with pytest.raises(ValueError, match="ambiguous"):
        unique_ordered_embedding(["A"], ["A", "B", "A"])
    with pytest.raises(ValueError, match="No complete"):
        unique_ordered_embedding(["B", "A"], ["A", "B"])


def test_human_mapping_preserves_original_gene_and_checks_center_coordinates():
    sites = pd.DataFrame([{"transcript_id": "T1", "transcript_position": 4,
                           "sequence": "AAGACCA", "gene_id": "G1"}])
    annotation = pd.DataFrame([{"transcript_id": "T1", "gene_id": "G1"}])
    reference = {"T1": ("T1.2 gene:G1.3", "CAAGACCATT")}
    result = map_human_sites(sites, annotation, reference)
    assert result.sequence_matches_reference.all()
    assert result.fasta_gene_matches_gtf.all()
    with pytest.raises(ValueError, match="conflict"):
        map_human_sites(sites.assign(gene_id="OTHER"), annotation, reference)
    with pytest.raises(ValueError, match="FASTA gene"):
        map_human_sites(sites, annotation, {"T1": ("T1 gene:OTHER", "CAAGACCATT")})


def test_synthetic_references_are_never_presented_as_human_gene_ids():
    sites = pd.DataFrame([
        {"transcript_id": tx, "transcript_position": p, "sequence": s}
        for tx in ("tx_id_0", "tx_id_6")
        for p, s in ((0, "AAGACCA"), (10, "TGGACTT"))
    ])
    geo = {"cc1": "AAGACCA", "cc2": "TGGACTT"}
    epinano = {k.upper(): "CCC" + s + "GG" for k, s in geo.items()}
    result, candidates, details = map_synthetic_sites(sites, geo, epinano)
    assert len(result) == 4
    assert result.gene_id.isna().all()
    assert result.gene_id_status.eq("not_applicable_synthetic_construct").all()
    assert result.geo_reference_center_0based.eq(3).all()
    assert result.epinano_reference_center_0based.eq(6).all()
    assert details["complete_ordered_mapping_is_unique"]
    assert len(candidates) == 2
    sites.loc[sites.transcript_id == "tx_id_6", "sequence"] = "AAAAAAA"
    with pytest.raises(ValueError, match="same ordered contexts"):
        map_synthetic_sites(sites, geo, epinano)


def test_mixture_matching_requires_all_nine_values_and_the_same_site(tmp_path):
    a, b = [1.0] * 9, [2.0] * 9
    changed = a.copy()
    changed[-1] += 0.00001
    records = [
        {"tx_id_0": {"0": {"AAGACCA": [a]}}},
        {"tx_id_6": {"0": {"AAGACCA": [b]}}},
        {"tx_id_1": {"0": {"AAGACCA": [a, b, changed]}}},
        {"tx_id_1": {"10": {"AAGACCA": [a]}}},
    ]
    path = tmp_path / "signals.json.gz"
    with gzip.open(path, "wt") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    result = {row["transcript_id"]: row for row in verify_mixture_signals(path)}
    assert result["tx_id_1"]["from_endpoint_1"] == 1
    assert result["tx_id_1"]["from_endpoint_0"] == 1
    assert result["tx_id_1"]["unmatched"] == 2
