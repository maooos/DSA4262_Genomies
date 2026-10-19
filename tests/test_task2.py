"""Scientific comparison invariants, including missingness and technical repeats."""
import json
import io
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from src.inference import predict_to_csv
from src.task2_analysis import cell_scores, compare_cells, harmonize_transcripts, replicate_agreement, replicate_scores, sample_summary
from src.sgnex import download_sample


def observation(cell, replicate, run, position, score, reads=25, sequence="AGAACAT"):
    return dict(sample_id=f"{cell}_{replicate}_{run}", cell_line=cell, replicate=replicate,
                run=run, transcript_id="tx", transcript_position=position,
                sequence=sequence, n_reads=reads, score=score)


def test_shared_sites_exclude_missing_and_low_coverage():
    rows = [observation(c, r, 1, 10, s) for c, s in (("A", .8), ("B", .1)) for r in (1, 2)]
    rows += [observation("A", r, 1, 20, .9) for r in (1, 2)]  # absent B
    rows += [observation(c, r, 1, 30, .9, 25 if c == "A" else 2) for c in ("A", "B") for r in (1, 2)]
    cells = cell_scores(replicate_scores(pd.DataFrame(rows), 20), .5)
    compared, summary = compare_cells(cells, "A", "B", .5)
    assert compared.transcript_position.tolist() == [10]
    assert summary["common_sites"] == 1
    assert summary["positive_a_only"] == 1
    assert summary["consistent_candidates"] == 1


def test_technical_runs_are_not_independent_replicates():
    data = pd.DataFrame([observation("A", 1, 1, 10, .9), observation("A", 1, 2, 10, .5)])
    reps = replicate_scores(data, 20)
    assert len(reps) == 1
    assert reps.n_runs.iloc[0] == 2
    assert reps.score.iloc[0] == .7
    assert cell_scores(reps, .5).empty


def test_equal_weight_replicates_and_sequence_matching():
    data = pd.DataFrame([observation("A", 1, 1, 10, .9), observation("A", 1, 2, 10, .9),
                         observation("A", 2, 1, 10, .1),
                         observation("B", 1, 1, 10, .1, sequence="TGAACAT"),
                         observation("B", 2, 1, 10, .1, sequence="TGAACAT")])
    cells = cell_scores(replicate_scores(data, 20), .5)
    assert cells.loc[cells.cell_line.eq("A"), "mean_score"].iloc[0] == .5
    compared, summary = compare_cells(cells, "A", "B", .5)
    assert compared.empty
    assert summary["jaccard"] is None
    assert summary["spearman"] is None


def test_empty_coverage_is_undefined_not_zero():
    summary = sample_summary(pd.DataFrame([observation("A", 1, 1, 10, .9, reads=2)]), {}, 20, .5)
    assert summary["eligible_sites"] == 0
    assert summary["positive_fraction"] is None


def test_single_replicate_pilot_has_explicit_empty_agreement_table():
    reps = replicate_scores(pd.DataFrame([observation("A", 1, 1, 10, .9)]), 20)
    agreement = replicate_agreement(reps, .5)
    assert agreement.empty
    assert list(agreement) == ["cell_line", "replicate_a", "replicate_b", "common_sites", "spearman", "jaccard"]


def test_identifier_versions_are_harmonized_but_conflicting_revisions_are_not():
    rows = [observation("A", i, 1, 10, .5) for i in range(1, 6)]
    frame = pd.DataFrame(rows)
    frame["transcript_id"] = ["ENST000001.2", "ENST000001", "ENST000002.1", "ENST000002.2", "ENST000002"]
    result, audit = harmonize_transcripts(frame)
    assert result.transcript_id.tolist() == ["ENST000001", "ENST000001", "ENST000002.1", "ENST000002.2", "ENST000002"]
    assert result.transcript_id_original.tolist() == frame.transcript_id.tolist()
    assert audit["changed_identifiers"] == 1 and audit["ambiguous_stable_ids"] == 1


def test_identifier_harmonization_rejects_within_sample_collisions():
    import pytest
    frame = pd.DataFrame([observation("A", 1, 1, 10, .5), observation("A", 1, 1, 10, .5)])
    frame["transcript_id"] = ["ENST000001.2", "ENST000001"]
    with pytest.raises(ValueError, match="duplicate sites"):
        harmonize_transcripts(frame)


def test_opt_in_metadata_preserves_original_scores(tmp_path):
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "scores.csv"
    predict_to_csv(root / "models/h04.joblib", root / "examples/sites.jsonl", output,
                   input_kind="reads", batch_size=3, threads=1, include_site_metadata=True)
    frame = pd.read_csv(output)
    original = pd.read_csv(root / "tests/fixtures/h04_scores.csv")
    np.testing.assert_allclose(frame.score, original.score, atol=1e-12, rtol=1e-12)
    assert frame.sequence.str.fullmatch("[ACGT]{7}").all()
    assert frame.n_reads.ge(1).all()
    assert json.loads(output.with_suffix(".metadata.json").read_text())["output_columns"] == list(frame)


def test_download_resume_checks_range_and_content(tmp_path, monkeypatch):
    import hashlib
    row = dict(sample_id="sample", url="https://example.test/data.json", size_bytes=6, etag="abc")
    directory = tmp_path / "sample"
    directory.mkdir()
    (directory / "data.json.partial").write_bytes(b"abc")
    (directory / "download.partial.json").write_text(json.dumps({k: row[k] for k in ("url", "size_bytes", "etag")}))

    class Response(io.BytesIO):
        status = 206
        headers = {"ETag": '"abc"', "Content-Range": "bytes 3-5/6"}

    def open_url(request, timeout):
        assert request.headers["Range"] == "bytes=3-"
        assert request.headers["If-match"] == '"abc"'
        return Response(b"def")

    monkeypatch.setattr("urllib.request.urlopen", open_url)
    target = download_sample(row, tmp_path)
    assert target.read_bytes() == b"abcdef"
    audit = json.loads((directory / "download.json").read_text())
    assert audit["sha256"] == hashlib.sha256(b"abcdef").hexdigest()
    assert not (directory / "data.json.partial").exists()


def test_partial_from_different_object_is_rejected(tmp_path):
    import pytest
    row = dict(sample_id="sample", url="https://example.test/data.json", size_bytes=6, etag="new")
    directory = tmp_path / "sample"
    directory.mkdir()
    (directory / "data.json.partial").write_bytes(b"old")
    (directory / "download.partial.json").write_text(json.dumps({**{k: row[k] for k in ("url", "size_bytes", "etag")}, "etag": "old"}))
    with pytest.raises(ValueError, match="unverified object"):
        download_sample(row, tmp_path)


def test_transcript_export_preserves_threshold_adjacent_scores(tmp_path):
    from scripts.task2_analyse import export_transcripts
    threshold = 0.20839324867970055
    values = [np.nextafter(threshold, 0), threshold, np.nextafter(threshold, 1)]
    frame = pd.DataFrame([observation("A", 1, 1, p, s) for p, s in enumerate(values)])
    export_transcripts(frame, tmp_path, ["A_1_1"])
    index = json.loads((tmp_path / "transcript_index.json").read_text())
    with gzip.open(tmp_path / "transcripts" / f"{index['tx']}.json.gz", "rt") as handle:
        rows = json.load(handle)["tx"]
    assert [r[-1] for r in rows] == values
    assert [r[-1] >= threshold for r in rows] == [False, True, True]


def test_docker_relocation_verifies_sha256_and_rejects_same_size_corruption(tmp_path):
    import hashlib
    import pytest
    from scripts.task2_predict import verify_existing_input
    source = tmp_path / 'data.json'
    source.write_bytes(b'original')
    audit = {'sample_id': 'sample', 'size_bytes': 8,
             'sha256': hashlib.sha256(b'original').hexdigest()}
    source.with_name('download.json').write_text(json.dumps(audit))
    metadata = {'input_path': '/old/machine/sample/data.json', 'input_size_bytes': 8,
                'input_mtime_ns': 1}
    verify_existing_input(metadata, source, 'sample')
    source.write_bytes(b'modified')
    with pytest.raises(ValueError, match='SHA256 mismatch'):
        verify_existing_input(metadata, source, 'sample')


def test_docker_relocation_rejects_missing_or_wrong_sample_audit(tmp_path):
    import pytest
    from scripts.task2_predict import verify_existing_input
    source = tmp_path / 'data.json'
    source.write_bytes(b'input')
    metadata = {'input_path': '/old/data.json', 'input_size_bytes': 5, 'input_mtime_ns': 1}
    with pytest.raises(ValueError, match='requires its download.json'):
        verify_existing_input(metadata, source, 'sample')
    source.with_name('download.json').write_text(json.dumps({'sample_id':'wrong', 'size_bytes':5}))
    with pytest.raises(ValueError, match='does not match sample'):
        verify_existing_input(metadata, source, 'sample')
