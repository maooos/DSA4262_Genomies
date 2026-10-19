import json
import threading
from functools import partial
from http.server import ThreadingHTTPServer
from urllib.error import HTTPError
from urllib.request import urlopen

import pandas as pd
import pytest

from scripts.task2_serve import DashboardHandler
from src.task2_analysis import (cell_scores, compare_cells, harmonize_transcripts,
                                replicate_agreement, replicate_scores, sample_summary)
from src.task2_query import CoverageQuery, MAX_READS, parse_reads


@pytest.fixture
def study(tmp_path):
    samples = [{"sample_id": f"run{i}", "cell_line": cell, "replicate": rep, "run": run}
               for i, (cell, rep, run) in enumerate([("A", 1, 1), ("A", 1, 2), ("A", 2, 1), ("B", 1, 1), ("B", 2, 1)])]
    rows = []
    for i, sample in enumerate(samples):
        for position, reads in enumerate([1, 24, 25, 37, 100], 1):
            rows.append({**sample, "transcript_id": "ENST000001.2" if i == 4 else "ENST000001",
                         "transcript_position": position, "sequence": "AAAAAAA" if not (i == 4 and position == 4) else "CCCCCCC",
                         "n_reads": reads if i != 1 else reads - 1 if reads > 1 else 1,
                         "score": [0.9, 0.1, 0.7, 0.2, 0.2][i] if position != 3 else 0.5})
    frame = pd.DataFrame(rows)
    canonical, _ = harmonize_transcripts(frame)
    metadata = {"created_at": "fixture", "model_sha256": "fixture-model", "threshold": 0.5,
                "samples": samples, "cell_lines": ["A", "B"]}
    data, assets = tmp_path / "queries", tmp_path / "assets"
    data.mkdir()
    assets.mkdir()
    canonical[["transcript_id_original", "transcript_id"]].drop_duplicates().rename(
        columns={"transcript_id_original": "original", "transcript_id": "canonical"}).to_parquet(data / "identifiers.parquet")
    bundled = []
    for sample in samples:
        subset = frame[frame.sample_id.eq(sample["sample_id"])]
        filename = sample["sample_id"] + ".parquet"
        subset.to_parquet(data / filename, index=False)
        bundled.append({**sample, "file": filename, "site_count": len(subset)})
    (assets / "summary.json").write_text(json.dumps(metadata))
    (data / "manifest.json").write_text(json.dumps({**metadata, "samples": bundled}))
    engine = CoverageQuery(data, assets)
    yield engine, canonical, samples, assets
    engine.close()


@pytest.mark.parametrize("reads", [1, 7, 25, 37, 51, MAX_READS])
def test_arbitrary_queries_match_original_scientific_analysis(study, reads):
    engine, frame, samples, _ = study
    reps = replicate_scores(frame, reads)
    cells = cell_scores(reps, 0.5)
    _, expected = compare_cells(cells, "A", "B", 0.5)
    result = engine.pair(reads, "A", "B")
    for key, value in expected.items():
        assert result["summary"][key] == pytest.approx(value) if isinstance(value, float) else result["summary"][key] == value
    overview = engine.overview(reads)
    assert overview["scenario"]["pairs"][0] == result["summary"]
    for sample in samples:
        for cutoff in [0.15, 0.3, 0.5]:
            expected = sample_summary(frame[frame.sample_id.eq(sample["sample_id"])], sample, reads, cutoff)
            actual = next(r for r in overview["summaries"] if r["sample_id"] == sample["sample_id"] and r["threshold"] == cutoff)
            assert actual == pytest.approx(expected)
    expected_agreement = json.loads(replicate_agreement(reps, 0.5).to_json(orient="records"))
    for expected, actual in zip(expected_agreement, overview["scenario"]["agreement"]):
        assert actual == pytest.approx(expected)


def test_transcripts_include_low_coverage_preserve_ids_and_do_not_round(study):
    engine, frame, _, _ = study
    result = engine.transcript("ENST000001.2")
    assert result["id"] == "ENST000001"
    assert len(result["rows"]) == len(frame)
    assert min(r[3] for r in result["rows"]) == 1
    assert result["rows"][-1][-1] == "ENST000001.2"
    assert engine.transcript("MISSING")["rows"] == []


@pytest.mark.parametrize("value", ["", "0", "-1", "2.5", "NaN", "1e3", MAX_READS + 1])
def test_invalid_coverage_rejected(value):
    with pytest.raises(ValueError, match="positive whole number"):
        parse_reads(value)


def test_http_validation_and_custom_comparison(study):
    engine, _, _, assets = study
    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(DashboardHandler, directory=str(assets)))
    server.engine, server.query_lock = engine, threading.Lock()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    try:
        with urlopen(base + "/api/pair?reads=25&a=A&b=B") as response:
            result = json.load(response)
            assert result["summary"]["common_sites"] == 2
        for query in ["reads=0", "reads=25.1", "reads=20&reads=50"]:
            with pytest.raises(HTTPError) as caught:
                urlopen(base + "/api/overview?" + query)
            assert caught.value.code == 400
        with urlopen(base + "/api/overview?reads=9007199254740991") as response:
            assert all(r["eligible_sites"] == 0 for r in json.load(response)["summaries"])
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
