"""Validate a completed static explorer against full-precision prediction tables."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("task2_site/dist/data"))
    parser.add_argument("--predictions", type=Path, default=Path("outputs/task2/predictions"))
    args = parser.parse_args()
    payload = json.loads((args.data / "summary.json").read_text())
    assert payload["status"] == "complete_processed_inventory", "An incomplete pilot cannot pass final validation."
    assert len(payload["samples"]) == payload["expected_samples"]
    threshold = payload["threshold"]
    total, count10, count20 = 0, 0, 0
    for sample in payload["samples"]:
        sample_id = sample["sample_id"]
        frame = pd.read_parquet(args.predictions / f"{sample_id}.parquet")
        assert not frame.duplicated(["transcript_id", "transcript_position"]).any()
        np.testing.assert_array_equal(frame.prediction, frame.score.ge(threshold).astype(int))
        total += len(frame)
        count10 += int(frame.n_reads.ge(10).sum())
        count20 += int(frame.n_reads.ge(20).sum())
        for coverage in (10, 20, 50):
            eligible = frame[frame.n_reads >= coverage]
            for summary in payload["summaries"]:
                if summary["sample_id"] == sample_id and summary["min_reads"] == coverage:
                    assert summary["eligible_sites"] == len(eligible)
                    assert summary["positive_sites"] == int(eligible.score.ge(summary["threshold"]).sum())
            histogram = next(h for h in payload["histograms"] if h["sample_id"] == sample_id and h["min_reads"] == coverage)
            assert sum(histogram["counts"]) == len(eligible)
    assert total == payload["total_scored_sites"]
    for scenario in payload["scenarios"].values():
        for pair in scenario["pairs"]:
            assert sum(pair[k] for k in ("positive_both", "positive_a_only", "positive_b_only", "negative_both")) == pair["common_sites"]
    index = json.loads((args.data / "transcript_index.json").read_text())
    aliases = json.loads((args.data / "transcript_aliases.json").read_text())
    assert set(aliases.values()) <= set(index)
    exported = 0
    for shard_path in (args.data / "transcripts").glob("*.json.gz"):
        with gzip.open(shard_path, "rt") as handle:
            shard = json.load(handle)
        for tx, rows in shard.items():
            assert index[tx] == shard_path.name.removesuffix(".json.gz")
            for row in rows:
                assert 0 <= row[0] < len(payload["samples"])
                assert row[3] >= 10 and 0 <= row[4] <= 1
                if len(row) > 5:
                    assert aliases.get(row[5], row[5]) == tx
            exported += len(rows)
    assert exported == count10
    print(json.dumps({"samples": len(payload["samples"]), "total_scored_observations": total,
                      "observations_min10": count10, "observations_min20": count20,
                      "transcripts_exported": len(index), "status": "validated"}, indent=2))


if __name__ == "__main__":
    main()
