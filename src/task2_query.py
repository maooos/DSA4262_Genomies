"""Exact coverage queries over the bundled predictions; no inference or network I/O."""
import json
import math
import re
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory

import duckdb

MAX_READS = 9007199254740991  # Largest exact integer shared with JavaScript.


def parse_reads(value):
    if not re.fullmatch(r"[0-9]+", str(value)) or not 1 <= int(value) <= MAX_READS:
        raise ValueError("Enter a positive whole number of reads (1 or more).")
    return int(value)


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


class CoverageQuery:
    """One connection, used under the HTTP server's query lock; bounded result caches."""
    def __init__(self, directory, site_data):
        directory, site_data = Path(directory), Path(site_data)
        self.manifest = json.loads((directory / "manifest.json").read_text())
        self.study = json.loads((site_data / "summary.json").read_text())
        if any(self.manifest[k] != self.study[k] for k in ("created_at", "model_sha256", "threshold")):
            raise ValueError("Query data and dashboard snapshot differ. Regenerate the query bundle.")
        self.samples = self.manifest["samples"]
        if [s["sample_id"] for s in self.samples] != [s["sample_id"] for s in self.study["samples"]]:
            raise ValueError("Query sample order differs from the dashboard.")
        self.threshold = self.study["threshold"]
        self.cells = self.study["cell_lines"]
        self.temp = TemporaryDirectory(prefix="sgnex-query-")
        self.db = duckdb.connect(config={"threads": 2, "memory_limit": "1GB", "temp_directory": self.temp.name})
        self.db.from_parquet([str(directory / s["file"]) for s in self.samples]).create_view("raw")
        self.db.from_parquet(str(directory / "identifiers.parquet")).create_view("identifiers")
        self.db.execute("""CREATE VIEW observations AS SELECT r.* EXCLUDE(transcript_id),
            r.transcript_id AS transcript_id_original, i.canonical AS transcript_id
            FROM raw r JOIN identifiers i ON r.transcript_id = i.original""")
        self.coverage = None

    def close(self):
        self.db.close()
        self.temp.cleanup()

    def records(self, sql, parameters=None):
        result = self.db.execute(sql, parameters or [])
        columns = [c[0] for c in result.description]
        return [dict(zip(columns, map(clean, row))) for row in result.fetchall()]

    def aggregate(self, reads):
        reads = parse_reads(reads)
        if self.coverage == reads:
            return
        self.coverage = None
        self.db.execute("""CREATE OR REPLACE TEMP TABLE replicates AS
            SELECT cell_line, replicate, transcript_id, transcript_position, sequence,
                avg(score) AS score, count(*) AS n_runs, sum(n_reads) AS n_reads
            FROM observations WHERE n_reads >= ? GROUP BY ALL""", [reads])
        self.db.execute("""CREATE OR REPLACE TEMP TABLE cells AS
            SELECT cell_line, transcript_id, transcript_position, sequence,
                avg(score) AS mean_score, min(score) AS min_score, max(score) AS max_score,
                count(*) AS n_replicates, count(*) FILTER(WHERE score >= ?) AS positive_replicates,
                sum(n_reads) AS total_reads
            FROM replicates GROUP BY ALL HAVING count(*) >= 2""", [self.threshold])
        self.coverage = reads

    def common(self, a, b):
        columns = ["mean_score", "min_score", "max_score", "n_replicates", "positive_replicates", "total_reads"]
        fields = ", ".join(f"{side}.{c} AS {c}_{side}" for side in ("a", "b") for c in columns)
        self.db.execute(f"""CREATE OR REPLACE TEMP TABLE common AS
            SELECT a.transcript_id, a.transcript_position, a.sequence, {fields},
                a.mean_score - b.mean_score AS delta_score,
                ((a.min_score >= $cutoff AND b.max_score < $cutoff) OR
                 (b.min_score >= $cutoff AND a.max_score < $cutoff)) AS consistent_difference
            FROM cells a JOIN cells b USING(transcript_id, transcript_position, sequence)
            WHERE a.cell_line = $a AND b.cell_line = $b""", {"a": a, "b": b, "cutoff": self.threshold})

    def stats(self, table="common"):
        # Average ranks preserve ties, matching pandas/scipy Spearman correlation.
        return self.records(f"""WITH ranked AS (
            SELECT *, rank() OVER(ORDER BY mean_score_a) +
                (count(*) OVER(PARTITION BY mean_score_a) - 1) / 2.0 AS rank_a,
                rank() OVER(ORDER BY mean_score_b) +
                (count(*) OVER(PARTITION BY mean_score_b) - 1) / 2.0 AS rank_b
            FROM {table})
            SELECT count(*) AS common_sites,
                count(*) FILTER(WHERE mean_score_a >= $t) AS positive_a,
                count(*) FILTER(WHERE mean_score_b >= $t) AS positive_b,
                count(*) FILTER(WHERE mean_score_a >= $t AND mean_score_b >= $t) AS positive_both,
                count(*) FILTER(WHERE mean_score_a >= $t AND mean_score_b < $t) AS positive_a_only,
                count(*) FILTER(WHERE mean_score_b >= $t AND mean_score_a < $t) AS positive_b_only,
                count(*) FILTER(WHERE mean_score_a < $t AND mean_score_b < $t) AS negative_both,
                positive_both * 1.0 / nullif(positive_a + positive_b - positive_both, 0) AS jaccard,
                corr(rank_a, rank_b) AS spearman,
                count(*) FILTER(WHERE consistent_difference) AS consistent_candidates
            FROM ranked""", {"t": self.threshold})[0]

    @lru_cache(maxsize=8)
    def overview(self, reads):
        reads = parse_reads(reads)
        self.aggregate(reads)
        summaries = []
        for cutoff in sorted({0.15, self.threshold, 0.3, 0.5}):
            rows = {r["sample_id"]: r for r in self.records("""SELECT sample_id,
                count(*) AS eligible_sites, count(*) FILTER(WHERE score >= ?) AS positive_sites,
                median(n_reads) AS median_reads, median(score) AS median_score
                FROM raw WHERE n_reads >= ? GROUP BY sample_id""", [cutoff, reads])}
            for sample in self.samples:
                row = rows.get(sample["sample_id"], {"eligible_sites": 0, "positive_sites": 0,
                                                    "median_reads": None, "median_score": None})
                summaries.append({**{k: sample[k] for k in ("sample_id", "cell_line", "replicate", "run")},
                                  **row, "scored_sites": sample["site_count"], "min_reads": reads,
                                  "threshold": cutoff, "positive_fraction": row["positive_sites"] / row["eligible_sites"] if row["eligible_sites"] else None})
        pairs = []
        for a, b in combinations(self.cells, 2):
            self.common(a, b)
            pairs.append({"cell_a": a, "cell_b": b, **self.stats()})
        agreement = []
        for cell in self.cells:
            groups = sorted({s["replicate"] for s in self.samples if s["cell_line"] == cell})
            for a, b in combinations(groups, 2):
                self.db.execute("""CREATE OR REPLACE TEMP TABLE agreement AS
                    SELECT a.score AS mean_score_a, b.score AS mean_score_b, false AS consistent_difference
                    FROM replicates a JOIN replicates b USING(cell_line, transcript_id, transcript_position, sequence)
                    WHERE a.cell_line = ? AND a.replicate = ? AND b.replicate = ?""", [cell, a, b])
                stats = self.stats("agreement")
                agreement.append({"cell_line": cell, "replicate_a": a, "replicate_b": b,
                                  **{k: stats[k] for k in ("common_sites", "jaccard", "spearman")}})
        return {"summaries": summaries, "scenario": {"pairs": pairs, "agreement": agreement},
                "min_reads": reads}

    @lru_cache(maxsize=12)
    def pair(self, reads, a, b):
        reads = parse_reads(reads)
        if a not in self.cells or b not in self.cells or a == b:
            raise ValueError("Choose two different cell lines from this study.")
        self.aggregate(reads)
        self.common(a, b)
        summary = {"cell_a": a, "cell_b": b, **self.stats()}
        candidates = self.records("""SELECT * FROM common WHERE consistent_difference
            ORDER BY abs(delta_score) DESC, transcript_id, transcript_position, sequence LIMIT 500""")
        scatter = self.records("""SELECT transcript_id, transcript_position, sequence, mean_score_a, mean_score_b
            FROM common ORDER BY hash(transcript_id, transcript_position, sequence, 4262),
                transcript_id, transcript_position, sequence LIMIT 2000""")
        return {"summary": summary, "candidates": candidates, "scatter": scatter,
                "scatter_total": summary["common_sites"], "candidate_limit": 500}

    @lru_cache(maxsize=16)
    def transcript(self, identifier):
        if not identifier or len(identifier) > 160:
            raise ValueError("Enter a transcript identifier.")
        found = self.records("SELECT canonical FROM identifiers WHERE original = ? OR canonical = ? LIMIT 1", [identifier, identifier])
        if not found:
            return {"id": identifier, "rows": []}
        canonical = found[0]["canonical"]
        result = self.records("""SELECT sample_id, transcript_position, sequence, n_reads, score, transcript_id_original
            FROM observations WHERE transcript_id = ? ORDER BY sample_id, transcript_position""", [canonical])
        indices = {s["sample_id"]: i for i, s in enumerate(self.samples)}
        return {"id": canonical, "rows": [[indices[r["sample_id"]], r["transcript_position"], r["sequence"],
                r["n_reads"], r["score"], r["transcript_id_original"]] for r in result]}
