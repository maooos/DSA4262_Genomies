"""Export reproducible Task 2 tables, figures and a static explorer data bundle."""
import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from src.task2_analysis import SITE, cell_scores, compare_cells, harmonize_transcripts, replicate_agreement, replicate_scores, sample_summary


def records(frame):
    return json.loads(frame.to_json(orient="records"))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".partial")
    encoded = json.dumps(value, allow_nan=False, separators=(",", ":"))
    if path.suffix == ".gz":
        with gzip.open(temp, "wt", encoding="utf-8") as handle:
            handle.write(encoded)
    else:
        temp.write_text(encoded)
    temp.replace(path)


def export_transcripts(frame, directory, sample_ids):
    """Partition all >=10-read observations for exact transcript lookup, no subsampling."""
    indexed = frame.copy()
    indexed["sample_index"] = indexed.sample_id.map({s: i for i, s in enumerate(sample_ids)})
    indexed["shard"] = indexed.transcript_id.map(lambda s: hashlib.sha256(s.encode()).hexdigest()[:2])
    index = {}
    for shard, group in indexed.groupby("shard", sort=True):
        transcripts = {}
        for tx, sites in group.groupby("transcript_id", sort=False):
            columns = ["sample_index", "transcript_position", "sequence", "n_reads", "score"]
            if "transcript_id_original" in sites:
                columns.append("transcript_id_original")
            rows = sites[columns].copy()
            # Preserve float precision: rounding before a threshold comparison can
            # change the call for a score immediately beside the saved threshold.
            transcripts[tx] = list(rows.itertuples(index=False, name=None))
            index[tx] = shard
        write_json(directory / "transcripts" / f"{shard}.json.gz", transcripts)
        (directory / "transcripts" / f"{shard}.json").unlink(missing_ok=True)
    write_json(directory / "transcript_index.json", index)
    aliases = {}
    if "transcript_id_original" in frame:
        identifiers = frame[["transcript_id_original", "transcript_id"]].drop_duplicates()
        aliases = {raw: canonical for raw, canonical in identifiers.itertuples(index=False, name=None) if raw != canonical}
    write_json(directory / "transcript_aliases.json", aliases)


def plot_results(directory, summaries, agreement, pair_summaries, threshold):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    lines = sorted(summaries.cell_line.unique())
    for i, cell in enumerate(lines):
        subset = summaries[summaries.cell_line.eq(cell)]
        offsets = np.linspace(-0.16, 0.16, len(subset))
        axes[0].scatter(i + offsets, subset.eligible_sites, color="#087d83", s=35)
        axes[1].scatter(i + offsets, subset.positive_fraction * 100, color="#d56b37", s=35)
    for ax in axes:
        ax.set_xticks(range(len(lines)), lines, rotation=25)
        ax.grid(axis="y", alpha=0.2)
    axes[0].set(ylabel="Sites with at least 20 reads", title="Coverage by sequencing run")
    axes[1].set(ylabel="Predicted-positive sites (%)", title=f"Run-level calls at threshold {threshold:.4f}")
    fig.savefig(directory / "sample_overview.png", dpi=220)
    fig.savefig(directory / "sample_overview.pdf")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    for i, cell in enumerate(lines):
        subset = agreement[agreement.cell_line.eq(cell)]
        ax.scatter([i] * len(subset), subset.spearman, s=45, color="#087d83")
    ax.set_xticks(range(len(lines)), lines, rotation=25)
    ax.set(ylabel="Spearman correlation on shared sites", ylim=(-1, 1), title="Within-cell-line replicate agreement (≥20 reads/run)")
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(directory / "replicate_agreement.png", dpi=220)
    fig.savefig(directory / "replicate_agreement.pdf")
    plt.close(fig)
    matrix = pd.DataFrame(np.nan, index=lines, columns=lines)
    for row in pair_summaries:
        matrix.loc[row["cell_a"], row["cell_b"]] = row["jaccard"]
        matrix.loc[row["cell_b"], row["cell_a"]] = row["jaccard"]
    fig, ax = plt.subplots(figsize=(7, 6), layout="constrained")
    im = ax.imshow(matrix, vmin=0, vmax=1, cmap="YlGnBu")
    ax.set_xticks(range(len(lines)), lines, rotation=35)
    ax.set_yticks(range(len(lines)), lines)
    for i in range(len(lines)):
        for j in range(len(lines)):
            if np.isfinite(matrix.iloc[i, j]):
                ax.text(j, i, f"{matrix.iloc[i, j]:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title("Positive-site overlap on pairwise shared eligible sites\n≥2 replicates/cell line; ≥20 reads/run")
    fig.colorbar(im, ax=ax, label="Jaccard similarity")
    fig.savefig(directory / "cell_line_overlap.png", dpi=220)
    fig.savefig(directory / "cell_line_overlap.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("data/sgnex/manifest.csv"))
    parser.add_argument("--predictions", type=Path, default=Path("outputs/task2/predictions"))
    parser.add_argument("--output", type=Path, default=Path("outputs/task2/analysis"))
    parser.add_argument("--site-data", type=Path, default=Path("task2_site/dist/data"))
    parser.add_argument("--allow-partial", action="store_true", help="Explicitly label an incomplete pilot; never imply all samples were scored")
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest)
    frames, summaries, included, missing, histogram = [], [], [], [], []
    threshold, model_sha = None, None
    for sample in manifest.to_dict("records"):
        sample_id = sample["sample_id"]
        path = args.predictions / f"{sample_id}.parquet"
        if not path.exists():
            missing.append(sample_id)
            continue
        meta = json.loads(path.with_suffix(".metadata.json").read_text())
        if threshold is not None and (threshold != meta["prediction_threshold"] or model_sha != meta["model_sha256"]):
            raise ValueError("Cannot compare outputs from different model versions or thresholds.")
        threshold, model_sha = meta["prediction_threshold"], meta["model_sha256"]
        frame = pd.read_parquet(path)
        if frame.sample_id.nunique() != 1 or frame.sample_id.iloc[0] != sample_id or len(frame) != meta["site_count"]:
            raise ValueError(f"Prediction metadata mismatch: {sample_id}")
        if frame.duplicated(["sample_id", *SITE[:2]]).any():
            raise ValueError(f"Duplicate observations: {sample_id}")
        if not np.isfinite(frame.score).all() or not frame.score.between(0, 1).all():
            raise ValueError(f"Invalid scores: {sample_id}")
        included.append({k: sample[k] for k in ("sample_id", "cell_line", "replicate", "run")})
        for coverage in (10, 20, 50):
            for cutoff in sorted(set([0.15, threshold, 0.3, 0.5])):
                summaries.append(sample_summary(frame, included[-1], coverage, cutoff))
            counts, edges = np.histogram(frame.loc[frame.n_reads >= coverage, "score"], bins=np.linspace(0, 1, 21))
            histogram.append(dict(sample_id=sample_id, min_reads=coverage, counts=counts.tolist()))
        frames.append(frame[frame.n_reads >= 10].copy())
        print(f"Loaded {sample_id}: {len(frame):,} sites", flush=True)
    if missing and not args.allow_partial:
        raise ValueError("Missing predictions (use --allow-partial for a labelled pilot): " + ", ".join(missing))
    if not frames:
        raise ValueError("No predictions available.")
    all_sites, identifier_audit = harmonize_transcripts(pd.concat(frames, ignore_index=True))
    args.output.mkdir(parents=True, exist_ok=True)
    args.site_data.mkdir(parents=True, exist_ok=True)
    all_sites.to_parquet(args.output / "observations_min10.parquet", index=False)
    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(args.output / "sample_summary.csv", index=False)
    cell_lines = sorted(all_sites.cell_line.unique())
    scenarios, agreement_default, default_pairs = {}, None, None
    for coverage in (10, 20, 50):
        print(f"Computing replicate and cell-line comparisons at ≥{coverage} reads/run...", flush=True)
        reps = replicate_scores(all_sites, coverage)
        reps.to_parquet(args.output / f"replicate_scores_min{coverage}.parquet", index=False)
        agreement = replicate_agreement(reps, threshold)
        agreement.to_csv(args.output / f"replicate_agreement_min{coverage}.csv", index=False)
        cells = cell_scores(reps, threshold)
        cells.to_parquet(args.output / f"cell_scores_min{coverage}.parquet", index=False)
        pairs = []
        for a, b in combinations(cell_lines, 2):
            compared, stats = compare_cells(cells, a, b, threshold)
            stem = f"{a}__{b}_min{coverage}"
            compared.to_csv(args.output / f"{stem}.csv.gz", index=False)
            candidates = compared[compared.consistent_difference].head(500)
            scatter = compared.sample(n=min(2000, len(compared)), random_state=4262)
            write_json(args.site_data / "pairs" / f"{stem}.json", {
                "summary": stats, "candidates": records(candidates),
                "scatter": records(scatter[[*SITE, "mean_score_a", "mean_score_b"]]),
                "scatter_total": len(compared), "candidate_limit": 500,
            })
            pairs.append(stats)
        pd.DataFrame(pairs, columns=["cell_a", "cell_b", "common_sites", "positive_a", "positive_b",
                                    "positive_both", "positive_a_only", "positive_b_only", "negative_both",
                                    "jaccard", "spearman", "consistent_candidates"]).to_csv(
                                        args.output / f"cell_comparisons_min{coverage}.csv", index=False)
        scenarios[str(coverage)] = {"agreement": records(agreement), "pairs": pairs,
                                    "cell_sites": cells.groupby("cell_line").size().to_dict()}
        if coverage == 20:
            agreement_default, default_pairs = agreement, pairs
    print(f"Exporting {len(all_sites):,} observations for transcript lookup...", flush=True)
    export_transcripts(all_sites, args.site_data, [s["sample_id"] for s in included])
    conflicts = all_sites.groupby(SITE[:2]).sequence.nunique().gt(1).sum()
    payload = dict(created_at=datetime.now(timezone.utc).isoformat(),
                   status="partial" if missing else "complete_processed_inventory",
                   samples=included, missing_samples=missing, expected_samples=len(manifest),
                   cell_lines=cell_lines, threshold=threshold, model_sha256=model_sha,
                   summaries=summaries, histograms=histogram, scenarios=scenarios,
                   total_scored_sites=sum(x["scored_sites"] for x in summaries if x["min_reads"] == 20 and x["threshold"] == threshold),
                   transcript_export_min_reads=10, sequence_conflict_sites=int(conflicts),
                   identifier_harmonization=identifier_audit,
                   inventory=records(manifest))
    write_json(args.site_data / "summary.json", payload)
    write_json(args.output / "analysis_metadata.json", {k: v for k, v in payload.items() if k not in ("summaries", "histograms", "scenarios", "inventory")})
    default_summary = summary_frame[(summary_frame.min_reads == 20) & (summary_frame.threshold == threshold)]
    plot_results(args.output, default_summary, agreement_default, default_pairs, threshold)
    print(f"Analysed {len(included)}/{len(manifest)} samples; {len(all_sites):,} observations at ≥10 reads.", flush=True)


if __name__ == "__main__":
    main()
