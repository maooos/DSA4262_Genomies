"""Produce a concise HTML results draft from completed Task 2 analysis tables."""
import argparse
import html
import json
import shutil
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, default=Path("outputs/task2/analysis"))
    parser.add_argument("--inventory", type=Path, default=Path("data/sgnex/inventory.json"))
    parser.add_argument("--publish-directory", type=Path,
                        help="Copy a compact, reviewable results bundle here (e.g. reports/task2)")
    args = parser.parse_args()
    path = args.analysis
    meta = json.loads((path / "analysis_metadata.json").read_text())
    inventory = json.loads(args.inventory.read_text())
    threshold = meta["threshold"]
    all_summary = pd.read_csv(path / "sample_summary.csv")
    samples = all_summary[(all_summary.min_reads == 20) & ((all_summary.threshold - threshold).abs() < 1e-12)]
    grouped = samples.groupby("cell_line").agg(runs=("sample_id", "count"), replicates=("replicate", "nunique"),
             eligible_min=("eligible_sites", "min"), eligible_max=("eligible_sites", "max"),
             positive_min=("positive_fraction", "min"), positive_max=("positive_fraction", "max"))
    display = pd.DataFrame({"Cell line": grouped.index, "Runs": grouped.runs,
        "Replicates": grouped.replicates,
        "Eligible sites/run": [f"{a:,}–{b:,}" for a, b in zip(grouped.eligible_min, grouped.eligible_max)],
        "Positive sites/run (%)": [f"{100*a:.2f}–{100*b:.2f}" for a, b in zip(grouped.positive_min, grouped.positive_max)]})
    agreement = pd.read_csv(path / "replicate_agreement_min20.csv")
    pairs = pd.read_csv(path / "cell_comparisons_min20.csv")
    if pairs.jaccard.notna().any():
        best = pairs.loc[pairs.jaccard.idxmax()]
        least = pairs.loc[pairs.jaccard.idxmin()]
        overlap_text = (f"The largest pairwise positive-site Jaccard similarity was {best.jaccard:.3f} "
                        f"for {best.cell_a} versus {best.cell_b}, based on {int(best.common_sites):,} shared eligible sites. "
                        f"The smallest was {least.jaccard:.3f} for {least.cell_a} versus {least.cell_b}, based on "
                        f"{int(least.common_sites):,} shared eligible sites. Pairwise eligible universes differ; these similarities are descriptive.")
    else:
        overlap_text = "No cell-line pair had an eligible nonempty positive union. Overlap cannot be estimated from this build."
    if agreement.spearman.notna().any():
        agreement_text = (f"Within-cell-line replicate Spearman correlations on shared covered sites ranged from "
                          f"{agreement.spearman.min():.3f} to {agreement.spearman.max():.3f} "
                          f"(median {agreement.spearman.median():.3f}).")
    else:
        agreement_text = "No replicate pair supported a score correlation in this build."
    candidate_rows = []
    for row in pairs.itertuples(index=False):
        compared = pd.read_csv(path / f"{row.cell_a}__{row.cell_b}_min20.csv.gz")
        top = compared[compared.consistent_difference].head(1)
        for candidate in top.to_dict("records"):
            candidate_rows.append({"Comparison": f"{row.cell_a} / {row.cell_b}",
                "Transcript": candidate["transcript_id"], "Position": candidate["transcript_position"],
                "Score A": candidate["mean_score_a"], "Score B": candidate["mean_score_b"],
                "Difference": candidate["delta_score"]})
    candidates = pd.DataFrame(candidate_rows, columns=["Comparison", "Transcript", "Position", "Score A", "Score B", "Difference"]).sort_values("Difference", key=abs, ascending=False).head(5)
    candidates.to_csv(path / "report_example_candidates.csv", index=False)
    sensitivity = []
    for coverage in (10, 20, 50):
        subset = all_summary[(all_summary.min_reads == coverage) & ((all_summary.threshold - threshold).abs() < 1e-12)]
        p = pd.read_csv(path / f"cell_comparisons_min{coverage}.csv")
        sensitivity.append({"Minimum reads/run": coverage, "Eligible sample–site observations": int(subset.eligible_sites.sum()),
            "Median pairwise Jaccard": p.jaccard.median(), "Median shared sites/pair": p.common_sites.median()})
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity.to_csv(path / "coverage_sensitivity.csv", index=False)
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Task 2 · SG-NEx results draft</title>
<style>body{{font:11pt Arial,sans-serif;color:#1e3437;max-width:950px;margin:40px auto;line-height:1.5}}h1{{font-size:22pt}}h2{{font-size:14pt;margin-top:28px}}table{{border-collapse:collapse;width:100%;font-size:9pt}}th,td{{text-align:left;padding:7px;border-bottom:1px solid #ddd}}th{{background:#eef3ee}}figure{{margin:22px 0;break-inside:avoid}}img{{width:100%}}figcaption,.note{{font-size:9pt;color:#596862}}code{{font-size:9pt;overflow-wrap:anywhere}}@media print{{body{{margin:0}}@page{{size:A4;margin:18mm}}}}</style>
<h1>Task 2: H04 predictions across SG-NEx</h1>
<p class="note">Results-section draft · generated from saved predictions · {html.escape(meta['created_at'])} · Status: {html.escape(meta['status'])}</p>
<h2>Data and analysis</h2>
<p>H04 was applied without retraining to {len(meta['samples'])} of {meta['expected_samples']} processed direct-RNA sequencing runs across {len(meta['cell_lines'])} cell lines. The live SG-NEx m6Anet inventory contains {inventory['processed_samples']} runs ({inventory['size_bytes']/1e9:.2f} GB); {inventory['directRNA_samples_without_processed_json']} additional direct-RNA runs in the metadata have no JSON input in that directory. Thus this analysis covers the available processed inventory, not every raw SG-NEx experiment.</p>
<p>The model generated {meta['total_scored_sites']:,} sample–site predictions. All sites were scored, with a primary analysis restricted to observations supported by at least 20 reads per run. Scores were thresholded at {threshold:.9f}. Technical-run scores passing the coverage filter were averaged within each named replicate, followed by an equally weighted mean across replicates. A site needed at least two observed replicates in each cell line for pairwise comparison. Sites were matched on transcript ID, position and 7-mer; missing observations were excluded, never treated as negative.</p>
<h2>Coverage and predicted-positive fractions</h2>
<p>Ensembl version suffixes were harmonized when at most one explicit version was observed for a stable transcript ID. Original IDs are preserved, and sites still require matching positions and 7-mers. The build harmonized {meta.get('identifier_harmonization', {}).get('changed_identifiers', 0):,} distinct input identifiers. The MCF7 naming-mismatch investigation is documented in the Task 2 guide. Known conflicting transcript versions are kept separate; unversioned identifiers still rely on the shared reference convention.</p>
<p>The primary analysis retained {int(samples.eligible_sites.sum()):,} sample–site observations. Run-level predicted-positive fractions ranged from {100*samples.positive_fraction.min():.2f}% to {100*samples.positive_fraction.max():.2f}%. These fractions use different eligible-site sets and therefore cannot establish a ranking of true cellular modification abundance.</p>
{display.to_html(index=False, border=0)}
<figure><img src="sample_overview.png" alt="Per-run coverage and predicted-positive fractions"><figcaption>Figure 1. Each point represents a sequencing run. Coverage filtering is applied consistently, but expression and sequencing depth still change which sites are observable.</figcaption></figure>
<h2>Reproducibility and between-cell-line overlap</h2>
<p>{agreement_text} Because sequence is an input to H04, score agreement may partly reflect shared sequence features and should not be interpreted as independent biological validation.</p>
<p>{overlap_text}</p>
<figure><img src="cell_line_overlap.png" style="max-width:600px" alt="Shared-site Jaccard heatmap"><figcaption>Figure 2. Jaccard similarity among positive predictions, restricted to sites meeting the replicate and coverage requirements in both cell lines. Diagonal cells are omitted.</figcaption></figure>
<h2>Candidate differences and sensitivity</h2>
<p>Candidate sites require all observed qualifying replicates in one cell line to be positive and all in the other to be negative, with at least two per group. The examples below are selected by absolute mean-score difference from the strongest candidate in each pair. They are hypotheses for follow-up, not statistically significant differential-modification calls.</p>
{candidates.to_html(index=False, border=0, float_format=lambda x: f'{x:.3f}')}
<p>Coverage sensitivity at the saved score threshold:</p>
{sensitivity.to_html(index=False, border=0, float_format=lambda x: f'{x:.3f}')}
<p>Additional score-threshold sensitivity results (0.15, saved threshold, 0.30 and 0.50) are in <code>sample_summary.csv</code> and the interactive explorer. Full pairwise tables retain every compared site; the explorer displays up to 500 candidates and a fixed random sample of 2,000 points per scatter plot, while computing statistics on all shared sites.</p>
<h2>Interpretation and limitations</h2>
<p>These are predictions from an uncalibrated classifier, not modification stoichiometry measurements. Differences may reflect expression, coverage, batch effects, sequence composition or transfer from the model's training distribution. Transcript isoforms are kept separate and may represent the same genomic nucleotide. The model metadata names data0/data1/data2 without source accessions; potential overlap with SG-NEx training-related samples remains unresolved. This analysis is therefore not an independent accuracy evaluation. Experimental independence of named replicate groups should be confirmed before inferential testing.</p>
<p>All raw-input SHA256 values are in download audit files, and the model checksum and prediction thresholds accompany each sample's output. Model SHA256: <code>{meta['model_sha256']}</code>.</p>
<p>Sources: <a href="https://github.com/GoekeLab/sg-nex-data">SG-NEx project and metadata</a>; <a href="https://github.com/GoekeLab/sg-nex-data/blob/master/docs/AWS_data_access_tutorial.md">processed-data access tutorial</a>. Code: this repository's <code>scripts/task2_*.py</code> and <code>src/task2_analysis.py</code>.</p></html>'''
    (path / "task2_results.html").write_text(body)
    if args.publish_directory:
        args.publish_directory.mkdir(parents=True, exist_ok=True)
        for name in ("task2_results.html", "sample_overview.png", "sample_overview.pdf",
                     "replicate_agreement.png", "replicate_agreement.pdf",
                     "cell_line_overlap.png", "cell_line_overlap.pdf", "analysis_metadata.json",
                     "sample_summary.csv", "cell_comparisons_min20.csv", "replicate_agreement_min20.csv",
                     "report_example_candidates.csv", "coverage_sensitivity.csv"):
            shutil.copy2(path / name, args.publish_directory / name)
        for name in ("manifest.csv", "inventory.json", "directRNA_without_processed_json.csv"):
            shutil.copy2(args.inventory.parent / name, args.publish_directory / name)
    print(path / "task2_results.html")


if __name__ == "__main__":
    main()
