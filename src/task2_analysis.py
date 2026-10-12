"""Coverage-aware, descriptive comparisons; missing sites never become negatives."""
from itertools import combinations

import numpy as np
import pandas as pd

SITE = ["transcript_id", "transcript_position", "sequence"]


def harmonize_transcripts(frame):
    """Match Ensembl versioned/unversioned identifiers without merging known revisions.

    Original identifiers remain in the observation table. If more than one
    explicit version is observed for a stable ID, that ID is NOT harmonized.
    Coordinates and the full 7-mer must still agree in every site comparison.
    """
    result = frame.copy()
    result["transcript_id_original"] = result.transcript_id
    identifiers = result.transcript_id.drop_duplicates()
    parsed = identifiers.str.extract(r"^(ENST\d+)(?:\.(\d+))?$")
    version_counts = parsed.dropna().groupby(0)[1].nunique()
    ambiguous = set(version_counts[version_counts > 1].index)
    mapping = {}
    for raw, base in zip(identifiers, parsed[0]):
        mapping[raw] = base if pd.notna(base) and base not in ambiguous else raw
    result["transcript_id"] = result.transcript_id.map(mapping)
    if result.duplicated(["sample_id", "transcript_id", "transcript_position"]).any():
        raise ValueError("Identifier harmonization creates duplicate sites within a sample; resolve reference versions explicitly.")
    audit = dict(changed_observations=int(result.transcript_id.ne(result.transcript_id_original).sum()),
                 changed_identifiers=sum(raw != base for raw, base in mapping.items()),
                 ambiguous_stable_ids=len(ambiguous),
                 rule="Strip ENST numeric version suffix only when at most one explicit version is observed; preserve original IDs and require matching position and 7-mer.")
    return result, audit


def replicate_scores(frame, min_reads):
    """Average qualifying run scores within a replicate, then weight replicates equally.

    Coverage is required in EACH included run. Technical run scores are averaged,
    not counted as independent replicates or interpreted as molecule proportions.
    """
    if min_reads < 1:
        raise ValueError("Minimum read coverage must be positive.")
    eligible = frame.loc[frame.n_reads >= min_reads]
    return eligible.groupby(["cell_line", "replicate", *SITE], observed=True, sort=False).agg(
        score=("score", "mean"), n_runs=("sample_id", "nunique"),
        n_reads=("n_reads", "sum"),
    ).reset_index()


def cell_scores(replicates, threshold, min_replicates=2):
    values = replicates.assign(positive=replicates.score.ge(threshold))
    result = values.groupby(["cell_line", *SITE], observed=True, sort=False).agg(
        mean_score=("score", "mean"), min_score=("score", "min"),
        max_score=("score", "max"), n_replicates=("replicate", "nunique"),
        positive_replicates=("positive", "sum"), total_reads=("n_reads", "sum"),
    ).reset_index()
    return result[result.n_replicates >= min_replicates].copy()


def compare_cells(cells, cell_a, cell_b, threshold):
    a = cells[cells.cell_line.eq(cell_a)].drop(columns="cell_line")
    b = cells[cells.cell_line.eq(cell_b)].drop(columns="cell_line")
    common = a.merge(b, on=SITE, how="inner", suffixes=("_a", "_b"), validate="one_to_one")
    common["delta_score"] = common.mean_score_a - common.mean_score_b
    common["consistent_difference"] = (
        (common.min_score_a.ge(threshold) & common.max_score_b.lt(threshold))
        | (common.min_score_b.ge(threshold) & common.max_score_a.lt(threshold))
    )
    pa, pb = common.mean_score_a.ge(threshold), common.mean_score_b.ge(threshold)
    union, intersection = int((pa | pb).sum()), int((pa & pb).sum())
    corr = common.mean_score_a.corr(common.mean_score_b, method="spearman") if len(common) > 1 else np.nan
    summary = dict(cell_a=cell_a, cell_b=cell_b, common_sites=len(common),
                   positive_a=int(pa.sum()), positive_b=int(pb.sum()),
                   positive_both=intersection, positive_a_only=int((pa & ~pb).sum()),
                   positive_b_only=int((pb & ~pa).sum()), negative_both=int((~pa & ~pb).sum()),
                   jaccard=intersection / union if union else None,
                   spearman=float(corr) if np.isfinite(corr) else None,
                   consistent_candidates=int(common.consistent_difference.sum()))
    return common.sort_values("delta_score", key=abs, ascending=False), summary


def replicate_agreement(replicates, threshold):
    records = []
    for cell, subset in replicates.groupby("cell_line", observed=True):
        for a, b in combinations(sorted(subset.replicate.unique()), 2):
            left = subset[subset.replicate.eq(a)][[*SITE, "score"]]
            right = subset[subset.replicate.eq(b)][[*SITE, "score"]]
            matched = left.merge(right, on=SITE, suffixes=("_a", "_b"), validate="one_to_one")
            pa, pb = matched.score_a.ge(threshold), matched.score_b.ge(threshold)
            corr = matched.score_a.corr(matched.score_b, method="spearman") if len(matched) > 1 else np.nan
            union = int((pa | pb).sum())
            records.append(dict(cell_line=cell, replicate_a=int(a), replicate_b=int(b),
                                common_sites=len(matched), spearman=float(corr) if np.isfinite(corr) else None,
                                jaccard=int((pa & pb).sum()) / union if union else None))
    return pd.DataFrame(records, columns=["cell_line", "replicate_a", "replicate_b", "common_sites", "spearman", "jaccard"])


def sample_summary(frame, sample, min_reads, threshold):
    qualified = frame[frame.n_reads >= min_reads]
    count = len(qualified)
    return {**sample, "min_reads": min_reads, "threshold": threshold,
            "scored_sites": len(frame), "eligible_sites": count,
            "positive_sites": int(qualified.score.ge(threshold).sum()),
            "positive_fraction": float(qualified.score.ge(threshold).mean()) if count else None,
            "median_reads": float(qualified.n_reads.median()) if count else None,
            "median_score": float(qualified.score.median()) if count else None}
