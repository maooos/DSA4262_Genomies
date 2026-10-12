# Task 2: SG-NEx analysis

The frozen H04 model scores each available processed SG-NEx direct-RNA run.
The workflow produces per-site predictions, coverage-aware comparisons, report
figures, a results draft, and a static interactive explorer.

## Reproduce the analysis

From the repository root, using Python 3.12:

```bash
source .venv/bin/activate
python -m pip install -r requirements-task2.txt

# Save a live S3 inventory and the matching sample metadata, then download.
python -m scripts.sgnex_download --download --workers 3

# Score all inventoried samples with the saved model.
python -m scripts.task2_predict

# Compute comparisons and create the explorer's data assets.
python -m scripts.task2_analyse
python -m scripts.task2_validate
python -m scripts.task2_report --publish-directory reports/task2

# Local explorer: visit http://localhost:8765
python -m http.server 8765 --directory task2_site/dist
```

Downloads resume with `--reuse-inventory --download`. They use public HTTPS,
require no AWS account, verify the inventoried ETag and byte count, and record a
SHA256 digest. The October 12, 2026 inventory contains 23 runs across seven cell
lines, totalling 33.37 GB. Allow roughly 40 GB of free disk space for inputs,
predictions, and analysis. An interrupted download retains `.partial` data.

`task2_predict` preserves existing outputs, verifies input/model metadata before
reusing predictions, and starts missing samples. To score while downloads are
still running, use `--wait-seconds 14400` in a second terminal. A failed sample
stops the prediction job and does not leave a completed CSV. Rerun the same
command after resolving the error. Use a new output directory when changing
the model or inputs; do not mix model versions in one analysis.

For a smaller pilot, `sgnex_download --samples SAMPLE_ID ...` selects inputs.
Create a corresponding manifest subset for `task2_predict --manifest PATH`,
or use `task2_analyse --allow-partial` to explicitly label incomplete results.
The default analysis refuses missing predictions. Keep pilot analysis outputs
in a separate directory with `--output outputs/task2/pilot_analysis`.

## Outputs

| Location | Contents |
| --- | --- |
| `data/sgnex/manifest.csv` | Actual JSON URLs, sample IDs, cell lines, replicates, runs, sizes, ETags |
| `data/sgnex/directRNA_without_processed_json.csv` | Metadata samples without available processed m6Anet JSON |
| `data/sgnex/raw/<sample>/download.json` | Download source and SHA256 provenance |
| `outputs/task2/predictions/` | All scored sites, including low-coverage sites; CSV, metadata JSON, Parquet |
| `outputs/task2/analysis/observations_min10.parquet` | Combined site/sample table retaining observations with ≥10 reads |
| `outputs/task2/analysis/sample_summary.csv` | Run summaries at coverage 10/20/50 and four score thresholds |
| `outputs/task2/analysis/replicate_scores_min*.parquet` | Technical-run scores averaged within named replicates |
| `outputs/task2/analysis/cell_scores_min*.parquet` | Cell-line summaries supported by at least two replicates |
| `outputs/task2/analysis/*__*_min*.csv.gz` | Every shared eligible site for each cell-line pair |
| `outputs/task2/analysis/replicate_agreement_min*.csv` | Within-cell-line correlations and overlap |
| `outputs/task2/analysis/*.png` and `*.pdf` | Report figures |
| `outputs/task2/analysis/task2_results.html` | Results-section draft; adapt for the combined five-page assignment report |
| `task2_site/dist/` | Static explorer and partitioned result assets |
| `reports/task2/` | Compact results draft, figures, summary tables and inventory for GitHub review |

## Comparison rules

1. Score each run using identical H04 features and the saved validation
   threshold (`0.20839324867970055`). Preserve all model outputs.
2. Filter observations by reads **per run**. The primary analysis uses 20;
   10 and 50 are sensitivity settings, not validated biological cutoffs.
3. Within each named replicate, average the scores of qualifying technical
   runs. Runs failing coverage are omitted. This is score aggregation, not
   pooling reads and rerunning the nonlinear model.
4. Average replicate scores equally to form cell-line scores. Require at
   least two observed replicates in each cell line for pairwise comparisons.
5. Harmonize numeric Ensembl transcript version suffixes only when at most one
   explicit version is observed for that stable ID. Keep original identifiers
   in `transcript_id_original`; never merge known conflicting versions or
   duplicate observations within a sample. Match sites by the harmonized
   transcript ID, transcript position, and 7-mer. Use only
   their intersection. Missing/low-coverage sites are never negative calls.
6. A consistent candidate has all observed qualifying replicate scores on
   one side of the threshold in one cell line and on the opposite side in
   the other. Rank by absolute mean-score difference; no p-values are claimed.

The overview varies score thresholds only for run-level sensitivity. Cell-line
comparisons consistently use the saved threshold. Jaccard is positive-set
intersection / union within the pair's common eligible universe; it is missing
when both positive sets are empty. Spearman is missing for too few observations
or constant scores. The transcript explorer preserves full-precision scores.

The explorer partitions all ≥10-read observations by transcript for fast lookup.
Partitions are gzip-compressed and decoded on demand in the browser; compression
preserves all observations and full score precision. Use a current browser with
the standard `DecompressionStream` API for transcript lookup.
Scatter plots show up to 2,000 deterministically sampled sites and candidate
tables show up to 500 rows; their summary statistics use the full comparison.
Full pair tables are retained locally. No gene mapping or genome-coordinate
conversion is inferred from transcript identifiers.

MCF7 replicate 3 uses versioned Ensembl IDs, while replicate 4 uses unversioned
IDs. Before harmonization their exact transcript overlap was zero. After removing
the naming difference, all 969,788 shared raw coordinate pairs had identical
7-mers; 88,680 pairs had at least 20 reads in both runs. This check supports the
coordinate convention for these inputs. Unversioned identifiers still depend on
the shared reference convention and cannot prove full transcript-sequence identity.

Generated explorer data under `task2_site/dist/data/` is ignored by the main
repository to keep it compact. Regenerate it with `task2_analyse`; include that
directory when publishing the static Site. The code and small report bundle
remain reviewable in the main repository.

## Scope and interpretation

The live processed-data inventory has 23 samples, whereas the sample metadata
contains 56 direct-RNA runs. The other 33 need upstream signal processing before
H04 can accept them. The older `samples_with_RNAmod_data.tsv` misses newer
processed H9 and HEYA8 runs and contains outdated filenames; the workflow uses
the live S3 keys instead. This is complete for the available processed inputs,
not a claim that every raw SG-NEx run was analysed.

Scores are uncalibrated model outputs, not measured modification fractions.
Raw positive fractions mix biological and coverage/expression effects. Sequence
is a model input, so agreement can partly reflect shared sequence features.
Named replicate groups follow SG-NEx naming; verify experimental independence
before formal testing. Different transcript isoforms can represent the same
genomic nucleotide and are not collapsed here.

The training metadata lists data0/data1/data2 without source accessions. Potential
training overlap is unresolved. Do not present this unlabelled application as an
independent accuracy evaluation or calculate AUROC/AUPRC without ground truth.

Source: [SG-NEx](https://github.com/GoekeLab/sg-nex-data) and its
[data access tutorial](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/AWS_data_access_tutorial.md).

## Verification and AI-use record

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
node --check task2_site/dist/app.js
```

The tests check original H04 scores, complete-site batch boundaries, optional
metadata, missingness, coverage filtering, sequence matching and replicate
aggregation. Actual SG-NEx inputs also exercise the pipeline end to end.

Record this Codex session in the assignment's AI-use table: it implemented
inventory/download, batch scoring, descriptive analysis, tests, report generation
and the explorer. Two mistakes corrected during implementation were (1) filtering
the metadata using `direct-RNA` instead of its actual `directRNA` protocol value,
and (2) rounding scores for browser export, which can change calls beside the
saved threshold. Verified metadata values and full-precision transcript scores
were used instead. See [the AI-use record](AI_USE_TASK2.md). Review the scientific
interpretation and adapt the generated draft before submission.
