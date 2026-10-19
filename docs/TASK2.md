# Task 2: SG-NEx analysis

The frozen H04 model scores each available processed SG-NEx direct-RNA run.
The workflow produces per-site predictions, coverage-aware comparisons, report
figures, a results draft, and an interactive explorer with a local query service.

## Start the dashboard (recommended for teammates)

The large data folders are excluded from Git. On a fresh clone, first restore
`task2_site/dist/data/` and `task2_site/query_data/` from a teammate, or regenerate
them with the full pipeline below. The existing local checkout already has them.
With both folders present and Docker Desktop running, from the repository root:

```bash
docker compose up --build -d
```

Open **http://localhost:8080**. The initial image build needs internet; using the
built dashboard does not. There are no external chart libraries, fonts, CDNs,
accounts or runtime package installs. `docker compose down` stops it. Add
`PORT=8081` to a `.env` file to use another host port. The service binds to your
own computer's loopback interface.

The images support Apple Silicon and Intel/AMD through multi-platform base
images, pinned by manifest digest. Numerical package versions match the saved
H04 model. The dashboard uses nginx plus a small Python/DuckDB service for
custom coverage queries. The larger scientific Python image is built only when
the optional analysis service is used. Queries run locally over saved predictions;
no database account or separate database installation is needed.

Keep the website source, Docker setup, model and compact reports in Git.
Share `task2_site/dist/data/` (about 112 MB) and `task2_site/query_data/`
(about 340 MB) separately, preserving their paths. The latter includes all read
coverages in 23 Parquet files, an identifier mapping and a provenance manifest.
Both folders are ignored by Git, along with raw `data/` (about 31 GB), `outputs/`
(about 2.5 GB), `.venv/` and credentials. Local files remain available to Docker:
`.dockerignore` includes the two dashboard data folders while excluding the raw
inputs and analysis outputs from the image build context.

To prepare a data archive for a teammate, run from the repository root:

```bash
tar -czf task2-data.tar.gz task2_site/dist/data task2_site/query_data
```

Transfer that archive separately from GitHub, then extract it from the root of
the teammate's clone with `tar -xzf task2-data.tar.gz`. The archive itself is also
ignored by Git. A clone alone cannot launch the dashboard until its data is restored
or regenerated.

## A short tour for researchers

1. **Study overview:** start at 20 reads per run. Enter any positive whole number
   and select **Apply** (or press Enter); 10/20/50 remain quick presets.
   Custom cutoffs recalculate exactly and may take longer at low coverage.
   Blank, zero, negative and fractional values leave the current filter unchanged.
   A cutoff above all observed coverages gives an explicit empty result.
   The threshold dropdown affects only run summaries;
   comparisons keep the saved model threshold.
2. **Compare cell lines:** pick two cell lines or click an overlap square. Hover
   over scatter points, filter and sort candidates, and download CSV. Swap A/B
   to reverse the comparison and its score differences.
3. **Explore a transcript:** click a candidate or enter an Ensembl transcript ID.
   Select a heatmap position for individual scores, read counts, original IDs
   and replicate-weighted means. Dashes mean absent or insufficient coverage.
4. **How to read the data:** definitions, aggregation rules, limitations, provenance
   and the report. Reviewable report figures and tables also live in `reports/task2/`.

Coverage, comparison, transcript and selected position are encoded in the URL.
Share the fragment with teammates running their own copy. Candidate table filters
and the transcript cell-line display filter are local view controls and are not
encoded. CSV exports include the selected filters and model checksum.

## Reproduce with Docker (optional, compute-intensive)

The original machine already has predictions. On a fresh clone the full pipeline
downloads about 33.4 GB and runs inference. Allow at least 45 GB free disk space
and start with 8 GB or more Docker memory; full analysis combines millions of
observations and may need more. Runtime depends on CPU, memory and disk speed.
The dashboard itself uses much less.

```bash
# Full workflow using the bundled inventory snapshot.
docker compose run --rm --build analysis -m scripts.task2_pipeline

# Regenerate analysis using existing local predictions.
docker compose run --rm --build analysis -m scripts.task2_pipeline \
  --stages analyse validate report

# Validate existing results without downloading or rescoring.
docker compose run --rm --build analysis -m scripts.task2_validate

# Rebuild the dashboard image to show regenerated results.
docker compose up --build -d
```

Run stages in one container at a time. `data/`, `outputs/`, `task2_site/query_data/`
and the dashboard data and report directories are mounted, so outputs survive
container removal. The pipeline refreshes the query bundle after analysis.
Relocated inputs (including original host paths seen inside Docker) must pass a
SHA256 check against prediction metadata or their original `download.json` audit.
Keep those audit files with the input JSON when transferring datasets.

`--threads 2 --batch-size 50000` reduces inference memory use, but does not reduce
the combined analysis memory requirement. `--refresh-inventory` explicitly
discovers newer public samples; otherwise the bundled 23-run snapshot is used.
Use a separate checkout/data directory when changing models or inventories.

For a small setup check without downloading the full dataset:

```bash
docker compose run --rm --build analysis -m scripts.predict \
  --input examples/sites.jsonl --output outputs/docker_example.csv
```

Use a new output filename when repeating this check; existing files are preserved.
On Linux, `docker compose run --user "$(id -u):$(id -g)" ...` makes generated files
belong to your host user; ensure mounted output folders are writable first.
Docker Desktop handles normal Mac/Windows file sharing.

## Reproduce the analysis without Docker

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
python -m scripts.task2_bundle_queries
python -m scripts.task2_validate
python -m scripts.task2_report --publish-directory reports/task2

# Local explorer: visit http://localhost:8765
python -m scripts.task2_serve
```

Downloads resume with `--reuse-inventory --download`. They use public HTTPS,
require no AWS account, verify the inventoried ETag and byte count, and record a
SHA256 digest. The October 12, 2026 inventory contains 23 runs across seven cell
lines, totalling 33.37 GB. Allow at least 45 GB of free disk space for inputs,
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
| `task2_site/query_data/` | All saved predictions, identifier mapping and checksummed manifest for arbitrary coverage queries |
| `reports/task2/` | Compact results draft, figures, summary tables and inventory for GitHub review |

## Comparison rules

1. Score each run using identical H04 features and the saved validation
   threshold (`0.20839324867970055`). Preserve all model outputs.
2. Filter observations by reads **per run**. The primary analysis uses 20;
   users may enter any positive whole number (up to JavaScript's largest exact
   integer, 9,007,199,254,740,991). These are analysis choices, not validated
   biological cutoffs. Lower coverage admits less-supported observations.
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

The explorer retains saved comparisons for 10/20/50 and partitions ≥10-read
observations by transcript for fast lookup. Other cutoffs use exact DuckDB
queries over every saved run–site prediction, without retraining or rescoring.
Transcript requests below 10 reads (or absent from the quick-lookup catalogue)
use the complete query bundle. The static catalogue count still describes only
transcripts exported at ≥10 reads, as labelled in the overview.
Partitions are gzip-compressed and decoded on demand in the browser; compression
preserves all observations and full score precision. Use a current browser with
the standard `DecompressionStream` API for transcript lookup.
Scatter plots show up to 2,000 deterministically sampled sites and candidate
tables show up to 500 rows; their summary statistics use the full comparison.
Full pair tables for 10/20/50 are retained in the analysis outputs. Custom
comparisons calculate statistics over all shared eligible sites, then return
the top 500 candidates and a deterministic sample of up to 2,000 scatter points.
The custom scatter sample can differ from the precomputed presets.
No gene mapping or genome-coordinate
conversion is inferred from transcript identifiers.

MCF7 replicate 3 uses versioned Ensembl IDs, while replicate 4 uses unversioned
IDs. Before harmonization their exact transcript overlap was zero. After removing
the naming difference, all 969,788 shared raw coordinate pairs had identical
7-mers; 88,680 pairs had at least 20 reads in both runs. This check supports the
coordinate convention for these inputs. Unversioned identifiers still depend on
the shared reference convention and cannot prove full transcript-sequence identity.

Explorer data under `task2_site/dist/data/` and `task2_site/query_data/` is shared
separately from Git and must be restored before launching. Regenerate it
with `task2_analyse` followed by `task2_bundle_queries`, or use the pipeline.
A static-only web host supports the saved presets and ≥10-read transcript data;
custom comparisons and low-coverage transcripts require the query service.
The Docker setup runs locally; it does not publish a public URL.

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

Verified on 19 October 2026 on Apple Silicon with Docker Desktop:

- The Docker dashboard, query and analysis images built; dashboard and query
  health checks passed.
- All 23 relocated source inputs passed SHA256 verification inside Docker.
- Full bundle validation checked 28,670,080 predictions, 5,750,947 exported
  observations at ≥10 reads, and 24,502 transcripts.
- A temporary four-run fixture exercised model inference, analysis/export and
  validation inside Docker; all six example scores matched the regression reference.
- All 60 Python tests passed; 15 additional query tests cover custom
  cutoffs, sequence matching, replicate means, low coverage and HTTP validation.
  All 10 JavaScript tests passed. The latter check aggregation,
  comparison reversal, filtering, full-precision CSV, and all 63 comparison files.
- Query results at 20 reads match all 21 original pair summaries and replicate
  agreements. Docker HTTP checks verified counts at 1, 25 and 1,000,000 reads;
  the 1-read overview took about 71 seconds on this machine and retained all
  28,670,080 predictions. Results are cached while the service is running.
- Chrome desktop and 400-pixel mobile checks covered navigation, coverage changes,
  comparison reversal, candidate-to-transcript links, missing transcript feedback,
  and CSV download. Intel/AMD execution has not been tested on this machine.

```bash
# Restore or regenerate the dashboard data before tests that inspect real assets.
python -m pip install -r requirements-dev.txt
python -m pip install -r requirements-dashboard.txt
python -m pytest -q
node --check task2_site/dist/app.js
node --test tests/dashboard.test.mjs
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
