# DSA4262-Genomies

## Repository contents

- `src/data_parser.py`: converts signal JSON files into per-read Parquet tables.
- `scripts/parse_data.py`: command-line interface for parsing.
- `scripts/prepare_datasets.py`: complete raw -> parsed -> gene-annotated dataset workflow.
- `scripts/run_modelling.py`: command-line execution of the pooled modelling notebook, with an optional preparation-only stop.
- `src/modelling_data.py`: pooled site identities, shared gene/construct splits and batched features.
- `AdvancedModel/`: four independent read-level MIL / ensemble notebooks with shared training utilities.
- `scripts/annotate_gene_ids.py`: attach recovered annotations to all rows of a parsed Parquet.
- `scripts/inspect_data.py`: inspects the first signal site and its label in data0.
- `notebooks/EDA_Findings.ipynb`: consolidated findings and next steps.
- `notebooks/inference.ipynb`: generate per-site scores with a saved model; defaults to H04 and its internal test split.
- `notebooks/indiv_data_eda/`: detailed EDA code and saved outputs for each dataset.
- `handout_project2_RNAModifications.html`: project requirements.

## Setup

Using Python 3.12:

```bash
git clone https://github.com/maooos/DSA4262_Genomies.git
cd DSA4262_Genomies

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Select `.venv` as the Python environment when opening notebooks in VS Code.

## Start here: raw data to modelling

Run these commands **from the repository root**. Set up the environment once:

```bash
python3.12 -m venv .venv  # Skip if .venv already exists
source .venv/bin/activate
python -m pip install -r requirements.txt
# macOS only, if the OpenMP runtime is not installed:
brew install libomp
```

**From raw inputs** (all six teacher-provided files in the layout below):

```bash
python -m scripts.prepare_datasets --download
```

This parses and annotates data0/data1/data2. `--download` allows public reference
retrieval; it does not download the raw inputs. Existing complete stages are
verified and reused. With cached references, omit `--download` for offline use.

**If annotation is already complete**, skip that command and check the inputs:

```bash
ls -lh data/processed/annotated/data{0,1,2}_reads.parquet
```

**Train one model on all three datasets**, including site feature preparation,
shared gene/construct splitting, grouped CV, model selection and evaluation:

```bash
python -m scripts.run_modelling
```

Alternatively open `notebooks/data0_modelling_pipeline.ipynb`, select `.venv`,
and choose **Restart Kernel and Run All**. The filename is retained, but it now
loads **data0 + data1 + data2**. Both entry points execute the same notebook code.
The CLI uses the active Python environment and saves plot PNGs; it does not
write executed outputs back into the source notebook.

Optional: **only prepare model inputs, without training**:

```bash
python -m scripts.run_modelling --prepare-only
```

This saves `X_train/val/test.parquet`, `y_train/val/test.parquet`, aligned
metadata and `split_manifest.parquet`, then stops before CV/model fitting.
The full command already includes this preparation, so you do not need both.
A later full run creates a new run directory and recomputes preparation.
Use `--threads 2` to limit model threads, or
`--processed-dir data/processed/rebuild_1` if preparation used another directory
(the notebook equivalent is changing `PROCESSED_DIR` in its configuration cell).

The full default plan is **21 configurations × 5 CV folds**, followed by final
refits. Its output directory is printed and has this form:

```text
data/processed/modelling/pooled_data0_data1_data2_group_split_seed42/person_c_runs/<run_id>/
```

Read `test_metrics_by_dataset.csv` for the separate data0/data1/data2 test
results (`selected_by_cv=True` identifies the chosen model). Other key files are
`split_dataset_summary.csv`, `cv_oof_metrics_by_dataset.csv`,
`validation_metrics_by_dataset.csv`, `cv_leaderboard.csv`,
`test_predictions.parquet`, and `selected_model.joblib`.

For subsequent scoring, set `RUN_DIR` in `notebooks/inference.ipynb` to that
completed run. `DATASET="internal_test"` scores the pooled held-out observations.
Choosing an entire data1/data2 input also scores its training observations;
those full-input scores are not a held-out evaluation of the pooled model.

## Advanced read-level models

[AdvancedModel/README.md](AdvancedModel/README.md) documents four independent
notebooks, each run with the `.venv` kernel and **Restart Kernel and Run All**:

| Notebook | Architecture |
|---|---|
| [01_noisy_or_mil.ipynb](AdvancedModel/01_noisy_or_mil.ipynb) | m6Anet-style read MLP and Noisy-OR pooling |
| [02_gated_attention_mil.ipynb](AdvancedModel/02_gated_attention_mil.ipynb) | Gated Attention MIL |
| [03_set_transformer_mil.ipynb](AdvancedModel/03_set_transformer_mil.ipynb) | Two set-attention blocks and gated pooling |
| [04_mil_h04_ensemble.ipynb](AdvancedModel/04_mil_h04_ensemble.ipynb) | Independently fitted MIL + H04, blended using training OOF |

Install the updated `requirements.txt` for PyTorch. These notebooks share an
all-read disk cache, use the same human gene / synthetic construct splits, and
fit k-mer normalization within each fitting partition. Early stopping uses
inner training groups, then each CV model is refitted on its complete fitting
partition. Validation selects the threshold and test is evaluated last.
The ensemble notebook trains its own components; it does not require the other
three notebooks to have run. No full advanced-model training results are claimed
by this implementation; only small synthetic tests have been executed.

Outputs go to `data/processed/advanced/runs/<architecture>/<run_id>/`.
The existing `scripts.run_modelling` command still runs the original pooled
notebook, not these advanced notebooks.

## Generate prediction scores

Open `notebooks/inference.ipynb`, edit its **Settings** cell, and choose **Restart
Kernel and Run All**. `RUN_DIR` selects a completed modelling run; `MODEL_PATH`
selects its `selected_model.joblib` saved by `data0_modelling_pipeline.ipynb`.
`DATASET` selects `internal_test`, `data1`,
`data2`, or `custom`. Custom inputs can be raw project JSON/JSONL (including
gzip), parsed read Parquet, or compatible site-feature CSV/Parquet tables.
The separate `newadjustmentfolder` model bundles use different preprocessing
and are not compatible with this notebook's loader.

The notebook defaults to H04 from run `20261005T103445_98fad904`, if that local
run is available, and generates scores for its 18,386 internal test sites. It
also verifies those scores against the saved historical predictions. For a new
training run, set `RUN_DIR` to your own output directory. It loads the fitted
model without retraining.

Each run writes a new CSV under `data/processed/inference/` with columns
`transcript_id,transcript_position,score` (with a leading `dataset` column for
pooled feature inputs), plus a `.metadata.json` file recording
the model and input. The notebook prints both paths and previews the scores.

### Intermediate submission: all three datasets with H04

From the repository root, use the completed pooled H04 run:

```bash
source .venv/bin/activate
python -m scripts.generate_submissions \
  --run-dir data/processed/modelling/pooled_data0_data1_data2_group_split_seed42/person_c_runs/20261006T085253_d868f3c6
```

This loads one fitted H04 model and scores every site using the union of its
saved training, validation and test features. It does not retrain or threshold
the probabilities. Original metadata is used only for identifiers and row order.
The command verifies complete site coverage, unique identifiers, finite scores
in `[0, 1]`, and the exact header `transcript_id,transcript_position,score`.

Upload these **three CSV files** from `data/processed/submissions/`:

- `Genomies_dataset0.csv` — 121,838 sites
- `Genomies_dataset1.csv` — 90,810 sites
- `Genomies_dataset2.csv` — 1,323 sites

`Genomies_submission_metadata.json` records model/input/output checksums for
reproducibility; it is not a submission file. Existing outputs are preserved;
use `--output-dir` to generate a new copy. Full-dataset predictions include
training observations, so their metrics are not a held-out evaluation.

## Local datasets

Raw datasets and generated Parquet files are excluded from Git because of their size. Obtain the original datasets separately and arrange them as follows:

```text
data/
├── raw/
│   ├── data0/
│   │   ├── dataset0.json.gz
│   │   └── data.info.labelled
│   ├── data1/
│   │   ├── dataset1.json.gz
│   │   └── data.info
│   └── data2/
│       ├── dataset2.json.gz
│       └── data.info
└── processed/              # Generated by the parsing commands below
```

## Generate the Parquet files

### Complete preparation workflow (recommended)

Run from the repository root, with the virtual environment activated:

```bash
python -m scripts.prepare_datasets --download
```

`--download` allows missing public reference files to be downloaded. With the
references/mappings already cached, `python -m scripts.prepare_datasets` works
offline. The workflow performs these stages in order:

1. Invoke `scripts.parse_data` for data0, data1 and data2 using their original
   raw JSON and metadata files. Write complete per-read Parquets to `parsed/`.
2. Validate the cached recovery mappings against the raw-file and mapping
   checksums, or run `scripts.recover_gene_ids` to create them.
3. Join the recovered annotations by `(transcript_id, transcript_position)`
   onto **every read**, writing the complete final datasets to `annotated/`.

```text
data/raw/dataN/datasetN.json.gz + data.info[.labelled]
    -> scripts.parse_data
data/processed/parsed/dataN_reads.parquet
    -> scripts.annotate_gene_ids + gene_id_recovery/*.csv
data/processed/annotated/dataN_reads.parquet   # Final dataset to load
```

| Final dataset | Read rows | Sites | Gene annotation |
|---|---:|---:|---|
| `data/processed/annotated/data0_reads.parquet` | 11,027,106 | 121,838 | Existing gene IDs verified and retained |
| `data/processed/annotated/data1_reads.parquet` | 7,907,952 | 90,810 | All missing gene IDs filled |
| `data/processed/annotated/data2_reads.parquet` | 1,171,940 | 1,323 | Synthetic reference IDs, positions and evidence status added |

These final files contain all nine signal features, sequences, coverage and read
indices from the parsed datasets, plus the annotations. data0/data1 labels are
unchanged. For data2, the final `label` is an integer: **0 if the original label
is 0, and 1 if the original label is greater than 0**. The original mixture
proportion is preserved in `label_original`; raw, parsed and recovery mapping
files keep their original labels. This binary target means a positive mixture
proportion at the site, not that every read is modified. These are
full datasets, not mapping-only tables. Each has an adjacent `.audit.json`
recording source/output checksums and validation counts. The annotation stage
verifies every site's read count and compares all original non-gene columns
against the written Parquet, including row order (using `label_original` for
data2). It also verifies the binary conversion and records `label_transformation`
and `original_label_column` in the audit. It rejects conflicting gene
IDs, missing mappings, duplicate mappings and inconsistent metadata.

The final data2 Parquet uses `gene_id="Synthetic Curlcake RNA: {geo_reference_id}"`
as a readable reference-specific description, for example
`Synthetic Curlcake RNA: cc6m_2244_t7_ecorv`, with
`gene_id_status="not_applicable_synthetic_construct"`.
There are four distinct descriptions, one for each mapped synthetic construct;
these are not biological gene IDs. Each of the seven `tx_id_*` mixture groups
contains sites from all four constructs, so a transcript ID does not determine
one unique construct. Assignments retain their evidence strength in
`reference_mapping_status`, including order-dependent inferences.
The audit records `genes=0` and `gene_id_unique_values=4`.
Raw, parsed and recovery mapping files retain their original null gene IDs.
Its `geo_reference_id`
and `reference_mapping_status` columns carry the recovered synthetic source
and its evidence strength, as explained below. These are not human gene IDs.

Successful stages are reused only when their checksums still match. To rebuild
after input changes, choose a new destination, for example
`--processed-dir data/processed/rebuild_1`. Raw files and earlier datasets are
preserved. The older `data/processed/dataN_reads.parquet` files are legacy
parsed outputs; load the files in **`annotated/`** for complete annotations.
The modelling notebook and the inference notebook's data1/data2 choices use
these annotated files by default. Older EDA notebooks may still use the legacy
path; point their data-loading setting to the corresponding final file.

```python
import pandas as pd
from src.feature_engineering import build_features

reads = pd.read_parquet("data/processed/annotated/data1_reads.parquet")
X = build_features(reads)  # Existing site-level feature engineering
```

The two stages can also be run separately (example for data1):

```bash
python -m scripts.parse_data \
  --signals data/raw/data1/dataset1.json.gz \
  --labels data/raw/data1/data.info \
  --output data/processed/parsed/data1_reads.parquet

# If recovery mappings have not been generated yet:
python -m scripts.recover_gene_ids --download

python -m scripts.annotate_gene_ids \
  --dataset data1 \
  --input data/processed/parsed/data1_reads.parquet \
  --mapping data/processed/gene_id_recovery/data1.gene_ids.csv \
  --output data/processed/annotated/data1_reads.parquet
```

The individual stage commands refuse existing output files; use
`scripts.prepare_datasets` to resume an existing preparation safely.

### Recover missing gene IDs

`scripts/recover_gene_ids.py` recovers **all 90,810 data1 sites** (4,451
transcripts, 3,149 gene IDs) using the full Ensembl release 91 annotation,
including patches and alternate haplotypes. The required file is
`Homo_sapiens.GRCh38.91.chr_patch_hapl_scaff.gtf.gz`; the smaller default GTF
omits 145 data1 transcripts. This reference also reproduces every original
data0 gene ID, with zero conflicts across 121,838 sites. Each site's seven-base
sequence matches the release 91 transcript FASTA at its supplied position in
both datasets, and the FASTA gene headers agree with the GTF.

The [SG-NEx annotation documentation](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/ANNOTATIONS.md)
identifies release 91. The script downloads the complete annotation and the
coding/noncoding transcript FASTAs from the official Ensembl archive and checks
their SHA-256 checksums before use.

```bash
python -m scripts.recover_gene_ids --download
# Subsequent runs use the cached references:
python -m scripts.recover_gene_ids
```

Outputs are under `data/processed/gene_id_recovery/`:

- `data1.gene_ids.csv`: all original data1 metadata plus `gene_id`, `gene_name`,
  sequence verification and reference provenance. This can be passed as
  `--labels` to `scripts.parse_data` when generating a **new** annotated Parquet.
- `data0.gene_ids.csv`: original data0 labels with their reference cross-checks.
- `human_transcript_to_gene.csv`: mapping for the union of data0/data1 transcripts.
- `data2.reference_mapping.csv`: synthetic reference assignments, with blank
  `gene_id` and an explicit `not_applicable_synthetic_construct` status.
- `data2.reference_candidates.csv`: every seven-base match in the Curlcake
  reference, retaining ambiguities instead of silently discarding alternatives.
- `data2.mixture_signal_evidence.csv`: exact endpoint-signal matching counts.
- `all_sites.gene_ids.csv`: all 213,971 sites with dataset and mapping status.
- `recovery.audit.json`: counts, checksums, source URLs and interpretation limits.

The script writes derived outputs; it does not replace raw data, existing read
Parquets or trained models. Rerunning replaces its own generated CSV/audit files.

**data2 is strongly supported as a synthetic Curlcake mixture dataset.** Its
189 ordered sequence contexts form a unique ordered subsequence of 192 DRACH
sites in the four [GEO GSE124309 Curlcake reference sequences](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124309).
The [original Curlcake study](https://www.nature.com/articles/s41467-019-11713-9)
describes computationally designed sequences covering all possible five-base
contexts; these have reference names rather than human Ensembl gene IDs.
All seven `tx_id_*` groups
share the same contexts. Exact comparison of position, sequence and all nine
numeric signal features assigns every one of the 1,171,940 read rows to the
`tx_id_0` or `tx_id_6` endpoints; the intermediate groups have approximately
95%, 75%, 70%, 50% and 25% reads from `tx_id_0`, matching their labels.
A [public copy of dataset2](https://github.com/bce99/m6A-RNA-Modification-Prediction/blob/main/dataset2.json.gz)
also has the same decompressed SHA-256 as our file.

The four inferred source constructs contribute 36, 47, 49 and 57 sites per
mixture, respectively: `cc6m_2244_t7_ecorv`, `cc6m_2459_t7_ecorv`,
`cc6m_2595_t7_ecorv` and `cc6m_2709_t7_ecorv`. Of the 189 contexts, 97 have a
unique seven-base match **within the Curlcake reference**. The remaining 92
assignments require the explicit assumption that anonymization preserved the
order of sites across the four references. The script checks that the complete
ordered alignment is unique, but this is still inference: the original renaming
key and dataset-generation code are unavailable.

Both the GEO and [EpiNano reference FASTAs](https://github.com/novoalab/EpiNano/blob/master/Reference_sequences/cc.fasta)
are retained. Their complete construct cores agree, but EpiNano includes a
31-base prefix and a three-base suffix. Output coordinates are explicitly
zero-based central-base positions in each named FASTA; they are not claimed to
be the unknown preprocessing pipeline's original coordinates. A synthetic
construct identifier or display description does not establish a biological
gene identity.

### Parse signal data

Run these commands from the repository root:

```bash
python -m scripts.parse_data \
  --signals data/raw/data0/dataset0.json.gz \
  --labels data/raw/data0/data.info.labelled \
  --output data/processed/parsed/data0_reads.parquet

python -m scripts.parse_data \
  --signals data/raw/data1/dataset1.json.gz \
  --labels data/raw/data1/data.info \
  --output data/processed/parsed/data1_reads.parquet

python -m scripts.parse_data \
  --signals data/raw/data2/dataset2.json.gz \
  --labels data/raw/data2/data.info \
  --output data/processed/parsed/data2_reads.parquet
```

Each command creates a Parquet table and a corresponding `dataN_reads.audit.json` file containing parsing counts and label-matching information. The output directory is created automatically.

Each Parquet row represents one read observation at a transcript-position site, with nine signal features, sequence information, coverage and available annotations. Labels are joined using transcript ID and position.

Existing output files are not overwritten. For a smaller trial/sample, add `--max-sites 100` and use a different output filename.

## EDA

The individual notebooks read from `data/processed/`. Set the notebook working directory to the repository root when executing them.

- `eda0.ipynb`: EDA for data0.
- `eda1.ipynb`: EDA for data1.
- `eda2.ipynb`: EDA for data2.
- `EDA_Findings.ipynb`: summary for feature engineering and modelling.

## Feature engineering

`src/feature_engineering.py` turns a per-read Parquet table (the output of
`src/data_parser.py`) into one feature row per `(transcript_id,
transcript_position)` site.

```python
import pandas as pd
from src.feature_engineering import build_features, FEATURE_COLUMNS

df = pd.read_parquet("data/processed/annotated/data0_reads.parquet")
X = build_features(df)          # one row per site, no label/gene_id
X[FEATURE_COLUMNS]              # feature matrix only, keys dropped
```

Required input columns: `transcript_id`, `transcript_position`,
`sequence`, `n_reads`, and the nine raw signal columns. `gene_id` and
`label` are ignored if present, so the same function works on labelled
training data and on label-free new data. The output always has the same
columns (`FEATURE_COLUMNS`, 153 columns) and no missing values,
regardless of how many reads or sites are in the input.

### Features

- **Dwell time**: log-transformed, since dwell times are strictly
  positive and right-skewed.
- **Center-flank contrasts**: for dwell, signal sd, and signal mean, the
  central position's value minus each flank's value
  (`{dwell,sd,mean}_contrast_{minus1,plus1}`).
- **Site-level aggregation**: the 15 per-read features above are
  summarised per site with mean, sd, median, p25, p75, min, and max (105
  columns). A site with a single read has its `_sd` columns set to 0
  rather than left undefined.
- **Coverage**: `n_reads`.
- **Sequence**: one-hot encoding per position of the 7-mer
  (`seq_pos{0-6}_{A,C,G,T}`, 28 columns).
- **Central motif** (`motif_<5-mer>` / `motif_other`, 19 columns): one-hot
  encoding of the central 5-mer as a single DRACH motif category. The
  central motif's identity alone is strongly associated with the
  modification label (see `notebooks/EDA_Findings.ipynb`, "Main
  findings"); this makes that signal directly available, rather than
  only recoverable from a combination of the position-wise one-hot
  columns above.

### Optional: k-mer-normalised signal

Nanopore current and dwell readings depend heavily on which k-mer is
passing through the pore, largely independent of modification state. For
this reason, m6Anet's own model uses signal intensity, standard
deviation, and dwell time normalised per k-mer rather than raw values
(Hendra et al., 2022, *Nature Methods*, ["Detection of m6A from direct
RNA sequencing using a multiple instance learning
framework"](https://www.nature.com/articles/s41592-022-01666-1)). Our
nine raw features are not normalised this way, so `fit_kmer_baselines`
computes each central motif's typical signal level, and
`build_features(df, kmer_baselines=...)` uses it to add three residual
features (`central_{log_dwell,signal_sd,signal_mean}_kmer_resid`),
aggregated the same way as the other features (21 extra columns,
`FEATURE_COLUMNS_WITH_KMER_RESID`).

This is opt-in rather than always on: `fit_kmer_baselines` must be fit on
the training split only, then reused — never refit — on validation,
test, and new data, so that no split's baseline is computed from its own
signal.

```python
from src.feature_engineering import build_features, fit_kmer_baselines

kmer_baselines = fit_kmer_baselines(train_df)
X_train = build_features(train_df, kmer_baselines=kmer_baselines)
X_test = build_features(test_df, kmer_baselines=kmer_baselines)
```

Run the tests with:

```bash
python -m pytest tests/
```

## Modelling and evaluation

The complete workflow is in
[`notebooks/data0_modelling_pipeline.ipynb`](notebooks/data0_modelling_pipeline.ipynb):
data loading, training-set EDA, feature engineering, grouped cross-validation,
model comparison, hyperparameter tuning, threshold selection and evaluation.

### Run the notebook or command

Use `python -m scripts.run_modelling`, or open the notebook with the `.venv`
kernel and **Run All**. See the exact setup/preparation commands above.
It loads all three annotated Parquets and fails early if one is missing or
has incompatible annotations. No raw parsing happens inside modelling.

### Data and evaluation design

The model uses **data0 + data1 + data2 together**. One observation is one
`(dataset, transcript_id, transcript_position)` site, aggregated over its reads.
Shared human sites remain separate observations rather than having their
signals averaged across datasets. The target is binary `label`; data2 retains
`label_original` and uses `label_original > 0`, meaning positive mixture at the
site rather than every read being modified.

`group_id` uses `human:<gene_id>` across both data0 and data1, and
`synthetic:<geo_reference_id>` for data2. The same human gene cannot cross
partitions or CV folds, even if observed in both datasets. For synthetic data,
all mixtures of a construct stay together, including reused endpoint reads.
`tx_id_*` identifies mixtures spanning multiple constructs, so it is not a
suitable gene split key. Construct assignments include inferred mapping
evidence; the split is conditional on that recovered mapping.

The outer split sorts and shuffles the union of human genes with `seed=42`,
allocating approximately 70% / 15% / 15% of genes to train / validation / test.
Synthetic constructs are shuffled separately with the same seed: four groups
allow **2 train / 1 validation / 1 test**, not exact 70/15/15 fractions.
Each dataset/split must contain both classes. Actual counts and rates are saved
in `split_dataset_summary.csv`; no label-based seed search is performed.
The historical data0-only partition membership does not apply to this new run.

The shared `X_train`, `X_val` and `X_test` tables contain 162 columns: the
153 existing engineered features plus nine arithmetic means of the raw
signals. Each model pipeline selects its own subset: nine raw means for the
mean-only baseline, 47 sequence/motif features for the sequence-only model,
or 153 engineered features for the full models. `y_train`, `y_val` and
`y_test` contain aligned labels. Identifiers remain in the index or metadata,
not in predictor columns. Optional k-mer residuals are disabled.

All candidates use the same **five-fold `GroupKFold` by `group_id`, within
the training partition only**. Scaling and class weights are fitted within
each fold's training portion. The report includes ROC-AUC, trapezoidal PR-AUC
and average precision (AP), each as **mean ± sample standard deviation
(`ddof=1`) across folds**. AP and trapezoidal PR-AUC are reported separately.

Only two training constructs are available for data2, so at most two CV
assessment folds include it. `cv_oof_metrics_by_dataset.csv` reports metrics on
each dataset's combined OOF predictions and its actual contributing fold count;
it does not claim five independent data2 folds. Pooled CV weights sites equally,
so larger datasets contribute more. Dataset identity is not a predictor.

Model selection maximises mean pooled CV AP, with ROC-AUC and then experiment ID
as tie-breakers. The selected configuration is fitted on all training data;
validation selects the threshold that maximises F1. The model is not refitted
on train + validation afterward. Test results do not determine model choice,
hyperparameters or the threshold in this workflow.

### Models, imbalance handling and tuning

The default experiment plan contains **21 configurations across four model
families**, giving **105 CV fits**:

| Model family | Configurations | Comparisons |
|---|---:|---|
| DummyClassifier | 1 | Constant training-prior reference |
| Logistic Regression | 2 | Sequence-only versus full engineered features |
| HistGradientBoosting | 10 | Raw means versus full features, class weighting, six search candidates |
| XGBoost | 8 | Positive-class weighting and six search candidates |

E01 is the required **unweighted gradient-boosted tree baseline on nine raw
mean features**. Paired experiments compare unweighted and balanced HGB,
and XGBoost with `scale_pos_weight=1` versus the fold-training
negative/positive ratio. The XGBoost search also considers the square root
of that ratio. Weights are calculated without assessment-fold labels.

Reproducible random searches vary learning rate, boosting iterations,
tree complexity, regularisation and weighting; XGBoost additionally varies
row and column sampling. Six configurations per boosting family are sampled
before inspecting scores. This is a bounded search, not an exhaustive one.

### Historical data0-only results (not pooled results)

The new pooled pipeline has not been run as part of this code change. The figures
below describe an earlier **data0-only** model and do not predict the new result.

The completed run `20261004T183731_5bbdbfb1` selected **H04:
HistGradientBoosting with all 153 engineered features and no class weights**.
Its settings are `learning_rate=0.03`, `max_iter=350`, `max_leaf_nodes=31`,
`min_samples_leaf=20`, `l2_regularization=2.0`, `early_stopping=False` and
`random_state=42`.

| Metric | Training CV: mean ± sd | Historical test |
|---|---:|---:|
| ROC-AUC | 0.9183 ± 0.0104 | 0.9279 |
| PR-AUC (trapezoidal) | 0.4825 ± 0.0320 | 0.4967 |
| Average precision | 0.4834 ± 0.0317 | 0.4973 |

The validation-selected threshold is **0.186744**. On the test partition,
precision is **0.4721**, recall **0.5997**, and F1 **0.5283**: 466 true
positives, 521 false positives, 311 false negatives and 17,088 true negatives.
The notebook includes PR/ROC curves, a confusion matrix and 95% percentile
intervals from 500 bootstrap samples of whole test genes.

Mean CV AP improves from **0.3991** for the raw-mean baseline to **0.4834**
for H04. The best XGBoost configuration, X05, is close at **0.4815**;
the scores do not establish a clear superiority of HGB over XGBoost.

CV scores were used for tuning and selection, so they are not unbiased
nested-CV estimates. The test partition was inspected in earlier work and
is therefore reported as a **historical internal evaluation**, not a newly
untouched test set. Gene-bootstrap intervals do not capture training or
selection uncertainty. External validation and the separate m6Anet benchmark
remain necessary for broader performance claims.

### Experiment logs and saved outputs

With `SAVE_OUTPUTS=True` (the default), each run writes to a new directory:

```text
data/processed/modelling/pooled_data0_data1_data2_group_split_seed42/person_c_runs/<run_id>/
```

Key outputs include:

- `experiment_plan.csv`, `experiment_log.csv`, and `fold_results.csv`:
  experiment IDs, changes from parent experiments, parameters, realised
  fold-training weights, timestamps, fold scores and aggregate metrics.
- `cv_fold_manifest.parquet`, `cv_leaderboard.csv`, `cv_mean_sd_report.csv`,
  and `paired_experiment_changes.csv`: fold membership and model comparisons.
- `oof_<experiment_id>.parquet`: out-of-fold predictions for each candidate.
- `split_dataset_summary.csv`, `cv_oof_metrics_by_dataset.csv`,
  `validation_metrics_by_dataset.csv`, `test_metrics_by_dataset.csv`: separate
  dataset counts and evaluation, using one shared validation-selected threshold.
- `test_group_bootstrap_intervals.csv`: pooled intervals stratified by group kind;
  the single held-out synthetic construct is fixed, so these do not quantify
  between-construct uncertainty. No data2-specific interval is claimed.
- `X_*.parquet`, `y_*.parquet`, and `metadata_*.parquet`: prepared inputs,
  labels and site metadata with aligned indices.
- `selected_model.joblib` and `recommendation.json`: the fitted pipeline,
  feature schema, threshold, selected configuration and recommendation.
- `validation_*.csv`, `test_*.csv`, and `test_predictions.parquet`:
  evaluation metrics, threshold results, bootstrap intervals and predictions.
- `cv_started.json` and `cv_complete.json`: run provenance, package versions,
  input/code checksums and CV completion details.

A cumulative `person_c_experiment_history.csv` in the parent
`pooled_data0_data1_data2_group_split_seed42/` directory preserves completed experiments across
runs. These generated files are under the Git-ignored `data/` directory;
the source modelling notebook has cleared outputs to avoid displaying stale data0-only results.
The tests in `tests/test_modelling_notebook.py`, `tests/test_modelling_data.py`
and `tests/test_inference.py` check metrics, cross-dataset gene separation,
construct/mixture separation, feature batching, inference identity, fold-local
scaling and weighting, and append-only logging using small synthetic fixtures.
