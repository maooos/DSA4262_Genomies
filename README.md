# H04 inference

Score RNA sites and predict binary labels with the saved pooled H04 model. This repository contains the
trained model and inference code. The notebook currently targets the local
`data/raw/data0/dataset0.json.gz` input. Training workflows, EDA notebooks, and
intermediate caches have been removed.

## Task 2: SG-NEx analysis and explorer

The Task 2 workflow discovers available processed SG-NEx samples, downloads them
with resume support, applies H04, compares covered sites across cell lines and
replicates, and builds an interactive explorer. See [the Task 2 guide](docs/TASK2.md)
for commands, outputs, comparison rules, and the processed-data scope.

```bash
python -m pip install -r requirements-task2.txt
python -m scripts.sgnex_download --download --workers 3
python -m scripts.task2_predict
python -m scripts.task2_analyse
python -m scripts.task2_validate
python -m scripts.task2_report --publish-directory reports/task2
python -m http.server 8765 --directory task2_site/dist
```

Open http://localhost:8765 after running the workflow. Predictions and the
results draft are saved under `outputs/task2/`; the explorer uses precomputed
data and does not retrain the model.

The completed [private SG-NEx explorer](https://sgnex-m6a-explorer.chicacorn256.chatgpt.site)
contains all 23 available processed runs. The [results draft](reports/task2/task2_results.html),
figures, inventory, and validation record are in `reports/task2/`.

## Setup

Use Python 3.12. From the repository root:

```bash
python3.12 -m venv .venv  # Skip if the existing environment is available
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The pinned numerical-library versions match the saved model. H04 uses
scikit-learn; PyTorch and XGBoost are not needed.

## Run the notebook

Open `notebooks/inference.ipynb`, select the `.venv` kernel, and edit its Settings
cell:

```python
MODEL_PATH = ROOT / "models/h04.joblib"
INPUT_PATH = ROOT / "data/new_dataset.json.gz"  # Absolute paths also work
INPUT_KIND = "reads"
```

Choose **Restart Kernel and Run All**. The notebook defaults to the local
`data/raw/data0/dataset0.json.gz` input. For a quick setup check without the full
dataset, set `INPUT_PATH = ROOT / "examples/sites.jsonl"`. Each run writes a
uniquely named CSV and companion metadata JSON to `outputs/`.

## Run from the terminal

```bash
python -m scripts.predict \
  --input data/new_dataset.json.gz \
  --output outputs/new_dataset_H04.csv
```

For the included six-site example:

```bash
python -m scripts.predict --input examples/sites.jsonl --output outputs/example.csv
```

Use `--input-kind features` for a prepared site-feature CSV or Parquet. Optional
flags include `--model`, `--threads` (default 4), and `--batch-size` (default
200,000). Run `python -m scripts.predict --help` for details. Existing output
files are preserved; choose a new filename when rerunning.

## Input format

`reads` accepts project-format JSON/JSONL (optionally gzip-compressed) or a
per-read Parquet table. JSON must have one complete site per line:

```json
{"transcript_A":{"123":{"AGAACAT":[[0.1,1.0,70.0,0.2,1.1,71.0,0.3,1.2,72.0]]}}}
```

Each read has nine finite measurements in this exact order:

1. Upstream dwell time, signal SD, signal mean.
2. Central dwell time, signal SD, signal mean.
3. Downstream dwell time, signal SD, signal mean.

Use the same units and signal preprocessing as the training inputs. Dwell times
must be positive, and sequences must contain seven A/C/G/T bases. Labels and
gene IDs are unnecessary. FAST5/POD5 files require upstream signal processing;
they are not accepted directly.

Per-read Parquet requires `transcript_id`, `transcript_position`, `sequence`,
`n_reads`, and the nine columns listed in `src.feature_engineering.RAW_FEATURE_NAMES`.
Keep all reads for each site contiguous and complete. Process separate datasets
as separate inputs so overlapping site identifiers do not collide.

`features` accepts a CSV/Parquet containing site identifiers and all 162 input
columns recorded in the model bundle. The saved pipeline selects its 153
engineered features. Raw reads automatically go through the original feature
engineering code, including log dwell times, signal summaries, contrasts, and
sequence encodings.

## Output

One row per site, with columns `transcript_id,transcript_position,score,prediction`. A pooled
feature table with a `dataset` key retains that extra identifier in its CSV.
The companion `.metadata.json` records the model checksum, input, site count,
output columns, and the threshold and rule used for prediction.

Scores are continuous class-1 model outputs, not guaranteed calibrated
probabilities. Inference does not retrain the model. The `prediction` column uses
the saved validation threshold (approximately 0.208393 for H04): `1` means
predicted modified (`score >= threshold`), and `0` means predicted unmodified.
Both the original score and the predicted label are saved. Model bundles must
include a finite threshold in [0, 1].

## Files to retain

- `models/h04.joblib`: the fitted model, feature selection, and threshold (~695 KB).
- `models/h04.metadata.json`: checksum and training provenance.
- `src/`: the inference, parsing, and feature-engineering modules.
- `notebooks/inference.ipynb` and/or `scripts/predict.py`: inference entry points.
- `requirements.txt`: runtime dependencies, including the notebook kernel.

`examples/` and `tests/` are small verification assets. `data/` is available for
future inputs and `outputs/` holds predictions; both are ignored by Git. Keep
the model and code together when transferring this project. The `.venv/`
environment can be recreated and need not be transferred.

## Verification

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

Before removing the historical data, the relocated model reproduced all 32,766
saved internal-test scores within `1e-12` tolerance. The included six raw data0
sites also reproduced their scores from the original saved features. The
regression test retains this small check without requiring the training data;
the example is not a new evaluation set.
