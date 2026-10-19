# SG-NEx prediction on Ronin

This pipeline reuses the bundled H04 model and existing inference CLI. It does
not train a model. Use one sample first, then expand the manifest.

## Prepare the instance

Use your Ronin SSH connection and storage path. Transfer this branch's code,
`models/`, and `examples/`; do not transfer a local virtual environment.
From the repository root on Ronin, with Python 3.12 available:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m pytest -q
python -m scripts.predict --input examples/sites.jsonl --output outputs/smoke/predictions.csv
```

These commands install dependencies and run tests/prediction when you execute
them. Pins must be available and compatible with the saved model; investigate
installation failures rather than silently upgrading scikit-learn.

## Download a sample directly onto Ronin

AWS CLI must be available. Use the official
[access tutorial](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/AWS_data_access_tutorial.md),
[processed sample list](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/samples_with_RNAmod_data.tsv),
and [sample metadata](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/samples.tsv).
Select direct RNA samples with processed m6Anet inputs. Processed files are
available for a subset of SG-NEx samples; record any unavailable samples.

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/processed_data/m6Anet/
sample_name="REPLACE_WITH_LISTED_SAMPLE"
aws s3 ls --no-sign-request "s3://sg-nex-data/data/processed_data/m6Anet/${sample_name}/" --recursive --human-readable --summarize
df -h .
aws s3 sync --no-sign-request "s3://sg-nex-data/data/processed_data/m6Anet/${sample_name}/" "data/raw/sgnex/${sample_name}/"
```

Inspect the downloaded filename and use its signal JSON in the manifest.
Metadata/index files are not labels and are not passed to this inference code.
Keep the original downloads. No FAST5/BAM processing is needed for compatible
processed JSON inputs. Measure one-sample RAM/runtime before choosing resource
limits or downloading the full set.

## Create a manifest and run

Before downloading data, verify the runner with the included example manifest:

```bash
python3 -B -m scripts.run_pipeline --manifest examples/pipeline_samples.csv --output-dir outputs/example_run --dry-run
```

This only checks paths and prints the plan; it runs no inference. The
`data/sgnex_samples.csv` below is a file you must create after downloading your
samples, not a bundled file.

Create `data/sgnex_samples.csv` with exactly this header. Relative input paths
are resolved from the manifest's directory; replace the example filenames with
actual downloaded files. Sample IDs must be unique letters, digits, `_` or `-`.

```csv
sample_id,input_path
sample1,raw/sgnex/sample1/data.json
sample2,raw/sgnex/sample2/data.json
```

```bash
python -m scripts.run_pipeline --manifest data/sgnex_samples.csv --output-dir outputs/sgnex_run1 --dry-run
python -m scripts.run_pipeline --manifest data/sgnex_samples.csv --output-dir outputs/sgnex_run1 --threads 4 --batch-size 200000
python -m scripts.run_pipeline --manifest data/sgnex_samples.csv --output-dir outputs/sgnex_run1 --resume
# Retry recorded failures only:
python -m scripts.run_pipeline --manifest data/sgnex_samples.csv --output-dir outputs/sgnex_run1 --retry-failed
```

Dry run checks paths only, writes nothing, and does not load the model or verify
signal compatibility. Actual runs process samples sequentially using the active
Python interpreter. Each sample directory contains:

- `predictions.csv`: existing score and thresholded prediction output.
- `predictions.metadata.json`: existing inference provenance.
- `submission.csv`: validated `transcript_id,transcript_position,score` only.
- `run.log`: prediction process output.
- `complete.json`: input/model/code/output hashes and site count, or
  `failed.json`: failure details.

The runner takes an exclusive output-directory lock. Each sample is built in a hidden
work directory and published only after prediction and submission validation succeed.
Completion records are written atomically. Resume skips only checksum-verified completed
samples. Hashing large inputs adds disk I/O. Failed sample work is preserved with a
traceback and failed.json; pass --retry-failed to rerun only recorded failures. Unknown
or changed existing sample folders are preserved and cause a refusal. Per-sample and
per-run JSON records capture timestamps, elapsed time, resource settings, Python/platform
details, and the Git revision when available.

Keep sample identity when combining outputs for cell-line comparisons. Scores
are model outputs, not guaranteed calibrated probabilities. Submission exports
validate identifiers and scores but do not independently compare coverage to
raw JSON; use the existing inference tests and small example for integration
verification. Record the code revision, environment, download source, sample
metadata, and resource measurements alongside results. Retrieve outputs and
shut down unused Ronin compute through your normal project procedure.

To export a single existing prediction CSV separately:

```bash
python -m scripts.export_submission --input outputs/predictions.csv --output outputs/submission.csv
```

Model and inference dependencies are unnecessary for exporter/runner unit tests:

```bash
python3 -B -m unittest tests.test_pipeline -v
```
