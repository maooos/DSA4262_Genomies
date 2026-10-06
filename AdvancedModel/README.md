# AdvancedModel

Four independent modelling notebooks share data and training utilities while producing separate models and reports.
All use **data0 + data1 + data2**. They can run in any order and do not depend on outputs from other AdvancedModel notebooks.

| Notebook | Model | Architecture |
|---|---|---|
| [01_noisy_or_mil.ipynb](01_noisy_or_mil.ipynb) | m6Anet-style Noisy-OR MIL | Read MLP → read score → Noisy-OR |
| [02_gated_attention_mil.ipynb](02_gated_attention_mil.ipynb) | Gated Attention MIL | Read MLP → gated attention → site classifier |
| [03_set_transformer_mil.ipynb](03_set_transformer_mil.ipynb) | Set Transformer MIL | Read MLP → 2 × SAB → gated pooling → site classifier |
| [04_mil_h04_ensemble.ipynb](04_mil_h04_ensemble.ipynb) | MIL + H04 ensemble | Independently trained components with blend weights selected from training OOF predictions |

The Noisy-OR notebook adapts the architecture to this project; it is not an official m6Anet reproduction.
The Set Transformer variant uses SABs with gated pooling rather than the original paper's complete PMA architecture.
No performance improvement is claimed without executing and evaluating these experiments.

## Running a notebook

From the repository root:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
```

All three `data/processed/annotated/data{0,1,2}_reads.parquet` files must exist.
If they are missing, first run `python -m scripts.prepare_datasets --download`.
Open the notebook you want to evaluate, select the `.venv` kernel, and choose **Restart Kernel and Run All**.

- `device='auto'` selects CUDA → Apple MPS → CPU. You can explicitly set `device='cpu'`.
- Defaults are at most 32 sampled reads per site, a batch size of 128, and a maximum of 25 inner training epochs.
- Each neural model requires five inner early-stopping fits, five complete fold refits, and one final fit on all outer training data.
  The ensemble notebook additionally fits five fold-specific H04 models and one final H04 model.
- If memory is insufficient, reduce `batch_size` first. Changing `max_reads` changes the model inputs and should be recorded as a new experiment.
- After an interrupted run, Run All creates a new output directory. Resuming training from an intermediate checkpoint is not currently supported.
  Completed read caches are reused.
- The existing `scripts.run_modelling` command runs the original pooled H04 notebook, not these advanced notebooks.

## Directory structure and outputs

```text
AdvancedModel/
  01_noisy_or_mil.ipynb
  02_gated_attention_mil.ipynb
  03_set_transformer_mil.ipynb
  04_mil_h04_ensemble.ipynb
  common/
    data.py          # All-read cache, groups, fitting-only normalization, bag sampling
    models.py        # Three networks, masks, site-level loss, optional auxiliary head
    training.py      # Inner early stopping, fold refits, model saving and loading
    evaluation.py    # Thresholds, metrics, plots, and grouped bootstrap intervals
    ensemble.py      # Independent H04 OOF fits, blend selection, ensemble saving and loading
```

Generated data is stored under the Git-ignored `data/processed/advanced/` directory:

```text
cache/<fingerprint>/                     # Shared by all four notebooks
runs/<architecture>/<timestamp_uuid>/    # Separate directory for each experiment
```

The cache retains every read. The nine float32 signal measurements require approximately 0.72 GB,
plus indices, site metadata, and moments. Cache creation does not fit normalization.
Each run checks source SHA-256 hashes and reuses only complete caches.
If raw-data processing changes, regenerate the annotated files; their changed contents produce a new fingerprint.

Each completed experiment includes:

- `config.json` and `split_and_fold_manifest.parquet`: configuration, versions, input/code fingerprints, and exact group membership.
- `fold_*_inner_groups.parquet`, `fold_*_history.csv`, and `training_history.csv`: early-stopping and refit records.
- `cv_fold_metrics.csv`, `cv_mean_sd.csv`, and `oof_predictions.parquet`: training CV metrics and predictions.
- `cv_oof_metrics_by_dataset.csv`: pooled OOF metrics within each dataset and the number of contributing folds.
- `validation_threshold_search.csv`: threshold selection using validation only.
- `val_metrics_by_dataset.csv` and `test_metrics_by_dataset.csv`: pooled and separate data0/data1/data2 results.
- `val_predictions.parquet`, `test_predictions.parquet`, plot PNGs, and `test_group_bootstrap_intervals.csv`.
- `model.pt`: model weights, architecture configuration, k-mer normalization parameters, read sampling configuration, threshold, and training provenance.
- `evaluation.json` and `complete.json`: evaluation results and confirmation that saving and reloading were verified.

Ensemble experiments additionally save `h04.joblib`, `ensemble.json`, `blend_weight_search.csv`, and component OOF/test results.
Their `model.pt` contains the MIL component with no standalone classification threshold;
`ensemble.json` contains the final threshold and blend weight.
Preserve the complete run directory when transferring an ensemble.

## Training and evaluation rules

1. Each observation is identified by `(dataset, transcript_id, transcript_position)`.
   The same human gene shares a group across data0/data1. All data2 mixtures from the same construct share a group;
   mixture transcript IDs are not treated as independent genes.
2. Outer split and CV policies match the original pooled notebook: seed 42 and five folds by default.
   When `REFERENCE_RUN` exists, source hashes, splits, and folds must match that run exactly.
   When its directory is absent, the notebook independently generates splits using the same policy;
   it does not automatically substitute another model or test result.
3. Assessment folds are not used for early stopping. Split each fold's fitting groups into inner fitting and stopping partitions,
   select epochs using inner stopping AP, then train from scratch on all fold-fitting sites with newly fitted normalization.
   The final model uses the median selected epoch count and is trained from scratch on outer training data.
4. Normalization statistics for log dwell, signal standard deviation, and signal mean at the three adjacent 5-mers
   come only from the current fitting reads. Unseen k-mers use that fitting partition's global statistics for the corresponding position.
   Validation and test data never update these statistics.
5. Each site contributes one classification loss. Site labels are not treated as individual read labels.
   Long bags are sampled without replacement, short bags use masks, and evaluation averages repeated samples.
6. Validation selects one shared classification threshold; test evaluation comes last.
   These datasets have been explored previously, so this is internal validation.
   CV supports development and selection; it does not provide an unbiased nested-CV estimate of generalization performance.
7. Ensemble weight α is selected only from aligned training OOF predictions from both models.
   In-sample training predictions, validation, and test results are not used to select the weight.
8. `dataset_weights` controls each dataset's contribution to training loss.
   Setting `auxiliary_weight>0` enables a mixture-fraction regression head supervised only by training data2 `label_original` values.
   This auxiliary task is disabled by default. `label_original`, dataset ID, gene ID, and construct ID are not network predictor inputs.
9. data2 has two training constructs, one validation construct, and one test construct.
   Its OOF predictions cover at most two folds, and one test construct cannot support between-construct confidence intervals.
   Reports include prior AP references and always-positive F1 values to contextualize scores under high positive prevalence.

Seeds and read sampling are reproducible; floating-point training across different PyTorch versions or hardware backends is not guaranteed to be bitwise identical.
If you change the seed, CV fold count, or outer grouping, disable the old `REFERENCE_RUN` comparison and record the change as a distinct experiment.

## Loading a model

Neural models use PyTorch checkpoints and cannot be loaded by the original joblib-based `src.inference.load_model_bundle` function.
Once the notebook has prepared `store`, use:

```python
from AdvancedModel.common.training import load_mil_bundle, predict_bags

model, normalizer, config, bundle = load_mil_bundle(run_dir / "model.pt", device="cpu")
scores = predict_bags(model, normalizer, store, store.ids("test"), config)
```

For ensembles, use `AdvancedModel.common.ensemble.predict_saved_ensemble` with read bags and tabular features in matching order.
These utilities support the current three annotated datasets. New raw inputs must first be adapted to the same bag and feature formats.

## Verification

```bash
python -m pytest tests/test_advanced_models.py -q
```

Tests use small synthetic fixtures to check padding, permutation invariance, Noisy-OR gradients, training-partition boundaries,
fold-specific normalization, OOF alignment, complete CV/refit/evaluation workflows, checkpoint and ensemble reloading, and notebook structure.

## References

- [m6Anet: multiple instance learning for m6A](https://www.nature.com/articles/s41592-022-01666-1)
- [Attention-based Deep Multiple Instance Learning](https://proceedings.mlr.press/v80/ilse18a.html)
- [Set Transformer](https://proceedings.mlr.press/v97/lee19d/lee19d.pdf)
- [PyTorch MPS backend](https://docs.pytorch.org/docs/stable/notes/mps.html)
