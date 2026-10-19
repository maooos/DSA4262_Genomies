"""Bundle all scored observations for exact, arbitrary dashboard coverage filters."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from src.task2_analysis import harmonize_transcripts


def bundle(predictions, site_data, output):
    study = json.loads((site_data / "summary.json").read_text())
    output.mkdir(parents=True, exist_ok=True)
    identifiers, samples = [], []
    for sample in study["samples"]:
        source = predictions / f"{sample['sample_id']}.parquet"
        meta = json.loads(source.with_suffix(".metadata.json").read_text())
        count = pq.read_metadata(source).num_rows
        if (meta["model_sha256"] != study["model_sha256"] or
                meta["prediction_threshold"] != study["threshold"] or count != meta["site_count"]):
            raise ValueError(f"Prediction provenance mismatch: {source}")
        identifiers.extend(pq.read_table(source, columns=["transcript_id"])
                           .column(0).unique().to_pylist())
        target = output / source.name
        temporary = target.with_suffix(".partial")
        shutil.copyfile(source, temporary)
        with temporary.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        temporary.replace(target)
        samples.append({**sample, "file": target.name, "site_count": count, "sha256": digest})
        print(f"Bundled {sample['sample_id']}: {count:,} observations", flush=True)
    # Resolve identifiers across ALL coverages, preserving conflicting revisions.
    ids = pd.DataFrame({"transcript_id": sorted(set(identifiers))})
    ids["sample_id"] = ids.transcript_id
    ids["transcript_position"] = 0
    resolved, audit = harmonize_transcripts(ids)
    mapping = resolved.rename(columns={"transcript_id_original": "original", "transcript_id": "canonical"})
    mapping[["original", "canonical"]].to_parquet(output / "identifiers.parquet", index=False)
    aliases = dict(mapping.loc[mapping.original.ne(mapping.canonical), ["original", "canonical"]].itertuples(index=False, name=None))
    old_aliases = json.loads((site_data / "transcript_aliases.json").read_text())
    if any(aliases.get(raw, raw) != canonical for raw, canonical in old_aliases.items()):
        raise ValueError("Full-coverage identifiers conflict with the saved comparison snapshot; regenerate analysis.")
    if sum(s["site_count"] for s in samples) != study["total_scored_sites"]:
        raise ValueError("Observation count does not match the dashboard snapshot.")
    manifest = {"schema_version": 1, "created_at": study["created_at"],
                "model_sha256": study["model_sha256"], "threshold": study["threshold"],
                "samples": samples, "identifier_harmonization": audit,
                "transcript_count": int(mapping.canonical.nunique())}
    temporary = output / "manifest.json.partial"
    temporary.write_text(json.dumps(manifest, indent=2))
    temporary.replace(output / "manifest.json")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=Path("outputs/task2/predictions"))
    parser.add_argument("--site-data", type=Path, default=Path("task2_site/dist/data"))
    parser.add_argument("--output", type=Path, default=Path("task2_site/query_data"))
    args = parser.parse_args()
    bundle(args.predictions, args.site_data, args.output)


if __name__ == "__main__":
    main()
