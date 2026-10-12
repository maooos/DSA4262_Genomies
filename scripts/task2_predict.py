"""Score each inventoried SG-NEx run with H04; safely resume verified results."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import pandas as pd

from src.inference import predict_to_csv


def run(manifest_path, raw_directory, output_directory, model_path,
        threads=4, batch_size=200_000, wait_seconds=0, refresh_parquet=False):
    manifest = pd.read_csv(manifest_path)
    if manifest.sample_id.duplicated().any():
        raise ValueError("Duplicate sample IDs in manifest.")
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    with open(model_path, "rb") as handle:
        model_sha = hashlib.file_digest(handle, "sha256").hexdigest()
    pending = manifest.to_dict("records")
    deadline = time.monotonic() + wait_seconds
    while pending:
        progressed = False
        for sample in pending[:]:
            sample_id = sample["sample_id"]
            source = Path(raw_directory) / sample_id / "data.json"
            output = output_directory / f"{sample_id}.csv"
            meta_path = output.with_suffix(".metadata.json")
            if not source.exists():
                continue
            if output.exists():
                old = json.loads(meta_path.read_text())
                stat = source.stat()
                if (old["model_sha256"] != model_sha
                        or old["input_path"] != str(source.resolve())
                        or old["input_size_bytes"] != stat.st_size
                        or old["input_mtime_ns"] != stat.st_mtime_ns
                        or not {"sequence", "n_reads"} <= set(old["output_columns"])):
                    raise ValueError(f"Existing result does not match this run: {output}")
                print(f"Verified existing predictions: {sample_id}", flush=True)
            else:
                predict_to_csv(model_path, source, output, input_kind="reads",
                               batch_size=batch_size, threads=threads, include_site_metadata=True)
            parquet = output.with_suffix(".parquet")
            if refresh_parquet or not parquet.exists():
                frame = pd.read_csv(output, float_precision="round_trip")
                for key in ("sample_id", "cell_line", "replicate", "run"):
                    frame[key] = sample[key]
                temp = parquet.with_suffix(".parquet.partial")
                frame.to_parquet(temp, index=False)
                temp.rename(parquet)
            pending.remove(sample)
            progressed = True
        if not pending:
            break
        if time.monotonic() >= deadline:
            raise FileNotFoundError("Missing downloads: " + ", ".join(s["sample_id"] for s in pending))
        if not progressed:
            print(f"Waiting for {len(pending)} downloads...", flush=True)
            time.sleep(20)
    return len(manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("data/sgnex/manifest.csv"))
    parser.add_argument("--raw-directory", type=Path, default=Path("data/sgnex/raw"))
    parser.add_argument("--output-directory", type=Path, default=Path("outputs/task2/predictions"))
    parser.add_argument("--model", type=Path, default=Path("models/h04.joblib"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=200_000)
    parser.add_argument("--wait-seconds", type=int, default=0)
    parser.add_argument("--refresh-parquet", action="store_true",
                        help="Rebuild derived Parquet tables from verified CSV outputs without rescoring")
    args = parser.parse_args()
    run(args.manifest, args.raw_directory, args.output_directory, args.model,
        args.threads, args.batch_size, args.wait_seconds, args.refresh_parquet)


if __name__ == "__main__":
    main()
