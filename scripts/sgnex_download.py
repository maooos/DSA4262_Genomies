"""Discover all public m6Anet JSON inputs and optionally download them."""
import argparse
from pathlib import Path
import pandas as pd
from src.sgnex import inventory, download_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("data/sgnex"))
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--reuse-inventory", action="store_true")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--samples", nargs="+")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    path = args.directory / "manifest.csv"
    frame = pd.read_csv(path) if args.reuse_inventory else inventory(args.directory)
    if args.samples:
        unknown = set(args.samples) - set(frame.sample_id)
        if unknown:
            parser.error(f"Unknown samples: {unknown}")
        frame = frame[frame.sample_id.isin(args.samples)]
    print(f"{len(frame)} samples, {frame.cell_line.nunique()} cell lines, {frame.size_bytes.sum()/1e9:.2f} GB", flush=True)
    if args.download:
        download_manifest(frame, args.directory / "raw", args.workers)


if __name__ == "__main__":
    main()
