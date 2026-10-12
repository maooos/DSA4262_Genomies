"""Generate site scores and binary predictions with the saved H04 model."""

import argparse
from pathlib import Path

from src.inference import predict_to_csv


DEFAULT_MODEL = Path(__file__).resolve().parents[1] / "models/h04.joblib"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Project JSON/JSONL (optionally gzip), reads Parquet, or site features")
    parser.add_argument("--output", type=Path, required=True,
                        help="New output CSV; existing files are preserved")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL,
                        help="Saved model bundle (default: models/h04.joblib)")
    parser.add_argument("--input-kind", choices=("reads", "features"), default="reads")
    parser.add_argument("--batch-size", type=int, default=200_000)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    predict_to_csv(args.model, args.input, args.output,
                   input_kind=args.input_kind, batch_size=args.batch_size,
                   threads=args.threads)
    print(f"Metadata: {args.output.with_suffix('.metadata.json')}")


if __name__ == "__main__":
    main()
