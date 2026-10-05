"""Attach gene annotations; binarize data2 labels while retaining label_original."""

import argparse
import json
from pathlib import Path

from src.gene_annotation import annotate_parquet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["data0", "data1", "data2"], required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=100_000)
    args = parser.parse_args()
    result = annotate_parquet(args.input, args.mapping, args.output,
                              dataset=args.dataset, batch_size=args.batch_size)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
