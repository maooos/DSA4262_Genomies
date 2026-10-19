"""Validate prediction scores and export the handout's three-column CSV."""

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile

FIELDS = ["transcript_id", "transcript_position", "score"]


def export_submission(input_path, output_path):
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count, seen = 0, set()
    descriptor, temporary = tempfile.mkstemp(dir=output_path.parent, suffix=".partial")
    try:
        with os.fdopen(descriptor, "w", newline="", encoding="utf-8") as target, open(input_path, newline="", encoding="utf-8") as source:
            reader = csv.DictReader(source)
            if not set(FIELDS) <= set(reader.fieldnames or []):
                raise ValueError("Missing submission columns")
            if "dataset" in reader.fieldnames:
                raise ValueError("Export each dataset separately before submission")
            writer = csv.DictWriter(target, fieldnames=FIELDS)
            writer.writeheader()
            for row in reader:
                transcript = row["transcript_id"]
                position = int(row["transcript_position"])
                score = float(row["score"])
                key = (transcript, position)
                if not transcript.strip() or key in seen or not math.isfinite(score) or not 0 <= score <= 1:
                    raise ValueError("Empty identifier, duplicate site, or invalid score")
                seen.add(key)
                writer.writerow({field: row[field] for field in FIELDS})
                count += 1
            if not count:
                raise ValueError("No prediction sites")
        # Hard link publishes the validated file without replacing an existing output.
        os.link(temporary, output_path)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(f"Exported {export_submission(args.input, args.output)} sites")


if __name__ == "__main__":
    main()
