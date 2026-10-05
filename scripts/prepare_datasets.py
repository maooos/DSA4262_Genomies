"""Raw -> parse_data.py -> recovered annotations -> full annotated datasets.

Run from the repository root. Existing stages are reused only when their
recorded checksums match; stale or partial outputs are never silently reused.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from src.data_parser import file_sha256
from src.gene_annotation import annotate_parquet


MAPPING_FILES = {"data0": "data0.gene_ids.csv", "data1": "data1.gene_ids.csv",
                 "data2": "data2.reference_mapping.csv"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--reference-dir", type=Path, default=Path("data/reference/gene_id_recovery"))
    parser.add_argument("--datasets", nargs="+", choices=list(MAPPING_FILES), default=list(MAPPING_FILES))
    parser.add_argument("--download", action="store_true", help="Download missing public references for recovery")
    args = parser.parse_args(argv)
    parsed_dir = args.processed_dir / "parsed"
    output_dir = args.processed_dir / "annotated"
    recovery_dir = args.processed_dir / "gene_id_recovery"
    inputs = {}

    # Invoke the existing parsing CLI, without giving it the recovered mapping:
    # the supplementary annotation is deliberately a separate stage.
    for dataset in args.datasets:
        directory = args.raw_dir / dataset
        signals = directory / f"dataset{dataset[-1]}.json.gz"
        labels = directory / ("data.info.labelled" if dataset == "data0" else "data.info")
        expected = {"signal_sha256": file_sha256(signals), "label_sha256": file_sha256(labels)}
        inputs[dataset] = expected
        parsed = parsed_dir / f"{dataset}_reads.parquet"
        audit_path = parsed.with_suffix(".audit.json")
        if not parsed.exists() and not audit_path.exists():
            print(f"[parse] {dataset}: raw -> {parsed}", flush=True)
            subprocess.run([sys.executable, "-m", "scripts.parse_data", "--signals", str(signals),
                            "--labels", str(labels), "--output", str(parsed)], check=True)
        if not parsed.exists() or not audit_path.exists():
            raise ValueError(f"Incomplete parse stage: {parsed}; choose a new --processed-dir")
        audit = json.loads(audit_path.read_text())
        if (any(audit.get(k) != v for k, v in expected.items()) or
                audit.get("output_sha256") != file_sha256(parsed) or audit.get("max_sites_requested") is not None):
            raise ValueError(f"Stale/partial parsed dataset: {parsed}; choose a new --processed-dir")
        print(f"[parse] {dataset}: verified {audit['parsed_read_rows']:,} reads", flush=True)

    recovery_audit_path = recovery_dir / "recovery.audit.json"
    recovery = json.loads(recovery_audit_path.read_text()) if recovery_audit_path.exists() else {}

    def recovery_is_current():
        outputs = {Path(item["path"]).name: item["sha256"] for item in recovery.get("outputs", [])}
        for name in args.datasets:
            entry = recovery.get("datasets", {}).get(name, {})
            path = recovery_dir / MAPPING_FILES[name]
            if (entry.get("signal_sha256") != inputs[name]["signal_sha256"] or
                    entry.get("metadata_sha256") != inputs[name]["label_sha256"] or
                    not path.exists() or outputs.get(path.name) != file_sha256(path)):
                return False
        return True

    if not recovery_is_current():
        print("[reference] Building verified gene/source mappings", flush=True)
        command = [sys.executable, "-m", "scripts.recover_gene_ids", "--raw-dir", str(args.raw_dir),
                   "--reference-dir", str(args.reference_dir), "--output-dir", str(recovery_dir)]
        if args.download:
            command.append("--download")
        subprocess.run(command, check=True)
        recovery = json.loads(recovery_audit_path.read_text())
        if not recovery_is_current():
            raise ValueError("Recovery outputs do not match the current raw inputs")

    results = []
    for dataset in args.datasets:
        parsed = parsed_dir / f"{dataset}_reads.parquet"
        mapping = recovery_dir / MAPPING_FILES[dataset]
        output = output_dir / f"{dataset}_reads.parquet"
        audit_path = output.with_suffix(".audit.json")
        if output.exists() or audit_path.exists():
            if not output.exists() or not audit_path.exists():
                raise ValueError(f"Incomplete annotation stage: {output}")
            result = json.loads(audit_path.read_text())
            expected = {"input_sha256": file_sha256(parsed), "mapping_sha256": file_sha256(mapping),
                        "output_sha256": file_sha256(output)}
            if any(result.get(k) != v for k, v in expected.items()):
                raise ValueError(f"Stale annotation output: {output}; choose a new --processed-dir")
        else:
            print(f"[annotate] {dataset}: merging annotations into all read rows", flush=True)
            result = annotate_parquet(parsed, mapping, output, dataset=dataset)
        results.append(result)
        print(f"[ready] {output}: {result['read_rows']:,} reads, {result['sites']:,} sites; "
              f"filled gene_id on {result['gene_id_filled_reads']:,} reads", flush=True)
    print("Final modelling inputs are the annotated/ Parquet files above.", flush=True)
    return results


if __name__ == "__main__":
    main()
