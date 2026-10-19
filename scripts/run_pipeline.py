"""Run existing inference sequentially for a CSV sample manifest."""

import argparse
import csv
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
import uuid

from scripts.export_submission import export_submission

ROOT = Path(__file__).resolve().parents[1]
CODE_FILES = [
    ROOT / "src/inference.py", ROOT / "src/feature_engineering.py",
    ROOT / "src/data_parser.py", ROOT / "scripts/predict.py",
    Path(__file__).resolve(), ROOT / "scripts/export_submission.py",
]


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def output_lock(output_dir):
    """Acquire an exclusive advisory lock shared by cooperating pipeline runs."""
    import fcntl

    output_dir.mkdir(parents=True, exist_ok=True)
    lock_path = output_dir / ".run_pipeline.lock"
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another pipeline already holds {lock_path}") from exc
        lock.seek(0)
        lock.truncate()
        lock.write(json.dumps({"pid": os.getpid(), "started_at": datetime.now(timezone.utc).isoformat()}))
        lock.flush()
        yield


def load_manifest(path):
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Manifest not found: {path}. Create a CSV with header "
            "sample_id,input_path, or try --manifest examples/pipeline_samples.csv "
            "for the included example."
        )
    samples, seen = [], set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["sample_id", "input_path"]:
            raise ValueError("Manifest header must be sample_id,input_path")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Manifest row {reader.line_num} must contain exactly two fields")
            sample = row["sample_id"]
            if not re.fullmatch(r"[A-Za-z0-9_-]+", sample or "") or sample in seen:
                raise ValueError("Sample IDs must be unique and filename-safe")
            if not row["input_path"]:
                raise ValueError("Missing input path")
            source = Path(row["input_path"]).expanduser()
            source = (path.parent / source).resolve() if not source.is_absolute() else source.resolve()
            if not source.is_file():
                raise FileNotFoundError(source)
            seen.add(sample)
            samples.append((sample, source))
    if not samples:
        raise ValueError("Manifest has no samples")
    return samples


def environment_info():
    revision = None
    try:
        git_dir = ROOT / ".git"
        head = (git_dir / "HEAD").read_text(encoding="utf-8").strip()
        if head.startswith("ref: "):
            ref = head[5:]
            loose_ref = git_dir / ref
            revision = loose_ref.read_text(encoding="utf-8").strip() if loose_ref.is_file() else None
        else:
            revision = head
    except (OSError, UnicodeError):
        pass
    return {"git_revision": revision, "python_version": platform.python_version(),
            "platform": "-".join((platform.system(), platform.release(), platform.machine()))}


def sample_identity(sample, source, model_hash):
    return {"sample_id": sample, "input_path": str(source), "input_sha256": digest(source),
            "model_sha256": model_hash,
            "code_sha256": {str(p.relative_to(ROOT)): digest(p) for p in CODE_FILES}}


def verified_complete(folder, identity):
    record_path = folder / "complete.json"
    if not record_path.is_file():
        return False
    try:
        record = json.loads(record_path.read_text(encoding="utf-8"))
        outputs = record["outputs"]
        return (record.get("identity") == identity and
                all((folder / name).is_file() and outputs.get(name) == digest(folder / name)
                    for name in ("predictions.csv", "predictions.metadata.json", "submission.csv")))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_pipeline(manifest, model, output_dir, *, threads=4, batch_size=200_000,
                 resume=False, retry_failed=False, dry_run=False):
    if threads < 1 or batch_size < 1:
        raise ValueError("Thread count and batch size must be positive")
    samples = load_manifest(manifest)
    model, output_dir = Path(model).resolve(), Path(output_dir).resolve()
    if not model.is_file():
        raise FileNotFoundError(model)
    if dry_run:
        for sample, source in samples:
            print(f"{sample}: {source} -> {output_dir / sample}")
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    with output_lock(output_dir):
        started = datetime.now(timezone.utc).isoformat()
        run_started = time.monotonic()
        model_hash = digest(model)
        environment = environment_info()
        run_id = uuid.uuid4().hex
        failures = 0
        for sample, source in samples:
            folder = output_dir / sample
            identity = sample_identity(sample, source, model_hash)
            if folder.exists():
                if (resume or retry_failed) and verified_complete(folder, identity):
                    print(f"{sample}: verified completed outputs; skipped", flush=True)
                    continue
                if retry_failed and (folder / "failed.json").is_file() and not (folder / "complete.json").exists():
                    shutil.rmtree(folder)
                else:
                    raise FileExistsError(
                        f"{folder}: existing, incomplete, or changed run; use --resume for verified outputs "
                        "or --retry-failed for a recorded failure"
                    )

            work = output_dir / f".{sample}.work-{run_id}"
            work.mkdir()
            outputs = [work / "predictions.csv", work / "predictions.metadata.json", work / "submission.csv"]
            command = [sys.executable, "-m", "scripts.predict", "--model", str(model),
                       "--input", str(source), "--input-kind", "reads", "--output", str(outputs[0]),
                       "--threads", str(threads), "--batch-size", str(batch_size)]
            sample_started = time.monotonic()
            started_at = datetime.now(timezone.utc).isoformat()
            print(f"{sample}: predicting (log: {folder / 'run.log'})", flush=True)
            log_path = work / "run.log"
            try:
                with log_path.open("x", encoding="utf-8") as log:
                    subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
                sites = export_submission(outputs[0], outputs[2])
                if digest(source) != identity["input_sha256"] or digest(model) != model_hash:
                    raise ValueError("Input or model changed during prediction")
                record = {
                    "identity": identity, "site_count": sites, "command": command,
                    "started_at": started_at, "completed_at": datetime.now(timezone.utc).isoformat(),
                    "duration_seconds": round(time.monotonic() - sample_started, 3),
                    "threads": threads, "batch_size": batch_size, "environment": environment,
                    "outputs": {name: digest(work / name) for name in
                                ("predictions.csv", "predictions.metadata.json", "submission.csv")},
                }
                atomic_json(work / "complete.json", record)
                os.replace(work, folder)
                print(f"{sample}: complete, {sites} sites", flush=True)
            except Exception as exc:
                failures += 1
                failure = {"sample_id": sample, "error": str(exc), "traceback": traceback.format_exc(),
                           "command": command, "started_at": started_at,
                           "failed_at": datetime.now(timezone.utc).isoformat(),
                           "duration_seconds": round(time.monotonic() - sample_started, 3),
                           "threads": threads, "batch_size": batch_size, "environment": environment}
                atomic_json(work / "failed.json", failure)
                # Keep the final sample path stable and inspectable; retry_failed explicitly replaces it.
                os.replace(work, folder)
                print(f"{sample}: failed: {exc}", file=sys.stderr, flush=True)

        run_record = {"run_id": run_id, "started_at": started,
                      "completed_at": datetime.now(timezone.utc).isoformat(),
                      "duration_seconds": round(time.monotonic() - run_started, 3),
                      "manifest": str(Path(manifest).resolve()), "model": str(model),
                      "threads": threads, "batch_size": batch_size, "resume": resume,
                      "retry_failed": retry_failed, "environment": environment,
                      "sample_count": len(samples), "failure_count": failures}
        atomic_json(output_dir / f"run-{run_id}.json", run_record)
        return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", type=Path, default=ROOT / "models/h04.joblib")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=200_000)
    parser.add_argument("--resume", action="store_true", help="Skip only checksum-verified completed samples")
    parser.add_argument("--retry-failed", action="store_true", help="Retry samples with failed.json and no complete.json")
    parser.add_argument("--dry-run", action="store_true", help="Validate paths and print samples; write nothing and run no predictions")
    sys.exit(run_pipeline(**vars(parser.parse_args())))


if __name__ == "__main__":
    main()
