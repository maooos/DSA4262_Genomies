"""Public SG-NEx inventory and verified, resumable downloads (no AWS account)."""

import hashlib
import http.client
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

BUCKET = "https://sg-nex-data.s3.ap-southeast-1.amazonaws.com"
PREFIX = "data/processed_data/m6Anet/"
METADATA_URL = "https://raw.githubusercontent.com/GoekeLab/sg-nex-data/master/docs/samples.tsv"
NS = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}


def read_url(url):
    with urllib.request.urlopen(url, timeout=90) as response:
        return response.read()


def inventory(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    metadata = read_url(METADATA_URL)
    (directory / "samples.tsv").write_bytes(metadata)
    samples = pd.read_csv(directory / "samples.tsv", sep="\t").set_index("sample_alias")
    objects, token, page = [], None, 0
    while True:
        query = {"list-type": "2", "prefix": PREFIX}
        if token:
            query["continuation-token"] = token
        raw = read_url(BUCKET + "/?" + urllib.parse.urlencode(query))
        (directory / f"listing_{page}.xml").write_bytes(raw)
        root = ET.fromstring(raw)
        for item in root.findall("s:Contents", NS):
            key = item.findtext("s:Key", namespaces=NS)
            if not key.endswith("/data.json"):
                continue
            sample_id = key.split("/")[-2]
            match = re.fullmatch(r"SGNex_(.+)_directRNA_replicate(\d+)_run(\d+)", sample_id)
            if not match:
                raise ValueError(f"Unrecognised sample ID: {sample_id}")
            cell_line, replicate, run = match.groups()
            if sample_id in samples.index:
                if samples.loc[sample_id, "protocol"] != "directRNA":
                    raise ValueError(f"Non-directRNA sample: {sample_id}")
                if samples.loc[sample_id, "cellLine"] != cell_line:
                    raise ValueError(f"Cell-line metadata mismatch: {sample_id}")
            objects.append(dict(
                sample_id=sample_id, cell_line=cell_line, replicate=int(replicate), run=int(run),
                url=BUCKET + "/" + urllib.parse.quote(key),
                size_bytes=int(item.findtext("s:Size", namespaces=NS)),
                etag=item.findtext("s:ETag", namespaces=NS).strip('"'),
                last_modified=item.findtext("s:LastModified", namespaces=NS),
                metadata_matched=sample_id in samples.index,
            ))
        if root.findtext("s:IsTruncated", namespaces=NS) != "true":
            break
        token = root.findtext("s:NextContinuationToken", namespaces=NS)
        if not token:
            raise ValueError("Truncated listing has no continuation token.")
        page += 1
    manifest = pd.DataFrame(objects).sort_values("sample_id")
    if manifest.empty or manifest.sample_id.duplicated().any():
        raise ValueError("Inventory is empty or contains duplicate samples.")
    manifest.to_csv(directory / "manifest.csv", index=False)
    missing = samples[samples.protocol.eq("directRNA") & ~samples.index.isin(manifest.sample_id)]
    missing.to_csv(directory / "directRNA_without_processed_json.csv")
    (directory / "inventory.json").write_text(json.dumps({
        "created_at": datetime.now(timezone.utc).isoformat(), "bucket": BUCKET,
        "prefix": PREFIX, "metadata_url": METADATA_URL, "processed_samples": len(manifest),
        "cell_lines": sorted(manifest.cell_line.unique()),
        "size_bytes": int(manifest.size_bytes.sum()),
        "directRNA_samples_without_processed_json": len(missing),
    }, indent=2))
    return manifest


def download_sample(row, directory):
    """Resume a byte-range download, verify length/ETag and record SHA256."""
    directory = Path(directory) / row["sample_id"]
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "data.json"
    audit = directory / "download.json"
    expected = int(row["size_bytes"])
    if target.exists():
        if not audit.exists():
            raise ValueError(f"Unverified existing download: {target}")
        old = json.loads(audit.read_text())
        if target.stat().st_size != expected or old["etag"] != row["etag"]:
            raise ValueError(f"Stale download: {target}")
        return target
    partial = directory / "data.json.partial"
    partial_audit = directory / "download.partial.json"
    identity = {key: row[key] for key in ("url", "size_bytes", "etag")}
    if partial.exists():
        if not partial_audit.exists() or json.loads(partial_audit.read_text()) != identity:
            raise ValueError(f"Partial download belongs to an unverified object: {partial}")
    partial_audit.write_text(json.dumps(identity))
    for attempt in range(4):
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > expected:
            raise ValueError(f"Oversized partial: {partial}")
        if offset == expected:
            break
        headers = {"If-Match": '"' + row["etag"] + '"'}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        try:
            request = urllib.request.Request(row["url"], headers=headers)
            with urllib.request.urlopen(request, timeout=90) as response:
                if response.headers.get("ETag", "").strip('"') != row["etag"]:
                    raise ValueError("Object ETag changed after inventory.")
                if offset and (response.status != 206 or not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-")):
                    raise ValueError("Server did not honour the requested download range.")
                with partial.open("ab" if offset else "wb") as handle:
                    while block := response.read(4 * 1024 * 1024):
                        handle.write(block)
            if partial.stat().st_size == expected:
                break
        except (OSError, TimeoutError, http.client.IncompleteRead):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    if partial.stat().st_size != expected:
        raise ValueError(f"Incomplete download: {partial}")
    with partial.open("rb") as handle:
        sha = hashlib.file_digest(handle, "sha256").hexdigest()
    audit.write_text(json.dumps({**row, "sha256": sha}, indent=2))
    partial.rename(target)
    partial_audit.unlink(missing_ok=True)
    print(f"Downloaded {row['sample_id']} ({expected / 1e9:.2f} GB)", flush=True)
    return target


def download_manifest(manifest, directory, workers=2):
    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(lambda row: download_sample(row, directory), manifest.to_dict("records")))
