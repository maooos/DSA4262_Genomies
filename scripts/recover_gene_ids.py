"""Recover human gene IDs and audit the synthetic origin of data2.

Run with ``python -m scripts.recover_gene_ids --download`` once to obtain the
public references, or omit --download when the references are already cached.
Raw inputs are never modified. All output coordinates explicitly say 0-based.
"""

import argparse
import gzip
import hashlib
import json
import re
import shutil
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ENSEMBL = "https://ftp.ensembl.org/pub/release-91/"
GTF = "Homo_sapiens.GRCh38.91.chr_patch_hapl_scaff.gtf.gz"
CDNA = "Homo_sapiens.GRCh38.cdna.all.fa.gz"
NCRNA = "Homo_sapiens.GRCh38.ncrna.fa.gz"
CURLCAKE = "GSE124309_FASTA_sequences_of_Curlcakes.txt.gz"
EPINANO = "EpiNano_cc.fasta"
PUBLIC_DATA2 = "public_dataset2.json.gz"
REFERENCES = {
    GTF: (ENSEMBL + "gtf/homo_sapiens/" + GTF,
          "6cf849cde2dc588ef3ce0a66dc3e6d5698eff833fd59fef6100531e07f9d69f1"),
    CDNA: (ENSEMBL + "fasta/homo_sapiens/cdna/" + CDNA,
           "a7f0022e884826d70e0a151f6c38547aea8d443eaa9653aaf56e698e23109201"),
    NCRNA: (ENSEMBL + "fasta/homo_sapiens/ncrna/" + NCRNA,
            "0ab3565714f88fb98eb7df7814fe3fed770e2b098d37cd7b83601258f87e75dd"),
    CURLCAKE: ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124309/suppl/" + CURLCAKE,
               "78855ac7309ca92641fe96ba62439973615a4b9d4e1a2f54ccde63fe773ec2cb"),
    EPINANO: ("https://raw.githubusercontent.com/novoalab/EpiNano/master/Reference_sequences/cc.fasta",
              "d522eed26123f5afb9e4a4959cd94dbc162cdf93cbe47e99bc1724cd32c6f381"),
    PUBLIC_DATA2: ("https://raw.githubusercontent.com/bce99/m6A-RNA-Modification-Prediction/main/dataset2.json.gz",
                   "97559d9cad6f433056b446d99dd352dccdae30448d38f4665f4724f6df265a9e"),
}
SITE_KEYS = ["transcript_id", "transcript_position"]
FEATURES = [f"{p}_{f}" for p in ("minus1", "central", "plus1")
            for f in ("dwell", "signal_sd", "signal_mean")]


def sha256(path, decompress=False):
    digest = hashlib.sha256()
    opener = gzip.open if decompress else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_references(directory, download):
    directory.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, (url, expected) in REFERENCES.items():
        path = directory / name
        if not path.exists():
            if not download:
                raise FileNotFoundError(f"Missing {path}; rerun with --download")
            temporary = path.with_suffix(path.suffix + ".part")
            print(f"Downloading {url}", flush=True)
            with urllib.request.urlopen(url, timeout=180) as response, temporary.open("wb") as out:
                shutil.copyfileobj(response, out)
            if sha256(temporary) != expected:
                raise ValueError(f"Reference checksum changed: {url}; inspect before using it")
            temporary.replace(path)
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"Reference checksum mismatch: {path}")
        manifest.append({"path": str(path), "url": url, "sha256": actual,
                         "bytes": path.stat().st_size})
    return manifest


def fasta_records(path):
    """GEO's Curlcake text includes a title before its first FASTA header."""
    opener = gzip.open if path.suffix == ".gz" else open
    header, pieces = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(pieces).upper().replace("U", "T")
                header, pieces = line[1:], []
            elif header is not None and line:
                pieces.append(line)
    if header is not None:
        yield header, "".join(pieces).upper().replace("U", "T")


def read_annotation(path):
    rows = []
    attributes = ("transcript_id", "transcript_version", "gene_id", "gene_version",
                  "gene_name", "transcript_biotype")
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "transcript":
                continue
            values = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            rows.append({**{k: values.get(k, "") for k in attributes}, "seqname": fields[0]})
    mapping = pd.DataFrame(rows).drop_duplicates()
    if mapping.transcript_id.duplicated().any():
        raise ValueError("Annotation has multiple records for a transcript")
    return mapping


def raw_records(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            if len(record) != 1:
                raise ValueError("Expected one transcript per JSON record")
            for tx, positions in record.items():
                for position, contexts in positions.items():
                    if len(contexts) != 1:
                        raise ValueError("Expected one sequence context per site")
                    sequence, reads = next(iter(contexts.items()))
                    yield tx, int(position), sequence, reads


def load_sites(raw_directory, dataset):
    directory = raw_directory / dataset
    info_path = directory / ("data.info.labelled" if dataset == "data0" else "data.info")
    signal_path = directory / f"dataset{dataset[-1]}.json.gz"
    info = pd.read_csv(info_path)
    info["__input_order"] = np.arange(len(info))
    sites = pd.DataFrame(
        ((tx, p, s, len(r)) for tx, p, s, r in raw_records(signal_path)),
        columns=[*SITE_KEYS, "sequence", "raw_n_reads"],
    )
    merged = info.merge(sites, on=SITE_KEYS, how="outer", validate="one_to_one", indicator=True)
    if not merged["_merge"].eq("both").all():
        raise ValueError(f"{dataset}: metadata and raw JSON have different site keys")
    if "n_reads" in merged and not merged.n_reads.eq(merged.raw_n_reads).all():
        raise ValueError(f"{dataset}: metadata and JSON read counts differ")
    # An outer merge sorts keys in recent pandas; restore the input metadata order.
    merged = merged.sort_values("__input_order").drop(columns=["_merge", "__input_order"])
    return merged.reset_index(drop=True), {
        "metadata_path": str(info_path), "metadata_sha256": sha256(info_path),
        "signal_path": str(signal_path), "signal_sha256": sha256(signal_path),
        "sites": len(sites), "transcripts": sites.transcript_id.nunique(),
        "read_rows": int(sites.raw_n_reads.sum()),
    }


def map_human_sites(sites, annotation, references):
    had_gene_ids = "gene_id" in sites
    if had_gene_ids:
        sites = sites.rename(columns={"gene_id": "original_gene_id"})
    result = sites.merge(annotation, on="transcript_id", how="left", validate="many_to_one")
    if result.gene_id.isna().any():
        raise ValueError("Unmapped human transcript IDs")
    if had_gene_ids and not result.original_gene_id.eq(result.gene_id).all():
        raise ValueError("Reference gene IDs conflict with original labels")
    matches, header_matches = [], []
    for row in result.itertuples():
        header, sequence = references[row.transcript_id]
        pos = row.transcript_position
        matches.append(pos >= 3 and sequence[pos - 3:pos + 4] == row.sequence)
        gene = re.search(r"(?:^| )gene:(\S+)", header)
        header_matches.append(gene is not None and gene.group(1).split(".")[0] == row.gene_id)
    result["sequence_matches_reference"] = matches
    result["fasta_gene_matches_gtf"] = header_matches
    if not all(header_matches):
        raise ValueError("FASTA gene headers disagree with GTF")
    result["gene_id_status"] = "confirmed_ensembl91_mapping"
    result["annotation_source"] = GTF
    return result


def unique_ordered_embedding(target, reference):
    """Return a unique subsequence alignment, conditional on preserved order.

    Repeated short contexts are NOT treated as unique sequence matches. The
    earliest and latest feasible embeddings must agree for every target site.
    """
    early, late = [], []
    j = 0
    for sequence in target:
        while j < len(reference) and reference[j] != sequence:
            j += 1
        if j == len(reference):
            raise ValueError("No complete ordered reference match")
        early.append(j)
        j += 1
    j = len(reference) - 1
    for sequence in reversed(target):
        while j >= 0 and reference[j] != sequence:
            j -= 1
        if j < 0:
            raise ValueError("No complete reverse ordered reference match")
        late.append(j)
        j -= 1
    if early != late[::-1]:
        raise ValueError("Ordered reference mapping is ambiguous")
    return early


def map_synthetic_sites(sites, geo_refs, epinano_refs):
    ordered = sites[sites.transcript_id == "tx_id_0"].sort_values("transcript_position")
    for _, group in sites.groupby("transcript_id"):
        observed = group.sort_values("transcript_position")[["transcript_position", "sequence"]]
        if not observed.reset_index(drop=True).equals(
            ordered[["transcript_position", "sequence"]].reset_index(drop=True)
        ):
            raise ValueError("Mixture groups do not share the same ordered contexts")
    all_sites, candidates, offsets = [], defaultdict(list), {}
    for name, sequence in geo_refs.items():
        alternative_name = next(k for k in epinano_refs if k.lower() == name.lower())
        offset = epinano_refs[alternative_name].find(sequence)
        if offset < 0 or epinano_refs[alternative_name].find(sequence, offset + 1) != -1:
            raise ValueError("GEO reference is not a unique substring of the EpiNano reference")
        offsets[name] = {"epinano_reference_id": alternative_name, "offset": offset}
        for start in range(len(sequence) - 6):
            context = sequence[start:start + 7]
            if re.fullmatch(r".[AGT][AG]AC[ACT].", context):
                entry = (name, start + 3, context)
                all_sites.append(entry)
                candidates[context].append(entry)
    embedding = unique_ordered_embedding(ordered.sequence.tolist(), [s for _, _, s in all_sites])
    mapping, candidate_rows = [], []
    for site, index in zip(ordered.itertuples(), embedding):
        name, position, sequence = all_sites[index]
        alt = offsets[name]
        mapping.append({
            "transcript_position": site.transcript_position,
            "geo_reference_id": name, "geo_reference_center_0based": position,
            "epinano_reference_id": alt["epinano_reference_id"],
            "epinano_reference_center_0based": position + alt["offset"],
            "sequence_only_candidate_count": len(candidates[sequence]),
            "reference_mapping_status": ("unique_7mer_in_curlcake_reference" if len(candidates[sequence]) == 1
                                         else "inferred_unique_ordered_alignment"),
        })
        for candidate_name, candidate_position, _ in candidates[sequence]:
            candidate_rows.append({"transcript_position": site.transcript_position,
                                   "sequence": sequence, "geo_reference_id": candidate_name,
                                   "geo_reference_center_0based": candidate_position})
    result = sites.merge(pd.DataFrame(mapping), on="transcript_position", validate="many_to_one")
    # Synthetic constructs have reference identifiers, not human Ensembl genes.
    result["gene_id"] = pd.NA
    result["gene_id_status"] = "not_applicable_synthetic_construct"
    result["reference_source"] = "GEO:GSE124309; EpiNano:Reference_sequences/cc.fasta"
    details = {
        "origin_assessment": "Strong sequence and mixture evidence for synthetic Curlcake RNA",
        "gene_id_applicability": "Synthetic constructs have no human ENSG gene_id",
        "coordinate_warning": "Positions refer to named downloaded FASTAs; the original processing FASTA is unknown",
        "ordered_mapping_assumption": "Original sites retained order within four concatenated Curlcake references",
        "complete_ordered_mapping_is_unique": True,
        "reference_drach_sites": len(all_sites), "target_contexts": len(ordered),
        "unique_7mer_contexts": int(sum(len(candidates[s]) == 1 for s in ordered.sequence)),
        "order_resolved_contexts": int(sum(len(candidates[s]) > 1 for s in ordered.sequence)),
        "reference_offsets": offsets,
        "excluded_reference_sites": [dict(zip(("reference_id", "center_0based", "sequence"), site))
                                     for i, site in enumerate(all_sites) if i not in embedding],
        "per_mixture_construct_counts": pd.DataFrame(mapping).geo_reference_id.value_counts().to_dict(),
    }
    return result, pd.DataFrame(candidate_rows), details


def verify_mixture_signals(path):
    blocks = []
    for tx, position, sequence, reads in raw_records(path):
        block = pd.DataFrame(np.asarray(reads, dtype=float), columns=FEATURES)
        block["transcript_id"], block["transcript_position"], block["sequence"] = tx, position, sequence
        blocks.append(block)
    reads = pd.concat(blocks, ignore_index=True)
    keys = ["transcript_position", "sequence", *FEATURES]
    endpoints = reads[reads.transcript_id.isin(["tx_id_0", "tx_id_6"])].rename(
        columns={"transcript_id": "endpoint_group"}
    ).drop_duplicates()
    if endpoints.duplicated(keys).any():
        raise ValueError("Endpoint groups share indistinguishable signals")
    # Exact equality of position, sequence and all nine numeric features, not a hash-only match.
    matched = reads.merge(endpoints, on=keys, how="left", validate="many_to_one")
    rows = []
    for tx, group in matched.groupby("transcript_id"):
        a = int(group.endpoint_group.eq("tx_id_0").sum())
        b = int(group.endpoint_group.eq("tx_id_6").sum())
        rows.append({"transcript_id": tx, "read_rows": len(group), "from_endpoint_1": a,
                     "from_endpoint_0": b, "unmatched": len(group) - a - b,
                     "observed_endpoint_1_fraction": a / len(group)})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--reference-dir", type=Path, default=Path("data/reference/gene_id_recovery"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/gene_id_recovery"))
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    manifest = ensure_references(args.reference_dir, args.download)
    annotation = read_annotation(args.reference_dir / GTF)
    datasets, summary = {}, {"generated_at_utc": datetime.now(timezone.utc).isoformat(), "datasets": {}}
    for name in ("data0", "data1", "data2"):
        print(f"Reading and checking raw {name}", flush=True)
        datasets[name], summary["datasets"][name] = load_sites(args.raw_dir, name)
    needed = set(datasets["data0"].transcript_id) | set(datasets["data1"].transcript_id)
    sequences = {}
    for filename in (CDNA, NCRNA):
        for header, sequence in fasta_records(args.reference_dir / filename):
            identifier = header.split()[0].split(".")[0]
            if identifier in needed:
                if identifier in sequences:
                    raise ValueError("Duplicate transcript in reference FASTAs")
                sequences[identifier] = (header, sequence)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for name in ("data0", "data1"):
        result = map_human_sites(datasets[name], annotation, sequences)
        result.to_csv(output / f"{name}.gene_ids.csv", index=False)
        summary["datasets"][name].update({
            "mapped_sites": len(result), "mapped_transcripts": result.transcript_id.nunique(),
            "gene_count": result.gene_id.nunique(),
            "sequence_matches": int(result.sequence_matches_reference.sum()),
            "sequence_mismatches": int((~result.sequence_matches_reference).sum()),
            "fasta_gene_gtf_matches": int(result.fasta_gene_matches_gtf.sum()),
            "original_gene_id_conflicts": 0 if name == "data0" else None,
        })
        if not result.sequence_matches_reference.all():
            result[~result.sequence_matches_reference].to_csv(output / f"{name}.sequence_mismatches.csv", index=False)
        results.append(result.assign(dataset=name))
        print(name, summary["datasets"][name], flush=True)
    geo = {h.split()[0]: s for h, s in fasta_records(args.reference_dir / CURLCAKE)}
    epinano = {h.split()[0]: s for h, s in fasta_records(args.reference_dir / EPINANO)}
    result, candidates, synthetic_details = map_synthetic_sites(datasets["data2"], geo, epinano)
    local_digest = sha256(args.raw_dir / "data2/dataset2.json.gz", decompress=True)
    public_digest = sha256(args.reference_dir / PUBLIC_DATA2, decompress=True)
    synthetic_details.update({"local_uncompressed_sha256": local_digest,
                              "public_uncompressed_sha256": public_digest,
                              "public_copy_identical": local_digest == public_digest,
                              "mixture_signal_evidence": verify_mixture_signals(args.raw_dir / "data2/dataset2.json.gz")})
    summary["datasets"]["data2"].update(synthetic_details)
    result.to_csv(output / "data2.reference_mapping.csv", index=False)
    candidates.to_csv(output / "data2.reference_candidates.csv", index=False)
    pd.DataFrame(synthetic_details["mixture_signal_evidence"]).to_csv(output / "data2.mixture_signal_evidence.csv", index=False)
    results.append(result.assign(dataset="data2"))
    all_sites = pd.concat(results, ignore_index=True)
    all_sites.to_csv(output / "all_sites.gene_ids.csv", index=False)
    annotation[annotation.transcript_id.isin(needed)].to_csv(output / "human_transcript_to_gene.csv", index=False)
    summary["human_gene_ids_recovered"] = int(all_sites.gene_id.notna().sum())
    summary["synthetic_sites_without_biological_gene_id"] = len(result)
    summary["sources"] = manifest
    for filename in ("data0.gene_ids.csv", "data1.gene_ids.csv", "data2.reference_mapping.csv",
                     "data2.reference_candidates.csv", "data2.mixture_signal_evidence.csv",
                     "all_sites.gene_ids.csv", "human_transcript_to_gene.csv"):
        summary.setdefault("outputs", []).append({"path": str(output / filename),
                                                  "sha256": sha256(output / filename)})
    (output / "recovery.audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"human_gene_ids_recovered": summary["human_gene_ids_recovered"],
                      "synthetic_sites": len(result), "output_directory": str(output)}, indent=2))


if __name__ == "__main__":
    main()
