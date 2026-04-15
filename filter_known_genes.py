#!/usr/bin/env python3
"""
filter_known_genes.py

Pre-filter host transcriptomics to keep only genes recognised by mygene.info
(i.e. genes that exist in Ensembl / NCBI Entrez / HGNC).  Unrecognised gene
names — typically non-standard identifiers, pseudoautosomal duplicates, or
retired symbols — are dropped.

By default keeps ALL matched genes regardless of type.  Use --keep-types to
restrict to specific biotypes (e.g. protein-coding lncRNA).

Offline / cached operation
--------------------------
If an annotation TSV already exists at --out-annotations (default
Data/gene_annotations.tsv) and --refresh is NOT set, the cached table is
used directly and no network calls are made. This lets the pipeline run
fully offline after the first successful annotation pass; commit the
generated Data/gene_annotations.tsv to ship a self-contained repo.

Outputs
-------
  Data/host_tx_counts_filtered.tsv   — filtered transcript counts
  Data/gene_annotations.tsv          — full annotation lookup table
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.request

import pandas as pd

MYGENE_URL  = "https://mygene.info/v3/query"
BATCH_SIZE  = 1000
RETRY_MAX   = 4
RETRY_DELAY = 10


def trim_version(gene: str) -> str:
    """Remove trailing .<integer> version suffix, e.g. AC000003.1 -> AC000003."""
    return re.sub(r"\.\d+$", "", gene)


def query_mygene(genes: list[str], batch_label: str = "") -> dict[str, dict]:
    """
    Query mygene.info for gene type and description.
    Returns {queried_gene: {"type": ..., "desc": ...}} for matched genes.
    """
    payload = json.dumps({
        "q":       genes,
        "scopes":  ["symbol", "alias"],
        "fields":  ["symbol", "type_of_gene", "name", "query"],
        "species": "human",
        "size":    len(genes),
    }).encode("utf-8")

    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for attempt in range(1, RETRY_MAX + 1):
        try:
            req = urllib.request.Request(
                MYGENE_URL, data=payload, headers=headers, method="POST"
            )
            with urllib.request.urlopen(req, timeout=90) as resp:
                body = json.loads(resp.read().decode("utf-8"))
            break
        except Exception as exc:
            print(f"  [{batch_label} attempt {attempt}/{RETRY_MAX}] Error: {exc}",
                  file=sys.stderr)
            if attempt == RETRY_MAX:
                print(f"  [{batch_label}] All attempts failed - batch skipped.",
                      file=sys.stderr)
                return {}
            time.sleep(RETRY_DELAY)

    result: dict[str, dict] = {}
    for hit in body:
        if hit.get("notfound"):
            continue
        query_key = hit.get("query", "")
        tog  = hit.get("type_of_gene", "unknown")
        name = hit.get("name", "")
        if not query_key:
            continue
        if query_key not in result or tog == "protein-coding":
            result[query_key] = {"type": tog, "desc": name}
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Pre-filter host transcriptomics to known/recognised genes.")
    parser.add_argument("--input", default="Data/host_tx_counts.tsv",
                        help="Raw host transcript counts (genes as rows)")
    parser.add_argument("--out-counts", default="Data/host_tx_counts_filtered.tsv",
                        help="Filtered transcript counts")
    parser.add_argument("--out-annotations", default="Data/gene_annotations.tsv",
                        help="Gene annotation lookup table")
    parser.add_argument("--keep-types", nargs="+", default=None,
                        help="Gene types to keep (default: all matched genes)")
    parser.add_argument("--refresh", action="store_true", default=False,
                        help="Force re-query of mygene.info even if the "
                             "annotation cache exists")
    args = parser.parse_args()

    keep_types = set(args.keep_types) if args.keep_types else None

    print(f"Reading {args.input} ...")
    tx = pd.read_csv(args.input, sep="\t", index_col=0)
    all_genes = list(tx.index)
    print(f"  {len(all_genes)} genes, {tx.shape[1]} samples")

    orig_to_trimmed = {g: trim_version(g) for g in all_genes}
    unique_trimmed = sorted(set(orig_to_trimmed.values()))
    print(f"  {len(unique_trimmed)} unique gene names after version trimming")

    annot_map: dict[str, dict] = {}

    cache_exists = os.path.exists(args.out_annotations)
    if cache_exists and not args.refresh:
        print(f"\nAnnotation cache found: {args.out_annotations}")
        print("  Using cached annotations (offline mode); "
              "pass --refresh to re-query mygene.info.")
        cache = pd.read_csv(args.out_annotations, sep="\t")
        for _, row in cache.iterrows():
            annot_map[row["Gene_trimmed"]] = {
                "type": row["Gene_type"],
                "desc": row.get("Gene_description", "") or "",
            }
        print(f"  Loaded {len(annot_map)} cached gene annotations.")
    else:
        batches = [unique_trimmed[i:i + BATCH_SIZE]
                   for i in range(0, len(unique_trimmed), BATCH_SIZE)]
        print(f"\nQuerying mygene.info in {len(batches)} batches "
              f"(Ensembl / NCBI Entrez / HGNC via BioThings) ...")
        for idx, batch in enumerate(batches, 1):
            label = f"batch {idx}/{len(batches)}"
            print(f"  {label}  ({len(batch)} genes) ... ", end="", flush=True)
            result = query_mygene(batch, label)
            if result:
                annot_map.update(result)
                pc = sum(1 for v in result.values() if v["type"] == "protein-coding")
                print(f"OK {len(result)} matched ({pc} protein-coding)")
            else:
                print("FAILED")
            time.sleep(0.3)

    from collections import Counter
    type_counts = Counter(v["type"] for v in annot_map.values())
    print(f"\nMatched genes by type:")
    for gtype, cnt in type_counts.most_common():
        print(f"  {gtype:30s} : {cnt}")

    annot_rows = []
    for gene in all_genes:
        trimmed = orig_to_trimmed[gene]
        info = annot_map.get(trimmed)
        if info:
            annot_rows.append({
                "Gene": gene,
                "Gene_trimmed": trimmed,
                "Gene_type": info["type"],
                "Gene_description": info["desc"],
            })

    annot_df = pd.DataFrame(annot_rows)
    os.makedirs(os.path.dirname(args.out_annotations) or ".", exist_ok=True)
    annot_df.to_csv(args.out_annotations, sep="\t", index=False)
    print(f"\nAnnotation table ({len(annot_df)} genes) -> {args.out_annotations}")

    if keep_types is None:
        genes_to_keep = [row["Gene"] for row in annot_rows]
        filter_label = "all matched"
    else:
        genes_to_keep = [row["Gene"] for row in annot_rows
                         if row["Gene_type"] in keep_types]
        filter_label = ", ".join(sorted(keep_types))
    tx_filtered = tx.loc[genes_to_keep]

    os.makedirs(os.path.dirname(args.out_counts) or ".", exist_ok=True)
    tx_filtered.to_csv(args.out_counts, sep="\t")

    print(f"Kept {len(genes_to_keep)} / {len(all_genes)} genes ({filter_label})")
    print(f"Filtered counts -> {args.out_counts}")
    print("\nDone.")


if __name__ == "__main__":
    main()
