#!/usr/bin/env python3
"""
enrichr_enrichment.py

Local Fisher's exact test enrichment against any Enrichr gene set library.

Downloads the GMT file from Enrichr on first run (cached locally in Data/).
Runs enrichment for each diagnosis, using genes with p < threshold and
selected correlation sign as the query set, and the full filtered gene list
(~30K) as background.

Outputs a table per condition with: Name, p-value, adjusted p-value (BH),
odds ratio, combined score, and overlapping genes.

Usage
-----
  python enrichr_enrichment.py --library "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" \
      --diagnoses HD --corr-sign pos --p-threshold 0.05

  python enrichr_enrichment.py --library "GO_Biological_Process_2023" \
      --corr-sign neg --p-threshold 0.05
"""

import argparse
import os
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_DIR = "Data"
RESULTS_ROOT = "Results"

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="

BACKGROUND_GENES_FILE = os.path.join(DATA_DIR, "host_tx_counts_filtered.tsv")

DIAGNOSES = ["HD", "CD", "UC"]

P_THRESHOLD = 0.05
CORR_SIGN   = "pos"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def download_gmt(library: str) -> str:
    """Download GMT file if not already cached. Returns local path."""
    local_path = os.path.join(DATA_DIR, f"{library}.gmt")
    if os.path.exists(local_path):
        print(f"  GMT already cached: {local_path}")
        return local_path
    url = ENRICHR_BASE + library
    print(f"  Downloading {library} from Enrichr ...")
    urllib.request.urlretrieve(url, local_path)
    print(f"  Saved to {local_path}")
    return local_path


def parse_gmt(path: str) -> dict[str, set[str]]:
    """Parse GMT file into {set_name: {gene1, gene2, ...}}."""
    gene_sets = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if genes:
                gene_sets[name] = genes
    return gene_sets


def load_background_genes(path: str) -> set[str]:
    """Load gene names from filtered transcript counts (genes are rows)."""
    df = pd.read_csv(path, sep="\t", usecols=[0])
    return set(df.iloc[:, 0])


def load_query_genes(results_path: str, p_thresh: float,
                     corr_sign: str) -> set[str]:
    """Load significant genes from correlation results."""
    df = pd.read_csv(results_path, sep="\t")
    mask = df["P_value"] < p_thresh
    if corr_sign == "pos":
        mask &= df["Correlation"] > 0
    elif corr_sign == "neg":
        mask &= df["Correlation"] < 0
    return set(df.loc[mask, "Gene"])


def enrichment_test(query: set[str], background: set[str],
                    gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    """
    Fisher's exact test for each gene set.

    For each gene set:
      a = query AND set  (overlap)
      b = query NOT set
      c = set NOT query  (in background)
      d = background NOT query NOT set
    """
    N = len(background)
    k = len(query & background)

    results = []
    for name, gs in gene_sets.items():
        gs_bg = gs & background
        overlap = query & gs_bg
        a = len(overlap)
        b = k - a
        c = len(gs_bg) - a
        d = N - a - b - c

        if a == 0:
            results.append({
                "Name": name, "P_value": 1.0, "Odds_ratio": 0.0,
                "Combined_score": 0.0, "Overlap_count": 0,
                "Set_size": len(gs_bg), "Overlap_genes": "",
            })
            continue

        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append({
            "Name": name, "P_value": pval, "Odds_ratio": odds,
            "Combined_score": 0.0, "Overlap_count": a,
            "Set_size": len(gs_bg), "Overlap_genes": ";".join(sorted(overlap)),
        })

    df = pd.DataFrame(results).sort_values("P_value")

    n_tests = len(df)
    ranks = np.arange(1, n_tests + 1)
    pvals = df["P_value"].values
    adj = np.minimum(1.0, pvals * n_tests / ranks)
    for i in range(n_tests - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    df["Adjusted_p_value"] = adj

    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = -np.log(np.clip(df["P_value"].values, 1e-300, 1.0))
        log_or = np.log(np.clip(df["Odds_ratio"].values, 1e-300, None))
        df["Combined_score"] = log_p * log_or
    df["Combined_score"] = df["Combined_score"].replace([np.inf, -np.inf], 0.0)

    df = df[["Name", "P_value", "Adjusted_p_value", "Odds_ratio",
             "Combined_score", "Overlap_count", "Set_size", "Overlap_genes"]]
    return df.sort_values("P_value").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Local enrichment (any Enrichr library) on Candida-correlated genes.")
    parser.add_argument("--library", required=True,
                        help="Enrichr library name (e.g. GO_Biological_Process_2023)")
    parser.add_argument("--output-tag", default=None,
                        help="Short tag for output filename (default: derived from library)")
    parser.add_argument("--diagnoses", nargs="+", default=DIAGNOSES,
                        help=f"Diagnosis groups (default: {' '.join(DIAGNOSES)})")
    parser.add_argument("--p-threshold", type=float, default=P_THRESHOLD,
                        help=f"P-value threshold for query genes (default: {P_THRESHOLD})")
    parser.add_argument("--corr-sign", choices=["pos", "neg", "any"], default=CORR_SIGN,
                        help=f"Correlation direction filter (default: {CORR_SIGN})")
    parser.add_argument("--num-weeks", type=int, default=4,
                        help="Which pm*_weeks result file to use (default: 4)")
    parser.add_argument("--background", default=BACKGROUND_GENES_FILE,
                        help="Background gene list file")
    args = parser.parse_args()

    tag = args.output_tag or args.library

    print(f"Loading library: {args.library} ...")
    os.makedirs(DATA_DIR, exist_ok=True)
    gmt_path = download_gmt(args.library)
    gene_sets = parse_gmt(gmt_path)
    print(f"  {len(gene_sets)} gene sets loaded")

    print(f"Loading background genes from {args.background} ...")
    background = load_background_genes(args.background)
    print(f"  {len(background)} background genes")

    suffix = f"pm{args.num_weeks}_weeks"
    for diag in args.diagnoses:
        results_file = os.path.join(
            RESULTS_ROOT, diag, "Tables",
            f"candida_gene_correlations_{suffix}_with_p-values.tsv")

        if not os.path.exists(results_file):
            print(f"\n[{diag}] {results_file} not found, skipping")
            continue

        query = load_query_genes(results_file, args.p_threshold, args.corr_sign)
        print(f"\n[{diag}] {len(query)} query genes "
              f"(p<{args.p_threshold}, corr {args.corr_sign})")

        if len(query) < 5:
            print(f"  Too few query genes, skipping")
            continue

        enrichment = enrichment_test(query, background, gene_sets)

        out_path = os.path.join(
            RESULTS_ROOT, diag, "Tables",
            f"enrichr_{tag}_{suffix}.tsv")
        enrichment.to_csv(out_path, sep="\t", index=False)

        sig = enrichment[enrichment["Adjusted_p_value"] < 0.05]
        print(f"  {len(sig)} terms significant (adj. p < 0.05)")
        if len(sig) > 0:
            print(sig[["Name", "P_value", "Adjusted_p_value",
                       "Odds_ratio", "Combined_score"]].head(15)
                  .to_string(index=False))
        print(f"  Results -> {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
