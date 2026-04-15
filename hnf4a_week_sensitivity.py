#!/usr/bin/env python3
"""
hnf4a_week_sensitivity.py

Test how robust the HNF4A enrichment signal is across different week-window
sizes (+-1 to +-52).  For each window, correlate Candida CLR vs host gene
expression (HD only, all matches), select positively correlated genes at
p < 0.05, run Fisher's exact enrichment against ENCODE/ChEA TF library,
and record the HNF4A ENCODE p-value.

Uses Kraken 2024-10 reference only.

Produces:
  - A plot of -log10(p-value) and -log10(adj. p-value) vs week window
  - A summary table with all values

Usage
-----
  python hnf4a_week_sensitivity.py
"""

import os
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact
from tqdm import tqdm

from _figutil import save_fig

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_DIR = "Data"
OUT_DIR  = "Results-HNF4A-sensitivity"

CLR_TABLE = os.path.join(DATA_DIR, "species_table_clr_202410.tsv")
METADATA  = os.path.join(DATA_DIR, "PRJNA398089_MGX_merged.tsv")
HMP2_META = os.path.join(DATA_DIR, "hmp2_metadata_2018-08-20.csv")
TX_COUNTS = os.path.join(DATA_DIR, "host_tx_counts_filtered.tsv")

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
LIBRARY      = "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
GMT_LOCAL    = os.path.join(DATA_DIR, f"{LIBRARY}.gmt")

TARGET           = "S__Candida albicans"
DIAGNOSES        = ["HD", "CD", "UC"]
CORR_P_THRESHOLD = 0.05
TF_OF_INTEREST   = "HNF4A ENCODE"

WEEK_RANGE = range(1, 53)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_transcriptomics() -> pd.DataFrame:
    """Build host-week averaged relative-abundance transcriptomics."""
    print("Loading host transcriptomics ...")
    hmp2_meta = pd.read_csv(HMP2_META, sep=",", low_memory=False)
    hmp2_meta["week_num"] = (pd.to_numeric(hmp2_meta["week_num"], errors="coerce")
                             .fillna(-999).astype(int))
    tx_meta = (hmp2_meta[hmp2_meta["data_type"] == "host_transcriptomics"]
               [["External ID", "Participant ID", "week_num", "diagnosis"]]
               .rename(columns={"External ID": "eid", "week_num": "week"})
               .copy())
    tx_meta["diagnosis"] = tx_meta["diagnosis"].replace({"nonIBD": "HD"})
    tx_meta["host_week"] = tx_meta["Participant ID"] + "__" + tx_meta["week"].astype(str)
    tx_meta = tx_meta.set_index("eid")

    tx = pd.read_csv(TX_COUNTS, sep="\t", index_col=0).T
    tx.index.name = "eid"
    tx = tx.merge(tx_meta[["host_week"]], left_index=True, right_index=True)
    tx = tx.groupby("host_week").mean()
    tx = tx.div(tx.sum(axis=1), axis=0)
    print(f"  {len(tx)} host-week rows, {tx.shape[1]} genes")
    return tx


def load_candida_clr() -> pd.DataFrame:
    """Extract Candida CLR per host_week, with diagnosis."""
    clr = pd.read_csv(CLR_TABLE, sep="\t").set_index("Run").replace(0, np.nan)
    mgx_meta = pd.read_csv(METADATA, sep="\t").set_index("Run")
    mgx_meta = mgx_meta[["diagnosis", "Participant ID", "week_num"]]
    mgx_meta["host_week"] = (mgx_meta["Participant ID"] + "__"
                             + mgx_meta["week_num"].astype(str))
    candida = (clr[[TARGET]]
               .merge(mgx_meta[["diagnosis", "host_week"]],
                      left_index=True, right_index=True))
    candida = (candida.groupby("host_week")
               .agg({TARGET: "mean", "diagnosis": "first"}))
    return candida


def load_gmt() -> dict[str, set[str]]:
    """Download (if needed) and parse GMT file."""
    if not os.path.exists(GMT_LOCAL):
        print(f"  Downloading {LIBRARY} ...")
        urllib.request.urlretrieve(ENRICHR_BASE + LIBRARY, GMT_LOCAL)
    gene_sets = {}
    with open(GMT_LOCAL) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if genes:
                gene_sets[name] = genes
    return gene_sets


def load_background() -> set[str]:
    """Load background gene names."""
    df = pd.read_csv(TX_COUNTS, sep="\t", usecols=[0])
    return set(df.iloc[:, 0])


# ---------------------------------------------------------------------------
# Per-window analysis
# ---------------------------------------------------------------------------

def correlate_and_enrich(candida: pd.DataFrame, tx: pd.DataFrame,
                         num_weeks: int, diagnosis: str, background: set[str],
                         gene_sets: dict[str, set[str]]) -> dict:
    """
    For a given week window and diagnosis:
    1. Window-match Candida CLR vs transcriptomics, all matches
    2. Linregress each gene
    3. Select positively correlated genes at p < threshold
    4. Fisher enrichment, return HNF4A ENCODE result
    """
    c2 = candida.copy()
    c2["host"] = c2.index.str.split("__").str[0]
    c2["week"] = c2.index.str.split("__").str[1].astype(int)

    t2 = tx.copy()
    t2["host"] = t2.index.str.split("__").str[0]
    t2["week"] = t2.index.str.split("__").str[1].astype(int)

    merged = c2.merge(t2, on="host", how="inner", suffixes=("_mgx", "_tx"))
    merged = merged[abs(merged["week_mgx"] - merged["week_tx"]) <= num_weeks]

    subset = merged[merged["diagnosis"] == diagnosis]
    subset = subset.drop(["diagnosis", "host", "week_mgx", "week_tx"], axis=1)
    subset = subset.apply(pd.to_numeric, errors="coerce")

    n_samples = len(subset)
    if n_samples < 2:
        return {"num_weeks": num_weeks, "n_samples": n_samples,
                "n_query": 0, "pval": np.nan, "adj_pval": np.nan,
                "odds_ratio": np.nan}

    query_genes = set()
    for col in subset.columns:
        if col == TARGET:
            continue
        x = pd.to_numeric(subset[col], errors="coerce")
        y = subset[TARGET]
        valid = x.notna() & y.notna()
        if valid.sum() < 2 or x[valid].nunique() < 2:
            continue
        slope, intercept, r, p, _ = stats.linregress(x[valid], y[valid])
        if p < CORR_P_THRESHOLD and r > 0:
            query_genes.add(col)

    if len(query_genes) < 5:
        return {"num_weeks": num_weeks, "n_samples": n_samples,
                "n_query": len(query_genes), "pval": np.nan,
                "adj_pval": np.nan, "odds_ratio": np.nan}

    N = len(background)
    k = len(query_genes & background)

    enrichment_results = []
    for name, gs in gene_sets.items():
        gs_bg = gs & background
        overlap = query_genes & gs_bg
        a = len(overlap)
        b = k - a
        c = len(gs_bg) - a
        d = N - a - b - c
        if a == 0:
            enrichment_results.append({"Name": name, "P_value": 1.0,
                                       "Odds_ratio": 0.0})
            continue
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        enrichment_results.append({"Name": name, "P_value": pval,
                                   "Odds_ratio": odds})

    edf = pd.DataFrame(enrichment_results).sort_values("P_value")

    n_tests = len(edf)
    ranks = np.arange(1, n_tests + 1)
    adj = np.minimum(1.0, edf["P_value"].values * n_tests / ranks)
    for i in range(n_tests - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    edf["Adjusted_p_value"] = adj

    hnf4a = edf[edf["Name"] == TF_OF_INTEREST]
    if len(hnf4a) == 0:
        return {"num_weeks": num_weeks, "n_samples": n_samples,
                "n_query": len(query_genes), "pval": np.nan,
                "adj_pval": np.nan, "odds_ratio": np.nan}

    row = hnf4a.iloc[0]
    return {"num_weeks": num_weeks, "n_samples": n_samples,
            "n_query": len(query_genes), "pval": row["P_value"],
            "adj_pval": row["Adjusted_p_value"],
            "odds_ratio": row["Odds_ratio"]}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=== HNF4A enrichment sensitivity to week window (Kraken 2024-10) ===\n")

    tx = load_transcriptomics()
    gene_sets = load_gmt()
    background = load_background()
    print(f"  {len(gene_sets)} TF gene sets, {len(background)} background genes\n")

    candida = load_candida_clr()
    for diag in DIAGNOSES:
        n = (candida["diagnosis"] == diag).sum()
        print(f"  {n} {diag} host-weeks with Candida data")

    all_results: dict[str, pd.DataFrame] = {}
    for diag in DIAGNOSES:
        results = []
        for w in tqdm(list(WEEK_RANGE), desc=f"  {diag}"):
            r = correlate_and_enrich(candida, tx, w, diag, background, gene_sets)
            results.append(r)
        df = pd.DataFrame(results)
        out_tsv = os.path.join(OUT_DIR, f"hnf4a_sensitivity_{diag}.tsv")
        df.to_csv(out_tsv, sep="\t", index=False)
        print(f"  {diag} saved -> {out_tsv}")
        all_results[diag] = df

    # Combined table
    combined = pd.DataFrame({"num_weeks": list(WEEK_RANGE)})
    for diag, df in all_results.items():
        combined[f"n_samples_{diag}"] = df["n_samples"].values
        combined[f"n_query_{diag}"]   = df["n_query"].values
        combined[f"pval_{diag}"]      = df["pval"].values
        combined[f"adj_pval_{diag}"]  = df["adj_pval"].values
        combined[f"odds_ratio_{diag}"] = df["odds_ratio"].values
    combined_path = os.path.join(OUT_DIR, "hnf4a_sensitivity_combined.tsv")
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"\nCombined -> {combined_path}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    diag_colors = {"HD": "#2171b5", "CD": "#cb181d", "UC": "#e6550d"}

    for diag, df in all_results.items():
        color = diag_colors[diag]
        neg_log_p = -np.log10(df["pval"].clip(lower=1e-300))
        ax1.plot(df["num_weeks"], neg_log_p,
                 "o-", color=color, markersize=4, linewidth=1.5, label=diag)
        neg_log_adj = -np.log10(df["adj_pval"].clip(lower=1e-300))
        ax2.plot(df["num_weeks"], neg_log_adj,
                 "s-", color=color, markersize=4, linewidth=1.5, label=diag)

    ax1.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8,
                label="p = 0.05")
    ax1.set_ylabel(r"$-\log_{10}(p\text{-value})$")
    ax1.set_title(f"{TF_OF_INTEREST} enrichment significance vs. week window")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8,
                label="adj. p = 0.05")
    ax2.set_xlabel("Week window (+-)")
    ax2.set_ylabel(r"$-\log_{10}(\mathrm{adj.}\ p\text{-value})$")
    ax2.set_title(f"{TF_OF_INTEREST} BH-adjusted significance vs. week window")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    plot_path = os.path.join(OUT_DIR, "hnf4a_sensitivity_plot.png")
    save_fig(fig, plot_path)
    plt.close(fig)
    print(f"Plot -> {plot_path} (+ .pdf)")

    print("\nDone.")


if __name__ == "__main__":
    main()
