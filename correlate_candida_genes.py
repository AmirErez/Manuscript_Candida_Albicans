#!/usr/bin/env python3
"""
correlate_candida_genes.py

Correlate Candida albicans CLR abundance with host gene expression,
separately for each diagnosis group (HD, CD, UC). Uses Kraken 2024-10
reference only.

Steps
-----
1. Build host-week averaged transcriptomics (relative-abundance normalised)
2. Extract Candida albicans CLR, add diagnosis, group by host_week
3. For each diagnosis group:
   a. Correlate: exact host_week match
   b. Correlate: +- NUM_WEEKS window (all matches per host by default)
   c. Compute p-values, slopes, scatter plots for genes of interest

Outputs are placed under Results/<diagnosis>/{Tables,Figures}.

Usage
-----
  python correlate_candida_genes.py
  python correlate_candida_genes.py --num-weeks 2
  python correlate_candida_genes.py --closest-only
  python correlate_candida_genes.py --diagnoses HD CD
"""

import argparse
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _figutil import save_fig
import pandas as pd
from scipy import stats
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_DIR = "Data"

CLR_TABLE = os.path.join(DATA_DIR, "species_table_clr_202410.tsv")
METADATA  = os.path.join(DATA_DIR, "PRJNA398089_MGX_merged.tsv")
HMP2_META = os.path.join(DATA_DIR, "hmp2_metadata_2018-08-20.csv")
TX_COUNTS = os.path.join(DATA_DIR, "host_tx_counts_filtered.tsv")

OUT_ROOT  = "Results"

TARGET = "S__Candida albicans"

DIAGNOSES = ["HD", "CD", "UC"]

GENES_OF_INTEREST = {
    "TNF", "TNFAIP2", "TNFAIP8", "IL10RA", "IL12RB1", "IL16",
    "IL17C", "IL17RA", "IL17RB", "IL17REL", "IL18R1", "IL1F10",
    "IL21", "IL21-AS1", "IL22RA2", "IL23A", "IL27RA", "IL2RA",
    "IL2RB", "IL2RG", "IL33", "IL36A", "IL4", "IL4I1", "IL6", "IL7R",
    "IL9RP3", "REL", "TLR1", "TLR10", "TLR2", "TLR4", "TLR5", "TLR6",
    "TLR7", "TLR8", "TLR8-AS1", "TLR9", "NOD2",
    "HNF4A", "HNF4G", "HNF4G-AS1",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_host_transcriptomics(out_dir: str) -> pd.DataFrame:
    """Build host-week averaged relative-abundance transcriptomics table."""
    print("Building host transcriptomics table ...")

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

    tx_path = os.path.join(out_dir, "Tables", "merged_host_transcriptomic_data.tsv")
    tx.to_csv(tx_path, sep="\t")
    print(f"  {len(tx)} host-week rows, {tx.shape[1]} genes -> {tx_path}")
    return tx


def build_candida_clr(out_dir: str) -> pd.DataFrame:
    """Extract Candida CLR per host_week, with diagnosis."""
    print("Extracting Candida CLR ...")

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

    candida_path = os.path.join(out_dir, "Tables", "candida_clr_withdiagnosis.tsv")
    candida.to_csv(candida_path, sep="\t")
    print(f"  {len(candida)} host-weeks with Candida data -> {candida_path}")
    return candida


def correlate_exact(candida: pd.DataFrame, tx: pd.DataFrame,
                    diagnosis: str, diag_dir: str) -> None:
    """Exact host_week merge correlations for one diagnosis group."""
    merged = candida.merge(tx, on="host_week", how="inner")
    subset = merged[merged["diagnosis"] == diagnosis].drop("diagnosis", axis=1)

    if len(subset) < 2:
        print(f"    Exact match: {len(subset)} samples, skipping")
        return

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        corrs = subset.corrwith(subset[TARGET], method="pearson")
    path = os.path.join(diag_dir, "Tables", "candida_gene_correlations_exact.tsv")
    corrs.to_csv(path, sep="\t")
    print(f"    Exact match: {len(subset)} samples -> {path}")


def correlate_window(candida: pd.DataFrame, tx: pd.DataFrame,
                     diagnosis: str, num_weeks: int, closest_only: bool,
                     diag_dir: str) -> None:
    """+-num_weeks window merge with full linregress and scatter plots."""

    c2 = candida.copy()
    c2["host"] = c2.index.str.split("__").str[0]
    c2["week"] = c2.index.str.split("__").str[1].astype(int)

    t2 = tx.copy()
    t2["host"] = t2.index.str.split("__").str[0]
    t2["week"] = t2.index.str.split("__").str[1].astype(int)

    merged = c2.merge(t2, on="host", how="inner", suffixes=("_mgx", "_tx"))
    merged = merged[abs(merged["week_mgx"] - merged["week_tx"]) <= num_weeks]

    if closest_only:
        merged["week_diff"] = abs(merged["week_mgx"] - merged["week_tx"])
        merged = (merged.sort_values("week_diff")
                  .drop_duplicates("host", keep="first")
                  .drop("week_diff", axis=1))

    subset = merged[merged["diagnosis"] == diagnosis]

    if len(subset) < 2:
        print(f"    +-{num_weeks}w: {len(subset)} samples, skipping")
        return

    suffix = f"pm{num_weeks}_weeks"
    meta_cols = subset[["diagnosis", "host", "week_mgx", "week_tx"]].copy()
    meta_cols.to_csv(
        os.path.join(diag_dir, "Tables",
                     f"candida_gene_correlations_{suffix}_matching_info.tsv"),
        sep="\t")

    subset = subset.drop(["diagnosis", "host", "week_mgx", "week_tx"], axis=1)
    subset = subset.apply(pd.to_numeric, errors="coerce")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        corrs = subset.corrwith(subset[TARGET], method="pearson")
    corrs.to_csv(
        os.path.join(diag_dir, "Tables",
                     f"candida_gene_correlations_{suffix}.tsv"),
        sep="\t")
    print(f"    +-{num_weeks}w: {len(subset)} matched pairs")

    fig_dir = os.path.join(diag_dir, "Figures")
    results = []

    for col in tqdm(subset.columns, desc=f"  linregress ({diagnosis})"):
        if col == TARGET:
            continue

        x = pd.to_numeric(subset[col], errors="coerce")
        y = subset[TARGET]
        valid = x.notna() & y.notna()

        if valid.sum() < 2 or x[valid].nunique() < 2:
            results.append({"Gene": col, "Correlation": np.nan,
                            "Slope": np.nan, "P_value": np.nan,
                            "N": valid.sum()})
            continue

        slope, intercept, r, p, _ = stats.linregress(x[valid], y[valid])

        if col.upper() in GENES_OF_INTEREST or col in GENES_OF_INTEREST:
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.scatter(x[valid], y[valid], alpha=0.5, color="steelblue", s=20)
            xfit = np.linspace(x[valid].min(), x[valid].max(), 100)
            ax.plot(xfit, slope * xfit + intercept, color="firebrick")
            ax.set_xlabel(col)
            ax.set_ylabel(TARGET)
            ax.set_title(f"{col} vs. Candida albicans "
                         f"(+-{num_weeks}w, {diagnosis}, n={valid.sum()})")
            fig.tight_layout()
            save_fig(fig,
                     os.path.join(fig_dir,
                                  f"{col}_vs_Candida_{suffix}.png"),
                     dpi=150)
            plt.close(fig)

            pd.DataFrame({col: x[valid], TARGET: y[valid]}).to_csv(
                os.path.join(diag_dir, "Tables",
                             f"{col}_vs_Candida_{suffix}.tsv"),
                sep="\t", index=False)

        results.append({"Gene": col, "Correlation": r,
                        "Slope": slope, "P_value": p, "N": valid.sum()})

    results_df = pd.DataFrame(results).sort_values("P_value")
    out_path = os.path.join(diag_dir, "Tables",
                            f"candida_gene_correlations_{suffix}_with_p-values.tsv")
    results_df.to_csv(out_path, sep="\t", index=False)

    print(f"    Top 10 correlated genes ({diagnosis}, +-{num_weeks}w):")
    print(results_df.dropna(subset=["P_value"]).head(10).to_string(index=False))
    print(f"    All results -> {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    global TARGET, TX_COUNTS

    parser = argparse.ArgumentParser(
        description="Correlate Candida albicans CLR with host gene expression.")
    parser.add_argument("--num-weeks", type=int, default=4,
                        help="Window size in weeks for temporal matching (default: 4)")
    parser.add_argument("--closest-only", action="store_true", default=False,
                        help="Keep only closest match per host (default: keep all)")
    parser.add_argument("--diagnoses", nargs="+", default=DIAGNOSES,
                        help=f"Diagnosis groups to analyse (default: {' '.join(DIAGNOSES)})")
    parser.add_argument("--target", default=TARGET,
                        help=f"Target species column (default: {TARGET})")
    parser.add_argument("--tx-counts", default=None,
                        help="Override path to transcript counts")
    args = parser.parse_args()

    TARGET = args.target
    if args.tx_counts:
        TX_COUNTS = args.tx_counts

    out_dir = OUT_ROOT

    os.makedirs(os.path.join(out_dir, "Tables"), exist_ok=True)
    for diag in args.diagnoses:
        os.makedirs(os.path.join(out_dir, diag, "Tables"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, diag, "Figures"), exist_ok=True)

    print(f"=== Candida-gene correlation (Kraken 2024-10) ===")
    print(f"CLR table:    {CLR_TABLE}")
    print(f"TX counts:    {TX_COUNTS}")
    print(f"Window:       +-{args.num_weeks} weeks")
    print(f"Match mode:   {'closest only' if args.closest_only else 'all within window'}")
    print(f"Diagnoses:    {', '.join(args.diagnoses)}")
    print(f"Output dir:   {out_dir}/")
    print()

    tx = build_host_transcriptomics(out_dir)
    candida = build_candida_clr(out_dir)

    for diag in args.diagnoses:
        diag_dir = os.path.join(out_dir, diag)
        n_samples = (candida["diagnosis"] == diag).sum()
        print(f"\n--- {diag} ({n_samples} host-weeks) ---")

        correlate_exact(candida, tx, diag, diag_dir)
        correlate_window(candida, tx, diag, args.num_weeks,
                         args.closest_only, diag_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
