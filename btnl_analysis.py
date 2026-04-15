#!/usr/bin/env python3
"""
btnl_analysis.py

Analyse the correlation between Candida albicans CLR and BTNL3/BTNL8 host
gene expression across diagnosis groups. Uses Kraken 2024-10 reference only.

1. For each diagnosis at +-4 weeks: scatter plots + correlation stats
2. Sensitivity analysis: -log10(p-value) vs week window (+-1 to +-52),
   one plot per gene, with curves for each diagnosis

Outputs in Results/BTNL/.

Usage
-----
  python btnl_analysis.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

from _figutil import save_fig

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_DIR = "Data"
OUT_DIR  = "Results/BTNL"

CLR_TABLE = os.path.join(DATA_DIR, "species_table_clr_202410.tsv")
METADATA  = os.path.join(DATA_DIR, "PRJNA398089_MGX_merged.tsv")
HMP2_META = os.path.join(DATA_DIR, "hmp2_metadata_2018-08-20.csv")
TX_COUNTS = os.path.join(DATA_DIR, "host_tx_counts_filtered.tsv")

TARGET     = "S__Candida albicans"
GENES      = ["BTNL3", "BTNL8"]
DIAGNOSES  = ["HD", "CD", "UC"]
DEFAULT_WEEKS = 4
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


# ---------------------------------------------------------------------------
# Window matching
# ---------------------------------------------------------------------------

def match_window(candida: pd.DataFrame, tx: pd.DataFrame,
                 diagnosis: str, num_weeks: int) -> pd.DataFrame:
    """Return matched subset for a given diagnosis and week window."""
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
    return subset


def gene_correlation(subset: pd.DataFrame, gene: str) -> dict:
    """Compute linregress for one gene vs TARGET in the matched subset."""
    if gene not in subset.columns or TARGET not in subset.columns:
        return {"r": np.nan, "slope": np.nan, "p": np.nan,
                "intercept": np.nan, "n": 0}

    x = pd.to_numeric(subset[gene], errors="coerce")
    y = subset[TARGET]
    valid = x.notna() & y.notna()

    if valid.sum() < 2 or x[valid].nunique() < 2:
        return {"r": np.nan, "slope": np.nan, "p": np.nan,
                "intercept": np.nan, "n": valid.sum()}

    slope, intercept, r, p, _ = stats.linregress(x[valid], y[valid])
    return {"r": r, "slope": slope, "p": p,
            "intercept": intercept, "n": valid.sum()}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(os.path.join(OUT_DIR, "Figures"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "Tables"), exist_ok=True)

    print("=== BTNL3/BTNL8 - Candida correlation analysis (Kraken 2024-10) ===\n")

    tx = load_transcriptomics()
    candida = load_candida_clr()
    print(f"  {len(candida)} host-weeks with Candida data\n")

    # ------------------------------------------------------------------
    # Part 1: Scatter plots at default window
    # ------------------------------------------------------------------
    print(f"--- Scatter plots (+-{DEFAULT_WEEKS} weeks) ---")
    summary_rows = []

    for diag in DIAGNOSES:
        subset = match_window(candida, tx, diag, DEFAULT_WEEKS)
        if len(subset) < 2:
            print(f"  {diag}: too few samples, skipping")
            continue

        for gene in GENES:
            res = gene_correlation(subset, gene)
            summary_rows.append({
                "Diagnosis": diag, "Gene": gene,
                "Correlation": res["r"], "Slope": res["slope"],
                "P_value": res["p"], "N": res["n"],
            })

            if np.isnan(res["r"]):
                continue

            x = pd.to_numeric(subset[gene], errors="coerce")
            y = subset[TARGET]
            valid = x.notna() & y.notna()

            fig, ax = plt.subplots(figsize=(6, 5))
            ax.scatter(x[valid], y[valid], alpha=0.5, color="steelblue", s=20)
            xfit = np.linspace(x[valid].min(), x[valid].max(), 100)
            ax.plot(xfit, res["slope"] * xfit + res["intercept"],
                    color="firebrick")
            ax.set_xlabel(gene)
            ax.set_ylabel(TARGET)
            ax.set_title(f"{gene} vs. Candida albicans\n"
                         f"({diag}, +-{DEFAULT_WEEKS}w, "
                         f"n={res['n']}, r={res['r']:.3f}, p={res['p']:.2e})")
            fig.tight_layout()
            fname = f"{gene}_vs_Candida_{diag}_pm{DEFAULT_WEEKS}w.png"
            save_fig(fig, os.path.join(OUT_DIR, "Figures", fname), dpi=150)
            plt.close(fig)

            pd.DataFrame({gene: x[valid], TARGET: y[valid]}).to_csv(
                os.path.join(OUT_DIR, "Tables",
                             f"{gene}_vs_Candida_{diag}_pm{DEFAULT_WEEKS}w.tsv"),
                sep="\t", index=False)

            print(f"  {diag}/{gene}: r={res['r']:.3f}, "
                  f"p={res['p']:.2e}, n={res['n']}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUT_DIR, "Tables", "btnl_correlation_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary -> {summary_path}")
    print(summary_df.to_string(index=False))
    print()

    # ------------------------------------------------------------------
    # Part 2: Sensitivity analysis
    # ------------------------------------------------------------------
    print("--- Sensitivity analysis ---")

    sensitivity = {gene: {} for gene in GENES}

    for diag in DIAGNOSES:
        for gene in GENES:
            sensitivity[gene][diag] = []

        for w in tqdm(list(WEEK_RANGE), desc=f"  {diag}"):
            subset = match_window(candida, tx, diag, w)
            for gene in GENES:
                res = gene_correlation(subset, gene)
                sensitivity[gene][diag].append({
                    "num_weeks": w, "r": res["r"],
                    "p": res["p"], "n": res["n"],
                })

    diag_colors = {"HD": "#2171b5", "CD": "#cb181d", "UC": "#e6550d"}

    for gene in GENES:
        combined = pd.DataFrame({"num_weeks": list(WEEK_RANGE)})
        for diag, rows in sensitivity[gene].items():
            df_tmp = pd.DataFrame(rows)
            combined[f"r_{diag}"] = df_tmp["r"].values
            combined[f"p_{diag}"] = df_tmp["p"].values
            combined[f"n_{diag}"] = df_tmp["n"].values
        tbl_path = os.path.join(OUT_DIR, "Tables",
                                f"{gene}_sensitivity_combined.tsv")
        combined.to_csv(tbl_path, sep="\t", index=False)
        print(f"  {gene} table -> {tbl_path}")

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

        for diag, rows in sensitivity[gene].items():
            df_tmp = pd.DataFrame(rows)
            neg_log_p = -np.log10(df_tmp["p"].clip(lower=1e-300))
            color = diag_colors[diag]

            ax1.plot(df_tmp["num_weeks"], neg_log_p,
                     color=color, marker="o", markersize=3, linewidth=1.2,
                     label=diag)

            ax2.plot(df_tmp["num_weeks"], df_tmp["r"],
                     color=color, marker="o", markersize=3, linewidth=1.2,
                     label=diag)

        ax1.axhline(-np.log10(0.05), color="gray", linestyle="--",
                    linewidth=0.8, label="p = 0.05")
        ax1.set_ylabel(r"$-\log_{10}(p\text{-value})$")
        ax1.set_title(f"{gene} - Candida albicans correlation significance "
                      f"vs. week window")
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)

        ax2.axhline(0, color="gray", linestyle="-", linewidth=0.5)
        ax2.set_xlabel("Week window (+-)")
        ax2.set_ylabel("Pearson r")
        ax2.set_title(f"{gene} - Candida albicans correlation coefficient "
                      f"vs. week window")
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        plot_path = os.path.join(OUT_DIR, "Figures",
                                 f"{gene}_sensitivity_plot.png")
        save_fig(fig, plot_path)
        plt.close(fig)
        print(f"  {gene} plot -> {plot_path} (+ .pdf)")

    print("\nDone.")


if __name__ == "__main__":
    main()
