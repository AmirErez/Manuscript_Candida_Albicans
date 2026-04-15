#!/usr/bin/env python3
"""
hnf4a_btnl_fungal_ranking.py

For each of the 13 curated prevalent human-associated fungal species, compute:
  1. HNF4A ENCODE enrichment p-value (Fisher's exact on positively correlated
     host genes, +-4 weeks, HD)
  2. Direct linregress p-value for BTNL3 vs taxon CLR
  3. Direct linregress p-value for BTNL8 vs taxon CLR

Then combine the three p-values using the geometric mean and rank the taxa
by the combined score.

Output
------
  Results-HNF4A-fungi/hnf4a_btnl_fungal_ranking.tsv
  Results-HNF4A-fungi/hnf4a_btnl_fungal_ranking_plot.png
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
OUT_DIR  = "Results-HNF4A-fungi"

CLR_TABLE = os.path.join(DATA_DIR, "species_table_clr_202410.tsv")
METADATA  = os.path.join(DATA_DIR, "PRJNA398089_MGX_merged.tsv")
HMP2_META = os.path.join(DATA_DIR, "hmp2_metadata_2018-08-20.csv")
TX_COUNTS = os.path.join(DATA_DIR, "host_tx_counts_filtered.tsv")

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
LIBRARY      = "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
GMT_LOCAL    = os.path.join(DATA_DIR, f"{LIBRARY}.gmt")

DIAGNOSIS        = "HD"
NUM_WEEKS        = 4
CORR_P_THRESHOLD = 0.05
TF_OF_INTEREST   = "HNF4A ENCODE"
BTNL_GENES       = ["BTNL3", "BTNL8"]
HIGHLIGHT_TAXON  = "S__Candida albicans"
PVAL_CAP         = 1e-10  # cap used for the "capped" combined score

FUNGAL_TAXA = [
    "S__Malassezia restricta",
    "S__Saccharomyces cerevisiae",
    "S__Candida albicans",
    "S__Candida dubliniensis",
    "S__Aspergillus fumigatus",
    "S__Cryptococcus neoformans",
    "S__Lodderomyces elongisporus",
    "S__Malassezia japonica",
    "S__Nakaseomyces glabratus",
    "S__Debaryomyces hansenii",
    "S__Candida orthopsilosis",
    "S__Pichia kudriavzevii",
    "S__Lodderomyces beijingensis",
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_transcriptomics() -> pd.DataFrame:
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


def load_species_clr() -> pd.DataFrame:
    clr = pd.read_csv(CLR_TABLE, sep="\t").set_index("Run").replace(0, np.nan)
    mgx_meta = pd.read_csv(METADATA, sep="\t").set_index("Run")
    mgx_meta = mgx_meta[["diagnosis", "Participant ID", "week_num"]]
    mgx_meta["host_week"] = (mgx_meta["Participant ID"] + "__"
                             + mgx_meta["week_num"].astype(str))
    merged = clr.merge(mgx_meta[["diagnosis", "host_week"]],
                       left_index=True, right_index=True)
    diag = merged.groupby("host_week")["diagnosis"].first()
    species = (merged.drop(columns=["diagnosis"])
               .groupby("host_week").mean())
    species["diagnosis"] = diag
    return species


def load_gmt() -> dict[str, set[str]]:
    if not os.path.exists(GMT_LOCAL):
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
    df = pd.read_csv(TX_COUNTS, sep="\t", usecols=[0])
    return set(df.iloc[:, 0])


# ---------------------------------------------------------------------------
# Per-taxon analysis
# ---------------------------------------------------------------------------

def analyse_taxon(taxon: str, species: pd.DataFrame, tx: pd.DataFrame,
                  background: set[str],
                  gene_sets: dict[str, set[str]]) -> dict:
    c2 = species[[taxon, "diagnosis"]].copy().rename(columns={taxon: "TAXON"})
    c2 = c2.dropna(subset=["TAXON"])
    c2["host"] = c2.index.str.split("__").str[0]
    c2["week"] = c2.index.str.split("__").str[1].astype(int)

    t2 = tx.copy()
    t2["host"] = t2.index.str.split("__").str[0]
    t2["week"] = t2.index.str.split("__").str[1].astype(int)

    merged = c2.merge(t2, on="host", how="inner", suffixes=("_mgx", "_tx"))
    merged = merged[abs(merged["week_mgx"] - merged["week_tx"]) <= NUM_WEEKS]
    subset = merged[merged["diagnosis"] == DIAGNOSIS]
    subset = subset.drop(["diagnosis", "host", "week_mgx", "week_tx"], axis=1)
    subset = subset.apply(pd.to_numeric, errors="coerce")

    n_samples = len(subset)
    out = {"Taxon": taxon, "n_samples": n_samples, "n_query": 0,
           "HNF4A_pval": np.nan, "HNF4A_odds_ratio": np.nan,
           "HNF4A_overlap": 0,
           "BTNL3_r": np.nan, "BTNL3_pval": np.nan,
           "BTNL8_r": np.nan, "BTNL8_pval": np.nan}
    if n_samples < 3:
        return out

    y = subset["TAXON"]

    # Direct BTNL3 / BTNL8 linregress
    for gene in BTNL_GENES:
        if gene not in subset.columns:
            continue
        x = pd.to_numeric(subset[gene], errors="coerce")
        valid = x.notna() & y.notna()
        if valid.sum() < 3 or x[valid].nunique() < 2:
            continue
        slope, intercept, r, p, _ = stats.linregress(x[valid], y[valid])
        out[f"{gene}_r"] = r
        out[f"{gene}_pval"] = p

    # HNF4A ENCODE enrichment
    query_genes = set()
    for col in subset.columns:
        if col == "TAXON":
            continue
        x = pd.to_numeric(subset[col], errors="coerce")
        valid = x.notna() & y.notna()
        if valid.sum() < 3 or x[valid].nunique() < 2:
            continue
        slope, intercept, r, p, _ = stats.linregress(x[valid], y[valid])
        if p < CORR_P_THRESHOLD and r > 0:
            query_genes.add(col)

    out["n_query"] = len(query_genes)
    if len(query_genes) < 5:
        return out

    N = len(background)
    k = len(query_genes & background)

    enrich_rows = []
    for name, gs in gene_sets.items():
        gs_bg = gs & background
        overlap = query_genes & gs_bg
        a = len(overlap)
        b = k - a
        c = len(gs_bg) - a
        d = N - a - b - c
        if a == 0:
            enrich_rows.append({"Name": name, "P_value": 1.0,
                                "Odds_ratio": 0.0, "Overlap": 0})
            continue
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        enrich_rows.append({"Name": name, "P_value": pval,
                            "Odds_ratio": odds, "Overlap": a})

    edf = pd.DataFrame(enrich_rows)
    hit = edf[edf["Name"] == TF_OF_INTEREST]
    if len(hit) == 0:
        return out
    row = hit.iloc[0]
    out["HNF4A_pval"] = row["P_value"]
    out["HNF4A_odds_ratio"] = row["Odds_ratio"]
    out["HNF4A_overlap"] = int(row["Overlap"])
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"=== HNF4A + BTNL3/BTNL8 per fungal species "
          f"({DIAGNOSIS}, +-{NUM_WEEKS}w, Kraken 2024-10) ===\n")

    tx = load_transcriptomics()
    gene_sets = load_gmt()
    background = load_background()
    print(f"  {len(gene_sets)} TF gene sets, {len(background)} background genes\n")

    species = load_species_clr()

    clr_raw = pd.read_csv(CLR_TABLE, sep="\t").set_index("Run").replace(0, np.nan)
    n_total = len(clr_raw)
    prevalence = clr_raw.notna().sum(axis=0) / n_total

    taxa = [t for t in FUNGAL_TAXA if t in species.columns]
    print(f"  Analysing {len(taxa)} fungal species\n")

    results = []
    for taxon in tqdm(taxa, desc="  taxa"):
        r = analyse_taxon(taxon, species, tx, background, gene_sets)
        r["prevalence"] = prevalence[taxon]
        results.append(r)

    df = pd.DataFrame(results)

    # ------------------------------------------------------------------
    # Combined score: geometric mean of HNF4A, BTNL3, BTNL8 p-values.
    # Two variants: uncapped (floor 1e-300) and capped at PVAL_CAP so that
    # a single extreme HNF4A hit cannot dominate.
    # ------------------------------------------------------------------
    floor = 1e-300
    p_cols = ["HNF4A_pval", "BTNL3_pval", "BTNL8_pval"]

    log_p = np.log(df[p_cols].clip(lower=floor))
    df["combined_pval_geomean"] = np.exp(log_p.mean(axis=1))
    df["neg_log10_combined"] = -np.log10(df["combined_pval_geomean"].clip(lower=floor))

    log_p_cap = np.log(df[p_cols].clip(lower=PVAL_CAP, upper=1.0))
    df["combined_pval_geomean_capped"] = np.exp(log_p_cap.mean(axis=1))
    df["neg_log10_combined_capped"] = -np.log10(
        df["combined_pval_geomean_capped"].clip(lower=floor))

    # ------------------------------------------------------------------
    # Write uncapped table + plot
    # ------------------------------------------------------------------
    df_unc = df.sort_values("combined_pval_geomean")
    unc_tsv = os.path.join(OUT_DIR, "hnf4a_btnl_fungal_ranking.tsv")
    df_unc.to_csv(unc_tsv, sep="\t", index=False)
    print(f"\nUncapped table -> {unc_tsv}\n")

    show_cols = ["Taxon", "prevalence", "n_samples",
                 "HNF4A_pval", "BTNL3_r", "BTNL3_pval",
                 "BTNL8_r", "BTNL8_pval", "combined_pval_geomean"]
    disp = df_unc[show_cols].copy()
    disp["prevalence"] = disp["prevalence"].map(lambda v: f"{v:.3f}")
    for c in ["HNF4A_pval", "BTNL3_pval", "BTNL8_pval", "combined_pval_geomean"]:
        disp[c] = disp[c].map(lambda v: f"{v:.3e}" if pd.notna(v) else "NA")
    for c in ["BTNL3_r", "BTNL8_r"]:
        disp[c] = disp[c].map(lambda v: f"{v:+.3f}" if pd.notna(v) else "NA")
    print("Uncapped ranking:")
    print(disp.to_string(index=False))

    fig, ax = plt.subplots(figsize=(10, max(5, 0.45 * len(df_unc))))
    labels = [t.replace("S__", "") for t in df_unc["Taxon"]]
    colors = ["#cb181d" if t == HIGHLIGHT_TAXON else "#2171b5"
              for t in df_unc["Taxon"]]
    y_pos = np.arange(len(df_unc))
    ax.barh(y_pos, df_unc["neg_log10_combined"].values, color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--",
               linewidth=0.8, label="p = 0.05")
    ax.set_xlabel(r"$-\log_{10}$(geometric mean of HNF4A, BTNL3, BTNL8 p-values)")
    ax.set_title(f"Combined HNF4A + BTNL3 + BTNL8 ranking per fungal species "
                 f"({DIAGNOSIS}, +-{NUM_WEEKS}w)")
    ax.legend(loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    unc_png = os.path.join(OUT_DIR, "hnf4a_btnl_fungal_ranking_plot.png")
    save_fig(fig, unc_png)
    plt.close(fig)
    print(f"Uncapped plot -> {unc_png} (+ .pdf)")

    # ------------------------------------------------------------------
    # Write capped table + plot
    # ------------------------------------------------------------------
    df_cap = df.sort_values("combined_pval_geomean_capped")
    cap_tsv = os.path.join(OUT_DIR, "hnf4a_btnl_fungal_ranking_capped1e-10.tsv")
    df_cap.to_csv(cap_tsv, sep="\t", index=False)
    print(f"\nCapped table -> {cap_tsv}")

    fig, ax = plt.subplots(figsize=(10, max(5, 0.45 * len(df_cap))))
    labels = [t.replace("S__", "") for t in df_cap["Taxon"]]
    colors = ["#cb181d" if t == HIGHLIGHT_TAXON else "#2171b5"
              for t in df_cap["Taxon"]]
    y_pos = np.arange(len(df_cap))
    ax.barh(y_pos, df_cap["neg_log10_combined_capped"].values, color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--",
               linewidth=0.8, label="p = 0.05")
    ax.set_xlabel(r"$-\log_{10}$(geometric mean, HNF4A/BTNL3/BTNL8 p-values capped at $10^{-10}$)")
    ax.set_title(f"Combined HNF4A + BTNL3 + BTNL8 ranking, capped at p=1e-10\n"
                 f"(13 curated prevalent human-associated fungi, "
                 f"{DIAGNOSIS}, +-{NUM_WEEKS}w)")
    ax.legend(loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    cap_png = os.path.join(OUT_DIR, "hnf4a_btnl_fungal_ranking_capped1e-10_plot.png")
    save_fig(fig, cap_png)
    plt.close(fig)
    print(f"Capped plot -> {cap_png} (+ .pdf)")

    print("\nDone.")


if __name__ == "__main__":
    main()
