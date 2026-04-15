#!/usr/bin/env python3
"""
hnf4a_btnl_taxa_ranking.py

For the top-N most abundant species (by HD mean CLR) plus Candida albicans,
compute:
  1. HNF4A ENCODE enrichment p-value (Fisher's exact on positively correlated
     host genes, +-4 weeks, HD)
  2. Direct linregress p-value for BTNL3 vs taxon CLR
  3. Direct linregress p-value for BTNL8 vs taxon CLR

Combine the three p-values using the geometric mean and rank taxa by the
combined score.

This version is vectorised: per-gene Pearson correlations are computed with
a single numpy matmul per taxon rather than 33k linregress calls.

Output
------
  Results-HNF4A-taxa/hnf4a_btnl_taxa_ranking_top{N}.tsv
  Results-HNF4A-taxa/hnf4a_btnl_taxa_ranking_top{N}_plot.png
"""

import os
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, t as tdist
from tqdm import tqdm

from _figutil import save_fig

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_DIR = "Data"
OUT_DIR  = "Results-HNF4A-taxa"

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
TOP_N            = 100
FORCED_TAXA      = ["S__Candida albicans"]
HIGHLIGHT_TAXON  = "S__Candida albicans"
PVAL_CAP         = 1e-10  # cap used for the "capped" combined score


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
# Vectorised per-taxon analysis
# ---------------------------------------------------------------------------

def build_pair_index(species_hw: pd.Index, tx_hw: pd.Index, num_weeks: int):
    """
    Build the list of (species_row_idx, tx_row_idx) pairs that satisfy
    same host & |week difference| <= num_weeks. Returns numpy arrays.
    """
    sp_host = species_hw.str.split("__").str[0].values
    sp_week = species_hw.str.split("__").str[1].astype(int).values
    tx_host = tx_hw.str.split("__").str[0].values
    tx_week = tx_hw.str.split("__").str[1].astype(int).values

    # Group tx rows by host
    tx_by_host: dict[str, list[tuple[int, int]]] = {}
    for i, h in enumerate(tx_host):
        tx_by_host.setdefault(h, []).append((tx_week[i], i))

    sp_pairs: list[int] = []
    tx_pairs: list[int] = []
    for si, (h, w) in enumerate(zip(sp_host, sp_week)):
        if h not in tx_by_host:
            continue
        for tw, ti in tx_by_host[h]:
            if abs(tw - w) <= num_weeks:
                sp_pairs.append(si)
                tx_pairs.append(ti)
    return np.asarray(sp_pairs, dtype=np.int64), np.asarray(tx_pairs, dtype=np.int64)


def vectorised_pearson(y: np.ndarray, X: np.ndarray):
    """
    y: (n,), X: (n, g). Returns (r, p) each of shape (g,).
    Pearson correlation + two-sided t-test p-value.
    """
    n = y.shape[0]
    yc = y - y.mean()
    Xm = X.mean(axis=0)
    Xc = X - Xm  # (n, g)
    y_ss = (yc * yc).sum()
    X_ss = (Xc * Xc).sum(axis=0)
    denom = np.sqrt(X_ss * y_ss)
    num = Xc.T @ yc  # (g,)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0, num / denom, np.nan)
        dof = n - 2
        r2 = np.clip(r * r, 0.0, 1.0 - 1e-12)
        t = r * np.sqrt(dof / (1.0 - r2))
        p = 2.0 * tdist.sf(np.abs(t), dof)
    return r, p


def analyse_taxon(taxon: str,
                  species_hd: pd.DataFrame,
                  tx_np: np.ndarray,
                  gene_names: np.ndarray,
                  btnl_idx: dict[str, int],
                  pair_cache: dict,
                  background: set[str],
                  gene_sets: dict[str, set[str]]) -> dict:
    out = {"Taxon": taxon, "n_samples": 0, "n_query": 0,
           "HNF4A_pval": np.nan, "HNF4A_odds_ratio": np.nan,
           "HNF4A_overlap": 0,
           "BTNL3_r": np.nan, "BTNL3_pval": np.nan,
           "BTNL8_r": np.nan, "BTNL8_pval": np.nan}

    col = species_hd[taxon]
    valid_rows_mask = col.notna().values
    if valid_rows_mask.sum() < 3:
        return out

    sp_pairs_all, tx_pairs_all = pair_cache["pairs"]
    # Restrict to species rows with non-NaN taxon value
    sel = valid_rows_mask[sp_pairs_all]
    sp_pairs = sp_pairs_all[sel]
    tx_pairs = tx_pairs_all[sel]
    n_samples = sp_pairs.size
    out["n_samples"] = n_samples
    if n_samples < 3:
        return out

    y = col.values[sp_pairs].astype(np.float64)
    X = tx_np[tx_pairs]  # (n_samples, n_genes)  view into float32 array

    # Cast to float64 for stable correlation math
    X64 = X.astype(np.float64, copy=False)
    r, p = vectorised_pearson(y, X64)

    # BTNL
    for gene in BTNL_GENES:
        j = btnl_idx.get(gene)
        if j is None:
            continue
        out[f"{gene}_r"] = float(r[j]) if np.isfinite(r[j]) else np.nan
        out[f"{gene}_pval"] = float(p[j]) if np.isfinite(p[j]) else np.nan

    # HNF4A enrichment: positively correlated genes at p < threshold
    sel_mask = np.isfinite(p) & (p < CORR_P_THRESHOLD) & (r > 0)
    query_genes = set(gene_names[sel_mask].tolist())
    out["n_query"] = len(query_genes)
    if len(query_genes) < 5:
        return out

    N = len(background)
    k = len(query_genes & background)

    hnf4a_gs = gene_sets.get(TF_OF_INTEREST)
    if hnf4a_gs is None:
        return out
    gs_bg = hnf4a_gs & background
    overlap = query_genes & gs_bg
    a = len(overlap)
    b = k - a
    c = len(gs_bg) - a
    d = N - a - b - c
    if a == 0:
        out["HNF4A_pval"] = 1.0
        out["HNF4A_odds_ratio"] = 0.0
        out["HNF4A_overlap"] = 0
        return out
    odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
    out["HNF4A_pval"] = pval
    out["HNF4A_odds_ratio"] = odds
    out["HNF4A_overlap"] = int(a)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"=== HNF4A + BTNL3/BTNL8 per top-{TOP_N} taxon "
          f"({DIAGNOSIS}, +-{NUM_WEEKS}w, Kraken 2024-10) ===\n")

    tx = load_transcriptomics()
    gene_sets = load_gmt()
    background = load_background()
    print(f"  {len(gene_sets)} TF gene sets, {len(background)} background genes\n")

    species = load_species_clr()
    species_hd = species[species["diagnosis"] == DIAGNOSIS].drop(columns=["diagnosis"])

    # Prebuild tx matrix (float32) & gene name array
    gene_names = tx.columns.values.astype(object)
    tx_np = tx.values.astype(np.float32, copy=False)
    btnl_idx = {g: int(np.where(gene_names == g)[0][0])
                for g in BTNL_GENES if g in tx.columns}

    # Prebuild species/tx window-match pair index (computed once, reused per taxon)
    sp_pairs, tx_pairs = build_pair_index(species_hd.index, tx.index, NUM_WEEKS)
    print(f"  {len(sp_pairs)} HD species-tx (|dw|<={NUM_WEEKS}) pair candidates")
    pair_cache = {"pairs": (sp_pairs, tx_pairs)}

    # Abundance ranking
    mean_clr = species_hd.mean(axis=0, skipna=True).sort_values(ascending=False)
    top_taxa = mean_clr.head(TOP_N).index.tolist()
    for forced in FORCED_TAXA:
        if forced in mean_clr.index and forced not in top_taxa:
            rank = list(mean_clr.index).index(forced) + 1
            top_taxa.append(forced)
            print(f"  Added forced taxon {forced} (HD mean-CLR rank {rank})")
    print(f"  Analysing {len(top_taxa)} taxa "
          f"(top {TOP_N} by HD mean CLR + forced)\n")

    results = []
    for taxon in tqdm(top_taxa, desc="  taxa"):
        r = analyse_taxon(taxon, species_hd, tx_np, gene_names, btnl_idx,
                          pair_cache, background, gene_sets)
        r["mean_clr"] = mean_clr[taxon]
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

    show_cols = ["Taxon", "mean_clr", "n_samples",
                 "HNF4A_pval", "BTNL3_r", "BTNL3_pval",
                 "BTNL8_r", "BTNL8_pval", "combined_pval_geomean"]

    def _format(d):
        disp = d[show_cols].copy()
        disp["mean_clr"] = disp["mean_clr"].map(lambda v: f"{v:+.3f}")
        for c in ["HNF4A_pval", "BTNL3_pval", "BTNL8_pval", "combined_pval_geomean"]:
            disp[c] = disp[c].map(lambda v: f"{v:.3e}" if pd.notna(v) else "NA")
        for c in ["BTNL3_r", "BTNL8_r"]:
            disp[c] = disp[c].map(lambda v: f"{v:+.3f}" if pd.notna(v) else "NA")
        return disp

    # ------------------------------------------------------------------
    # Uncapped table + plot
    # ------------------------------------------------------------------
    df_unc = df.sort_values("combined_pval_geomean")
    unc_tsv = os.path.join(OUT_DIR, f"hnf4a_btnl_taxa_ranking_top{TOP_N}.tsv")
    df_unc.to_csv(unc_tsv, sep="\t", index=False)
    print(f"\nUncapped table -> {unc_tsv}")
    print("\nUncapped top 30:")
    print(_format(df_unc).head(30).to_string(index=False))

    fig, ax = plt.subplots(figsize=(11, max(6, 0.22 * len(df_unc))))
    labels = [t.replace("S__", "") for t in df_unc["Taxon"]]
    colors = ["#cb181d" if t == HIGHLIGHT_TAXON else "#2171b5"
              for t in df_unc["Taxon"]]
    y_pos = np.arange(len(df_unc))
    ax.barh(y_pos, df_unc["neg_log10_combined"].values, color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--",
               linewidth=0.8, label="p = 0.05")
    ax.set_xlabel(r"$-\log_{10}$(geometric mean of HNF4A, BTNL3, BTNL8 p-values)")
    ax.set_title(f"Combined HNF4A + BTNL3 + BTNL8 ranking "
                 f"(top {TOP_N} taxa + Candida, {DIAGNOSIS}, +-{NUM_WEEKS}w)")
    ax.legend(loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    unc_png = os.path.join(OUT_DIR, f"hnf4a_btnl_taxa_ranking_top{TOP_N}_plot.png")
    save_fig(fig, unc_png)
    plt.close(fig)
    print(f"\nUncapped plot -> {unc_png} (+ .pdf)")

    # ------------------------------------------------------------------
    # Capped table + plot
    # ------------------------------------------------------------------
    df_cap = df.sort_values("combined_pval_geomean_capped")
    cap_tsv = os.path.join(
        OUT_DIR, f"hnf4a_btnl_taxa_ranking_top{TOP_N}_capped1e-10.tsv")
    df_cap.to_csv(cap_tsv, sep="\t", index=False)
    print(f"Capped table -> {cap_tsv}")

    fig, ax = plt.subplots(figsize=(11, max(6, 0.22 * len(df_cap))))
    labels = [t.replace("S__", "") for t in df_cap["Taxon"]]
    colors = ["#cb181d" if t == HIGHLIGHT_TAXON else "#2171b5"
              for t in df_cap["Taxon"]]
    y_pos = np.arange(len(df_cap))
    ax.barh(y_pos, df_cap["neg_log10_combined_capped"].values, color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--",
               linewidth=0.8, label="p = 0.05")
    ax.set_xlabel(r"$-\log_{10}$(geometric mean, HNF4A/BTNL3/BTNL8 p-values capped at $10^{-10}$)")
    ax.set_title(f"Combined HNF4A + BTNL3 + BTNL8 ranking, capped at p=1e-10\n"
                 f"(top {TOP_N} HD-abundance taxa + Candida, "
                 f"{DIAGNOSIS}, +-{NUM_WEEKS}w)")
    ax.legend(loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    cap_png = os.path.join(
        OUT_DIR, f"hnf4a_btnl_taxa_ranking_top{TOP_N}_capped1e-10_plot.png")
    save_fig(fig, cap_png)
    plt.close(fig)
    print(f"Capped plot -> {cap_png} (+ .pdf)")

    print("\nDone.")


if __name__ == "__main__":
    main()
