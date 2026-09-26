#!/usr/bin/env python3
"""
plot_figure5_enrichment_panels.py

Bar plots of the ten most significant terms for Figure 5 panels C, E and F,
read from the HD enrichment tables written by run_figure5_enrichments.sh:

  5C  GO_Biological_Process_2026   positively correlated genes
  5E  Reactome_2022                negatively correlated genes
  5F  PhenGenI_Association_2021    negatively correlated genes

The tables list every term with at least one overlapping gene, as Enrichr
does. Terms are ranked by BH-adjusted p-value (smallest first), with exact
ties broken by raw p-value and then term name; no term is selected, excluded,
renamed or reordered by hand. Bars show -log10(adjusted p-value).

Output
------
  Results/HD/Figures/Fig5C_GO_BP_positive_pm4_weeks.png      (+ .pdf)
  Results/HD/Figures/Fig5E_Reactome_negative_pm4_weeks.png   (+ .pdf)
  Results/HD/Figures/Fig5F_PhenGenI_negative_pm4_weeks.png   (+ .pdf)
"""

import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _figutil import save_fig

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DIAGNOSIS = "HD"
NUM_WEEKS = 4
TOP_N     = 10
LABEL_WRAP = 60  # characters per y-tick line (display only; names unchanged)

TABLE_DIR  = os.path.join("Results", DIAGNOSIS, "Tables")
FIGURE_DIR = os.path.join("Results", DIAGNOSIS, "Figures")

SPECIES = r"$\it{C.\,albicans}$"

# (panel, output tag used by enrichr_enrichment.py, title)
PANELS = [
    ("C", "GO_BP_positive",
     f"GO biological processes enriched among transcripts positively "
     f"associated with {SPECIES}"),
    ("E", "Reactome_negative",
     f"Reactome pathways enriched among transcripts negatively "
     f"associated with {SPECIES}"),
    ("F", "PhenGenI_negative",
     f"PhenGenI disease terms enriched among transcripts negatively "
     f"associated with {SPECIES}"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def top_terms(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Top `n` terms by adjusted p, ties broken by raw p then term name."""
    ranked = df.sort_values(["Adjusted_p_value", "P_value", "Name"],
                            kind="mergesort")
    return ranked.head(n).reset_index(drop=True)


def plot_panel(top: pd.DataFrame, title: str, n_terms: int,
               png_path: str) -> None:
    labels = [textwrap.fill(name, LABEL_WRAP) for name in top["Name"]]
    n_lines = sum(label.count("\n") + 1 for label in labels)
    fig, ax = plt.subplots(figsize=(11, max(4.5, 1.6 + 0.3 * n_lines)),
                           layout="constrained")
    y_pos = np.arange(len(top))
    ax.barh(y_pos, -np.log10(top["Adjusted_p_value"].values),
            color="#2171b5", height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--",
               linewidth=0.8, label="adjusted $P$ = 0.05")
    ax.set_xlabel(r"$-\log_{10}$(adjusted $P$ value)")
    ax.legend(loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    ax.set_axisbelow(True)
    fig.suptitle(f"{title}\n"
                 f"({DIAGNOSIS}, \u00b1{NUM_WEEKS} weeks; top {len(top)} of {n_terms} "
                 f"terms with overlapping genes, by BH-adjusted $P$ value)",
                 fontsize=11)
    save_fig(fig, png_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(FIGURE_DIR, exist_ok=True)
    suffix = f"pm{NUM_WEEKS}_weeks"

    for panel, tag, title in PANELS:
        table = os.path.join(TABLE_DIR, f"enrichr_{tag}_{suffix}.tsv")
        df = pd.read_csv(table, sep="\t")
        top = top_terms(df, TOP_N)

        print(f"\n[Fig 5{panel}] {table}: top {len(top)} of {len(df)} terms")
        disp = top[["Name", "P_value", "Adjusted_p_value", "Overlap_count",
                    "Set_size"]].copy()
        for c in ["P_value", "Adjusted_p_value"]:
            disp[c] = disp[c].map(lambda v: f"{v:.3e}")
        print(disp.to_string(index=False))

        png_path = os.path.join(FIGURE_DIR, f"Fig5{panel}_{tag}_{suffix}.png")
        plot_panel(top, title, len(df), png_path)
        print(f"  Plot -> {png_path} (+ .pdf)")

    print("\nDone.")


if __name__ == "__main__":
    main()
