#!/usr/bin/env bash
# runme.sh — Manuscript_Candida_Albicans pipeline (Kraken Plus PFP 2024-10 reference only)
#
# Orchestrates the full analysis end-to-end:
#
#   Step 1  Filter host transcriptomics to known genes
#   Step 2  Per-diagnosis Candida vs. host-gene correlations
#   Step 3  ENCODE/ChEA TF enrichment on positively correlated genes
#   Step 4  GO (BP / MF / CC) enrichment on negatively correlated genes
#   Step 4b Figure 5C/E/F (HD): GO BP on positively correlated genes,
#           Reactome 2022 and PhenGenI 2021 on negatively correlated
#           genes, plus panel plots  (run_figure5_enrichments.sh)
#   Step 5  BTNL3 / BTNL8 direct analysis (scatter + week sensitivity)
#   Step 6  HNF4A ENCODE enrichment sensitivity to the +-week window
#   Step 7  HNF4A + BTNL3/BTNL8 combined ranking of 13 curated prevalent
#           human-associated fungi  (uncapped + capped)
#   Step 8  HNF4A + BTNL3/BTNL8 combined ranking of the top-100 most
#           abundant HD taxa plus Candida albicans  (uncapped + capped)
#
# Every figure is written as both .png and .pdf (via _figutil.save_fig).
#
# Usage:
#   bash runme.sh
#
# Prerequisites:
#   Data/ must contain:
#     species_table_clr_202410.tsv
#     PRJNA398089_MGX_merged.tsv
#     hmp2_metadata_2018-08-20.csv
#     host_tx_counts.tsv
#     ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt
#     GO_Biological_Process_2026.gmt
#     GO_Molecular_Function_2026.gmt
#     GO_Cellular_Component_2026.gmt
#     Reactome_2022.gmt
#     PhenGenI_Association_2021.gmt
#     enrichr_library_manifest.tsv
#   Python with: pandas, numpy, scipy, matplotlib, tqdm

set -euo pipefail

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
NUM_WEEKS=4                              # +- temporal matching window (weeks)
TARGET="S__Candida albicans"             # target species
DIAGNOSES="HD CD UC"                     # diagnosis groups for steps 2-4
KEEP_TYPES=""                            # gene-type filter (empty = all)
TX_COUNTS="Data/host_tx_counts.tsv"
TX_FILTERED="Data/host_tx_counts_filtered.tsv"
TX_ANNOTATIONS="Data/gene_annotations.tsv"
BACKGROUND="${TX_FILTERED}"
ENRICHR_P_THRESHOLD=0.05

echo "=== Manuscript_Candida_Albicans pipeline (Kraken Plus PFP 2024-10) ==="
echo "Date:        $(date)"
echo "Window:      +-${NUM_WEEKS} weeks"
echo "Target:      ${TARGET}"
echo "Diagnoses:   ${DIAGNOSES}"
echo ""

# ---------------------------------------------------------------------------
# Step 1: Pre-filter transcriptomics to known genes
# ---------------------------------------------------------------------------
if [ ! -f "${TX_FILTERED}" ]; then
    if [ -f "${TX_ANNOTATIONS}" ]; then
        echo "[Step 1] Filtering transcriptomics to known genes (offline, cached) ..."
    else
        echo "[Step 1] Filtering transcriptomics to known genes (querying mygene.info) ..."
    fi
    FILTER_CMD=(python filter_known_genes.py
        --input "${TX_COUNTS}"
        --out-counts "${TX_FILTERED}"
        --out-annotations "${TX_ANNOTATIONS}")
    if [ -n "${KEEP_TYPES}" ]; then
        FILTER_CMD+=(--keep-types ${KEEP_TYPES})
    fi
    "${FILTER_CMD[@]}"
    echo ""
else
    echo "[Step 1] Filtered transcriptomics already exists, skipping."
    echo "         (delete ${TX_FILTERED} to re-run; pipeline uses"
    echo "          ${TX_ANNOTATIONS} as an offline mygene.info cache)"
    echo ""
fi

# ---------------------------------------------------------------------------
# Step 2: Correlate Candida vs host gene expression, per diagnosis
# ---------------------------------------------------------------------------
echo "[Step 2] Correlating Candida CLR with gene expression ..."
python correlate_candida_genes.py \
    --num-weeks "${NUM_WEEKS}" \
    --target "${TARGET}" \
    --diagnoses ${DIAGNOSES} \
    --tx-counts "${TX_FILTERED}"
echo ""

# ---------------------------------------------------------------------------
# Step 3: ENCODE/ChEA TF enrichment (positively correlated genes)
# ---------------------------------------------------------------------------
echo "[Step 3] TF enrichment (positive correlation, ENCODE/ChEA) ..."
for diag in ${DIAGNOSES}; do
    echo "  --- ${diag} ---"
    python enrichr_enrichment.py \
        --library "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" \
        --output-tag "ENCODE_ChEA_TFs" \
        --diagnoses "${diag}" \
        --p-threshold "${ENRICHR_P_THRESHOLD}" \
        --corr-sign "pos" \
        --num-weeks "${NUM_WEEKS}" \
        --background "${BACKGROUND}"
done
echo ""

# ---------------------------------------------------------------------------
# Step 4: GO enrichment (negatively correlated genes)
# ---------------------------------------------------------------------------
GO_LIBRARIES=(
    "GO_Biological_Process_2026"
    "GO_Molecular_Function_2026"
    "GO_Cellular_Component_2026"
)
GO_TAGS=(
    "GO_BP"
    "GO_MF"
    "GO_CC"
)

echo "[Step 4] GO enrichment (negative correlation) ..."
for i in "${!GO_LIBRARIES[@]}"; do
    LIB="${GO_LIBRARIES[$i]}"
    TAG="${GO_TAGS[$i]}"
    echo "  Library: ${TAG}"
    for diag in ${DIAGNOSES}; do
        echo "    --- ${diag} ---"
        python enrichr_enrichment.py \
            --library "${LIB}" \
            --output-tag "${TAG}" \
            --diagnoses "${diag}" \
            --p-threshold "${ENRICHR_P_THRESHOLD}" \
            --corr-sign "neg" \
            --num-weeks "${NUM_WEEKS}" \
            --background "${BACKGROUND}"
    done
done
echo ""

# ---------------------------------------------------------------------------
# Step 4b: Figure 5C/E/F enrichments (HD) + panel plots
# ---------------------------------------------------------------------------
echo "[Step 4b] Figure 5 enrichments (HD: GO BP pos, Reactome neg, PhenGenI neg) ..."
bash run_figure5_enrichments.sh
echo ""

# ---------------------------------------------------------------------------
# Step 5: BTNL3 / BTNL8 direct analysis
# ---------------------------------------------------------------------------
echo "[Step 5] BTNL3 / BTNL8 direct correlation and week sensitivity ..."
python btnl_analysis.py
echo ""

# ---------------------------------------------------------------------------
# Step 6: HNF4A enrichment sensitivity to the +-week window
# ---------------------------------------------------------------------------
echo "[Step 6] HNF4A enrichment sensitivity to +-week window (HD/CD/UC) ..."
python hnf4a_week_sensitivity.py
echo ""

# ---------------------------------------------------------------------------
# Step 7: HNF4A + BTNL combined ranking of the 13 curated fungi
# ---------------------------------------------------------------------------
echo "[Step 7] HNF4A + BTNL3/BTNL8 ranking of 13 curated fungi ..."
python hnf4a_btnl_fungal_ranking.py
echo ""

# ---------------------------------------------------------------------------
# Step 8: HNF4A + BTNL combined ranking of top-100 HD taxa + Candida
# ---------------------------------------------------------------------------
echo "[Step 8] HNF4A + BTNL3/BTNL8 ranking of top-100 HD taxa + Candida ..."
python hnf4a_btnl_taxa_ranking.py
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "=== Pipeline complete ==="
echo ""
echo "Outputs:"
echo "  Results/<HD|CD|UC>/Tables/   correlations + Enrichr results"
echo "  Results/<HD|CD|UC>/Figures/  scatter plots for genes of interest"
echo "  Results/BTNL/                BTNL3/BTNL8 scatter + sensitivity"
echo "  Results-HNF4A-sensitivity/   HNF4A enrichment vs. +-week window"
echo "  Results-HNF4A-fungi/         13 curated fungi, HNF4A+BTNL3+BTNL8"
echo "                               (uncapped + capped at p=1e-10)"
echo "  Results-HNF4A-taxa/          top-100 HD taxa + Candida, same"
echo "                               (uncapped + capped at p=1e-10)"
echo ""
echo "Every .png plot has a matching .pdf alongside it."
echo "See methods.txt in the root folder and in each results folder for details."
