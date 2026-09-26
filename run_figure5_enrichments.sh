#!/usr/bin/env bash
# run_figure5_enrichments.sh — Figure 5C / 5E / 5F enrichment panels (HD only)
#
# Runs the three additional HD enrichments behind Figure 5 panels C, E and F
# with enrichr_enrichment.py (one-sided Fisher's exact test computed as
# Enrichr does, BH correction over the terms with >=1 overlapping gene,
# filtered-transcript background), on the existing +-4-week HD correlation
# table, then plots the panels:
#
#   5C  GO_Biological_Process_2026   positively correlated genes
#   5E  Reactome_2022                negatively correlated genes
#   5F  PhenGenI_Association_2021    negatively correlated genes
#
# Libraries are read from the bundled Data/*.gmt files and checked against
# the SHA-256 recorded in Data/enrichr_library_manifest.tsv, so this step
# runs offline and cannot silently pick up a different Enrichr release.
#
# Usage:
#   bash run_figure5_enrichments.sh      (also called from runme.sh, Step 4b)
#
# Prerequisites:
#   Results/HD/Tables/candida_gene_correlations_pm4_weeks_with_p-values.tsv
#   (runme.sh Step 2), Data/host_tx_counts_filtered.tsv, the three GMT files
#   and Data/enrichr_library_manifest.tsv

set -euo pipefail

# ---------------------------------------------------------------------------
# Parameters (identical to runme.sh)
# ---------------------------------------------------------------------------
NUM_WEEKS=4
DIAGNOSIS="HD"
BACKGROUND="Data/host_tx_counts_filtered.tsv"
ENRICHR_P_THRESHOLD=0.05
MANIFEST="Data/enrichr_library_manifest.tsv"

LIBRARIES=(
    "GO_Biological_Process_2026"
    "Reactome_2022"
    "PhenGenI_Association_2021"
)
TAGS=(
    "GO_BP_positive"
    "Reactome_negative"
    "PhenGenI_negative"
)
SIGNS=(
    "pos"
    "neg"
    "neg"
)

# ---------------------------------------------------------------------------
# Verify the bundled libraries (enrichr_enrichment.py would otherwise
# download a missing GMT from the live Enrichr endpoint)
# ---------------------------------------------------------------------------
for LIB in "${LIBRARIES[@]}"; do
    GMT="Data/${LIB}.gmt"
    if [ ! -f "${GMT}" ]; then
        echo "ERROR: ${GMT} not found; refusing to download a replacement." >&2
        exit 1
    fi
    EXPECTED=$(awk -F'\t' -v lib="${LIB}" '$1 == lib {print $6}' "${MANIFEST}")
    ACTUAL=$(python -c 'import hashlib, sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "${GMT}")
    if [ -z "${EXPECTED}" ] || [ "${EXPECTED}" != "${ACTUAL}" ]; then
        echo "ERROR: ${GMT} SHA-256 ${ACTUAL} does not match ${MANIFEST} (${EXPECTED:-missing})." >&2
        exit 1
    fi
    echo "  ${GMT}: SHA-256 OK"
done

# ---------------------------------------------------------------------------
# Enrichment
# ---------------------------------------------------------------------------
for i in "${!LIBRARIES[@]}"; do
    LIB="${LIBRARIES[$i]}"
    TAG="${TAGS[$i]}"
    SIGN="${SIGNS[$i]}"
    echo "  Library: ${LIB} (corr ${SIGN}, ${DIAGNOSIS})"
    python enrichr_enrichment.py \
        --library "${LIB}" \
        --output-tag "${TAG}" \
        --diagnoses "${DIAGNOSIS}" \
        --p-threshold "${ENRICHR_P_THRESHOLD}" \
        --corr-sign "${SIGN}" \
        --num-weeks "${NUM_WEEKS}" \
        --background "${BACKGROUND}"
done

# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
python plot_figure5_enrichment_panels.py
