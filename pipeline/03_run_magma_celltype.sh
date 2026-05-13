#!/usr/bin/env bash
# =============================================================================
# 03_run_magma_celltype.sh — MAGMA gene-property (cell-type) analysis
# =============================================================================
# Mirrors duncan_repo/MAGMA/2.genePropertyAnalysis.sh with three changes:
#   1. Parameterized folder/file (was hardcoded "SCZ")
#   2. Output tee'd to results/<folder>/run.log for provenance
#   3. magma binary path is /usr/local/bin/magma (Docker-installed)
#
# Run from project root via Docker:
#   ./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh
#
# Requires:
#   - results/<FOLDER>/<FILE>.step2.genes.raw (produced by 02_run_magma_step1and2.sh)
#   - gene-level/Siletti_l2_conti_specificity_matrix.txt
#
# Output:
#   - results/<FOLDER>/Siletti_l2_conti-spe_<FOLDER>.gsa.out
#     (461 rows, one per Siletti cluster: VARIABLE TYPE NGENES BETA BETA_STD SE P)
#
# Expected runtime: ~5 minutes
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# 1. Read config
# -----------------------------------------------------------------------------
read_config() {
    python3 -c "
import yaml, sys
with open('pipeline/00_config.yaml') as f: c = yaml.safe_load(f)
g = c['gwas'][c['active_gwas']]
print(f\"FOLDER={g['output_folder']}\")
print(f\"DIRECTION={c['magma']['model_direction']}\")
print(f\"SPEC_MATRIX={c['paths']['specificity_matrix']}\")
print(f\"REQ_VERSION={c['magma']['required_version']}\")
"
}

eval "$(read_config)"

FILE="insomnia_ukb.no_heading"
OUTFILE_NAME="Siletti_l2_conti-spe_${FOLDER}"

GENES_RAW="results/${FOLDER}/${FILE}.step2.genes.raw"
RUN_LOG="results/${FOLDER}/run.log"
OUT_PREFIX="results/${FOLDER}/${OUTFILE_NAME}"

# -----------------------------------------------------------------------------
# 2. Pre-flight checks
# -----------------------------------------------------------------------------
echo ""
echo "=================================================================="
echo "03_run_magma_celltype.sh — MAGMA cell-type analysis"
echo "=================================================================="
echo "Active GWAS folder:    ${FOLDER}"
echo "Gene results file:     ${GENES_RAW}"
echo "Specificity matrix:    ${SPEC_MATRIX}"
echo "Test direction:        ${DIRECTION} (one-sided positive)"
echo "Output prefix:         ${OUT_PREFIX}"
echo ""

for f in "${GENES_RAW}" "${SPEC_MATRIX}"; do
    if [[ ! -e "${f}" ]]; then
        echo "ERROR: required input not found: ${f}"
        exit 1
    fi
done

# Quick check on specificity matrix shape (should have 17,097 genes + header = 17,098 lines
# and 462 columns: GENE + 461 cluster columns)
SPEC_LINES=$(wc -l < "${SPEC_MATRIX}")
SPEC_COLS=$(head -1 "${SPEC_MATRIX}" | tr '\t' '\n' | wc -l)
echo "Specificity matrix:    ${SPEC_LINES} lines × ${SPEC_COLS} columns"
echo "                       (expected: 17098 × 462)"

# -----------------------------------------------------------------------------
# 3. Provenance logging
# -----------------------------------------------------------------------------
TIMESTAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)
MAGMA_VERSION=$(/usr/local/bin/magma --version 2>&1 | head -1)

{
    echo ""
    echo "================================================================="
    echo "[$TIMESTAMP] 03_run_magma_celltype.sh"
    echo "================================================================="
    echo "MAGMA: $MAGMA_VERSION"
    echo "Folder: ${FOLDER}"
    echo "Specificity matrix: ${SPEC_MATRIX} (${SPEC_LINES} lines × ${SPEC_COLS} cols)"
    echo "Test direction: ${DIRECTION}"
} | tee -a "${RUN_LOG}"

if ! echo "$MAGMA_VERSION" | grep -q "${REQ_VERSION}"; then
    echo "ERROR: MAGMA version mismatch! Required ${REQ_VERSION}, got: $MAGMA_VERSION"
    exit 1
fi

# -----------------------------------------------------------------------------
# 4. Run cell-type analysis
# -----------------------------------------------------------------------------
echo ""
echo "=== Cell-type analysis (461 Siletti clusters) ==="
echo "Expected runtime: ~5 minutes"
date -u +'Started: %Y-%m-%dT%H:%M:%SZ'

/usr/local/bin/magma \
    --gene-results "${GENES_RAW}" \
    --gene-covar "${SPEC_MATRIX}" \
    --model direction=${DIRECTION} \
    --out "${OUT_PREFIX}" \
    2>&1 | tee -a "${RUN_LOG}"

date -u +'Finished: %Y-%m-%dT%H:%M:%SZ'

# -----------------------------------------------------------------------------
# 5. Quick result summary
# -----------------------------------------------------------------------------
GSA_OUT="${OUT_PREFIX}.gsa.out"
if [[ ! -e "${GSA_OUT}" ]]; then
    echo "ERROR: expected output ${GSA_OUT} not found"
    exit 1
fi

echo ""
echo "=================================================================="
echo "✓ 03_run_magma_celltype.sh — DONE"
echo "=================================================================="
echo ""
echo "Output: ${GSA_OUT}"
echo ""
echo "=== File header (first 8 lines) ==="
head -8 "${GSA_OUT}"
echo ""
echo "=== Top 10 cell types by p-value ==="
# Skip 4 header lines + column header, then sort by P column (column 7) ascending
tail -n +6 "${GSA_OUT}" | sort -k7 -g | head -10

echo ""
echo "Next: open pipeline/04_celltype_results.py to plot results,"
echo "      identify Bonferroni-significant cell types (P < 0.05/461 = 1.08e-4),"
echo "      and map Cluster# → splatter / Siletti supercluster"
