#!/usr/bin/env bash
# =============================================================================
# 03_run_magma_celltype.sh — MAGMA gene-property (cell-type) analysis
# =============================================================================
# Mirrors duncan_repo/MAGMA/2.genePropertyAnalysis.sh, with three changes:
#   1. ALL GWAS-specific values are read from 00_config.yaml via the
#      _active_gwas.env helper (no hardcoded paths or file names in this script).
#      To switch GWAS: edit `active_gwas:` in pipeline/00_config.yaml.
#   2. Output tee'd to ${RUN_LOG} for provenance.
#   3. magma binary path is /usr/local/bin/magma (Docker-installed).
#
# Run from project root via Docker:
#   ./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh
#
# Requires:
#   - ${STEP2_GENES_RAW} (produced by 02_run_magma_step1and2.sh)
#   - ${SPECIFICITY_MATRIX}
#
# Output:
#   - ${CELLTYPE_OUT_PREFIX}.gsa.out
#     (461 rows, one per Siletti cluster: VARIABLE TYPE NGENES BETA BETA_STD SE P)
#
# Expected runtime: ~5 minutes
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# 1. Load active-GWAS variables from config
# -----------------------------------------------------------------------------
ENV_FILE="${1:-pipeline/_active_gwas.env}"
if [[ "$ENV_FILE" == "pipeline/_active_gwas.env" ]]; then
    python3 pipeline/_emit_active_env.py
fi
source "$ENV_FILE"

GSA_OUT="${CELLTYPE_OUT_PREFIX}.gsa.out"

# -----------------------------------------------------------------------------
# 2. Pre-flight checks
# -----------------------------------------------------------------------------
echo ""
echo "=================================================================="
echo "03_run_magma_celltype.sh — MAGMA cell-type analysis"
echo "=================================================================="
echo "Active GWAS:           ${ACTIVE_GWAS}  (${GWAS_NAME})"
echo "Output folder:         ${OUTPUT_FOLDER}"
echo "Gene results file:     ${STEP2_GENES_RAW}"
echo "Specificity matrix:    ${SPECIFICITY_MATRIX}"
echo "Test direction:        ${MODEL_DIRECTION} (one-sided positive)"
echo "Output prefix:         ${CELLTYPE_OUT_PREFIX}"
echo ""

for f in "${STEP2_GENES_RAW}" "${SPECIFICITY_MATRIX}"; do
    if [[ ! -e "${f}" ]]; then
        echo "ERROR: required input not found: ${f}"
        exit 1
    fi
done

# Quick check on specificity matrix shape (should have 17,097 genes + header = 17,098 lines
# and 462 columns: GENE + 461 cluster columns)
SPEC_LINES=$(wc -l < "${SPECIFICITY_MATRIX}")
SPEC_COLS=$(head -1 "${SPECIFICITY_MATRIX}" | tr '\t' '\n' | wc -l)
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
    echo "GWAS: ${ACTIVE_GWAS} (${GWAS_NAME})"
    echo "Folder: ${OUTPUT_FOLDER}"
    echo "Specificity matrix: ${SPECIFICITY_MATRIX} (${SPEC_LINES} lines × ${SPEC_COLS} cols)"
    echo "Test direction: ${MODEL_DIRECTION}"
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
    --gene-results "${STEP2_GENES_RAW}" \
    --gene-covar "${SPECIFICITY_MATRIX}" \
    --model direction=${MODEL_DIRECTION} \
    --out "${CELLTYPE_OUT_PREFIX}" \
    2>&1 | tee -a "${RUN_LOG}"

date -u +'Finished: %Y-%m-%dT%H:%M:%SZ'

# -----------------------------------------------------------------------------
# 5. Quick result summary
# -----------------------------------------------------------------------------
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
echo "Next steps:"
echo "  python pipeline/snapshot_results.py ${ACTIVE_GWAS}"
echo "  python pipeline/04_annotate_results.py"
echo "  python pipeline/05_make_plots.py"
