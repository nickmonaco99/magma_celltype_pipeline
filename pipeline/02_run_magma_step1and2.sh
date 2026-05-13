#!/usr/bin/env bash
# =============================================================================
# 02_run_magma_step1and2.sh — MAGMA annotation + gene-level analysis
# =============================================================================
# Mirrors duncan_repo/MAGMA/1.annotationAndGeneAnalysis.sh, with three changes:
#   1. ALL GWAS-specific values are read from 00_config.yaml via the
#      _active_gwas.env helper (no hardcoded paths or file names in this script).
#      To switch GWAS: edit `active_gwas:` in pipeline/00_config.yaml.
#      No edits to this script needed when adding new GWAS.
#   2. Output tee'd to ${RUN_LOG} for provenance.
#   3. magma binary path is /usr/local/bin/magma (Docker-installed).
#
# Run from project root via Docker:
#   ./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh
#
# Stages:
#   1. Annotation: maps SNPs → genes using window from config (default 35/10 kb)
#   2. Gene-level analysis: aggregates SNP p-values to gene-level using
#      1000G EUR LD reference panel
#
# Expected runtime: ~45-60 min (most time is stage 2 with ~10M SNPs)
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

# -----------------------------------------------------------------------------
# 2. Pre-flight: verify all inputs exist
# -----------------------------------------------------------------------------
echo ""
echo "=================================================================="
echo "02_run_magma_step1and2.sh — MAGMA annotation + gene analysis"
echo "=================================================================="
echo "Active GWAS:         ${ACTIVE_GWAS}  (${GWAS_NAME})"
echo "Output folder:       ${OUTPUT_FOLDER}"
echo "File prefix:         ${FILE_PREFIX}"
echo "Sumstats file:       ${SUMSTATS_FILE}"
echo "SNP loc file:        ${SNPLOC_FILE}"
echo "Gene window:         ${WIN_UP}kb upstream, ${WIN_DN}kb downstream"
echo "Reference panel:     ${G1000_PREFIX}"
echo "Required version:    ${REQ_VERSION}"
echo ""

for f in "${SUMSTATS_FILE}" "${SNPLOC_FILE}" "${GENE_LOC}" \
         "${G1000_PREFIX}.bed" "${G1000_PREFIX}.bim" "${G1000_PREFIX}.fam"; do
    if [[ ! -e "${f}" ]]; then
        echo "ERROR: required input not found: ${f}"
        exit 1
    fi
done
echo "All required inputs present."

# -----------------------------------------------------------------------------
# 3. Provenance logging
# -----------------------------------------------------------------------------
TIMESTAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)
MAGMA_VERSION=$(/usr/local/bin/magma --version 2>&1 | head -1)

{
    echo ""
    echo "================================================================="
    echo "[$TIMESTAMP] 02_run_magma_step1and2.sh"
    echo "================================================================="
    echo "MAGMA: $MAGMA_VERSION"
    echo "GWAS: ${ACTIVE_GWAS} (${GWAS_NAME})"
    echo "Folder: ${OUTPUT_FOLDER}"
    echo "File prefix: ${FILE_PREFIX}"
    echo "Sumstats: ${SUMSTATS_FILE}"
    echo "Window: ${WIN_UP},${WIN_DN}"
} | tee -a "${RUN_LOG}"

# Sanity check version matches our pin
if ! echo "$MAGMA_VERSION" | grep -q "${REQ_VERSION}"; then
    echo "ERROR: MAGMA version mismatch! Required ${REQ_VERSION}, got: $MAGMA_VERSION"
    exit 1
fi

# -----------------------------------------------------------------------------
# 4. Stage 1: annotation (SNPs → genes)
# -----------------------------------------------------------------------------
echo ""
echo "=== STAGE 1: annotation (window=${WIN_UP},${WIN_DN}) ==="
echo "Expected runtime: ~1 minute"
date -u +'Started: %Y-%m-%dT%H:%M:%SZ'

/usr/local/bin/magma \
    --annotate window=${WIN_UP},${WIN_DN} \
    --snp-loc "${SNPLOC_FILE}" \
    --gene-loc "${GENE_LOC}" \
    --out "${STEP1_PREFIX}" \
    2>&1 | tee -a "${RUN_LOG}"

date -u +'Finished: %Y-%m-%dT%H:%M:%SZ'

# -----------------------------------------------------------------------------
# 5. Stage 2: gene-level analysis
# -----------------------------------------------------------------------------
echo ""
echo "=== STAGE 2: gene-level analysis ==="
echo "Expected runtime: ~30-60 minutes"
date -u +'Started: %Y-%m-%dT%H:%M:%SZ'

/usr/local/bin/magma \
    --bfile "${G1000_PREFIX}" \
    --pval "${SUMSTATS_FILE}" use=${SNP_COL},${P_COL} ncol=${N_COL} \
    --gene-annot "${STEP1_ANNOT}" \
    --out "${STEP2_PREFIX}" \
    2>&1 | tee -a "${RUN_LOG}"

date -u +'Finished: %Y-%m-%dT%H:%M:%SZ'

# -----------------------------------------------------------------------------
# 6. Summary
# -----------------------------------------------------------------------------
echo ""
echo "=================================================================="
echo "✓ 02_run_magma_step1and2.sh — DONE"
echo "=================================================================="
echo "Outputs:"
ls -la "${STEP1_ANNOT}" \
       "${STEP2_GENES_RAW}" \
       "${STEP2_GENES_OUT}" \
       "${STEP2_LOG}" 2>&1
echo ""
echo "Next: open pipeline/02_qc_gene_results.py to validate gene-level output"
echo "      before running cell-type analysis (03_run_magma_celltype.sh)"
