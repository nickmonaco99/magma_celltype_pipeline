#!/usr/bin/env bash
# =============================================================================
# 02_run_magma_step1and2.sh — MAGMA annotation + gene-level analysis
# =============================================================================
# Mirrors duncan_repo/MAGMA/1.annotationAndGeneAnalysis.sh verbatim, with three
# changes:
#   1. Parameterized folder/file (was hardcoded "SCZ", now reads from config)
#   2. Output tee'd to results/<folder>/run.log for provenance
#   3. magma binary path is /usr/local/bin/magma (Docker-installed)
#
# Run from project root via Docker:
#   ./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh
#
# Stages:
#   1. Annotation: maps SNPs → genes using 35kb up / 10kb down window
#   2. Gene-level analysis: aggregates SNP p-values to gene-level using
#      1000G EUR LD reference panel
#
# Expected runtime: ~45-60 min (most time is stage 2 with ~10M SNPs)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# 1. Read config (yq parses YAML; installed in Docker image? if not, use python)
# -----------------------------------------------------------------------------
# We use a Python one-liner since Python is in the Docker image but yq may not be.
# Reads the active GWAS metadata into shell variables.

read_config() {
    python3 -c "
import yaml, sys
with open('pipeline/00_config.yaml') as f: c = yaml.safe_load(f)
g = c['gwas'][c['active_gwas']]
print(f\"FOLDER={g['output_folder']}\")
print(f\"WIN_UP={c['magma']['gene_window_upstream_kb']}\")
print(f\"WIN_DN={c['magma']['gene_window_downstream_kb']}\")
print(f\"REQ_VERSION={c['magma']['required_version']}\")
"
}

eval "$(read_config)"

# Hardcoded per our convention (Duncan also hardcoded these in his scripts).
# These match what 01_prep_sumstats.py writes.
FILE="insomnia_ukb.no_heading"
SNP_COL="ID"
P_COL="PVAL"
N_COL="NEFF"

# Paths (relative to project root = /work inside Docker)
SUMSTATS_FILE="results/${FOLDER}/${FILE}"
SNPLOC_FILE="results/${FOLDER}/snploc_insomnia_ukb"
RUN_LOG="results/${FOLDER}/run.log"
GENE_LOC="auxfiles/NCBI37.3.gene.loc"
G1000_PREFIX="auxfiles/g1000_eur"

# -----------------------------------------------------------------------------
# 2. Pre-flight: verify all inputs exist
# -----------------------------------------------------------------------------
echo ""
echo "=================================================================="
echo "02_run_magma_step1and2.sh — MAGMA annotation + gene analysis"
echo "=================================================================="
echo "Active GWAS folder:  ${FOLDER}"
echo "Sumstats file:       ${SUMSTATS_FILE}"
echo "SNP loc file:        ${SNPLOC_FILE}"
echo "Gene window:         ${WIN_UP}kb upstream, ${WIN_DN}kb downstream"
echo "Reference panel:     ${G1000_PREFIX}"
echo "Required version:    ${REQ_VERSION}"
echo ""

for f in "${SUMSTATS_FILE}" "${SNPLOC_FILE}" "${GENE_LOC}" "${G1000_PREFIX}.bed" "${G1000_PREFIX}.bim" "${G1000_PREFIX}.fam"; do
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
    echo "Folder: ${FOLDER}"
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
    --out "results/${FOLDER}/${FILE}.step1" \
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
    --gene-annot "results/${FOLDER}/${FILE}.step1.genes.annot" \
    --out "results/${FOLDER}/${FILE}.step2" \
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
ls -la "results/${FOLDER}/${FILE}.step1.genes.annot" \
       "results/${FOLDER}/${FILE}.step2.genes.raw" \
       "results/${FOLDER}/${FILE}.step2.genes.out" \
       "results/${FOLDER}/${FILE}.step2.log" 2>&1
echo ""
echo "Next: open pipeline/02_qc_gene_results.py in Cursor/your editor"
echo "      to validate the gene-level output before cell-type analysis"
