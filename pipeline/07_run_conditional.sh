#!/usr/bin/env bash
# ============================================================================
# Forward-selection conditional analysis using Duncan 2025's exact algorithm.
#
# Workflow:
#   1. Filter Siletti specificity matrix to Bonferroni-significant clusters
#   2. Run MAGMA --model joint-pairs direction=greater in Docker
#   3. Forward selection per Duncan's R (with a/b swap adaptation)
#
# Usage:
#   ./pipeline/07_run_conditional.sh scz_pgc3_2022
#   ./pipeline/07_run_conditional.sh chronotype_loh2018_ukb
#   ./pipeline/07_run_conditional.sh bmi_yengo2018
# ============================================================================

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <gwas_key>"
    echo ""
    echo "Available gwas_keys:"
    echo "  scz_pgc3_2022"
    echo "  chronotype_loh2018_ukb"
    echo "  bmi_yengo2018"
    exit 1
fi

GWAS_KEY="$1"

case "$GWAS_KEY" in
    scz_pgc3_2022)
        RESULTS_DIR="results/SCZ_PGC3"
        MARGINAL_GSA="${RESULTS_DIR}/Siletti_l2_conti-spe_SCZ_PGC3.gsa.out"
        GENES_RAW="${RESULTS_DIR}/scz_pgc3_2022.no_heading.step2.genes.raw"
        ;;
    chronotype_loh2018_ukb)
        RESULTS_DIR="results/Chronotype_Loh2018"
        MARGINAL_GSA="${RESULTS_DIR}/Siletti_l2_conti-spe_Chronotype_Loh2018.gsa.out"
        GENES_RAW="${RESULTS_DIR}/chronotype_loh2018_ukb.no_heading.step2.genes.raw"
        ;;
    bmi_yengo2018)
        RESULTS_DIR="results/BMI_Yengo2018"
        MARGINAL_GSA="${RESULTS_DIR}/Siletti_l2_conti-spe_BMI_Yengo2018.gsa.out"
        GENES_RAW="${RESULTS_DIR}/bmi_yengo2018.no_heading.step2.genes.raw"
        ;;
    *)
        echo "ERROR: Unknown gwas_key: $GWAS_KEY"
        exit 1
        ;;
esac

COND_DIR="${RESULTS_DIR}/conditional"
mkdir -p "$COND_DIR"
SPEC_FILTERED="${COND_DIR}/spec_matrix_sig_only.txt"
JOINT_PREFIX="${COND_DIR}/joint_pairs"
JOINT_OUT="${JOINT_PREFIX}.gsa.out"
INDEP_OUT="${COND_DIR}/independent_clusters.txt"

for f in "$MARGINAL_GSA" "$GENES_RAW" "gene-level/Siletti_l2_conti_specificity_matrix.txt"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: missing required input: $f"
        exit 1
    fi
done

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found in PATH."
    exit 1
fi

echo "================================================================"
echo "Conditional analysis for ${GWAS_KEY}"
echo "================================================================"
echo "Marginal:   ${MARGINAL_GSA}"
echo "Genes raw:  ${GENES_RAW}"
echo "Output dir: ${COND_DIR}"
echo ""

# Step 1: Filter spec matrix
echo "[1/3] Filtering specificity matrix to Bonferroni-significant clusters..."
python pipeline/07_filter_spec_matrix.py \
    --marginal-gsa "$MARGINAL_GSA" \
    --spec-matrix  "gene-level/Siletti_l2_conti_specificity_matrix.txt" \
    --output       "$SPEC_FILTERED"
echo ""

# Step 2: MAGMA joint-pairs in Docker
echo "[2/3] Running MAGMA --model joint-pairs direction=greater (in Docker)..."
./run_in_docker.sh bash pipeline/07_magma_joint_pairs.sh \
    "${GENES_RAW}" \
    "${SPEC_FILTERED}" \
    "${JOINT_PREFIX}"
echo ""

# Step 3: Duncan's R forward selection (with a/b swap fix)
echo "[3/3] Running Duncan's R forward selection..."
Rscript pipeline/07_forward_selection.R \
    "$MARGINAL_GSA" \
    "$JOINT_OUT" \
    "$INDEP_OUT"
echo ""

echo "================================================================"
echo "DONE. Independent clusters for ${GWAS_KEY}:"
echo "================================================================"
cat "$INDEP_OUT"
