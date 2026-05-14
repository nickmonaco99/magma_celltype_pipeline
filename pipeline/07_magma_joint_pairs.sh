#!/usr/bin/env bash
# MAGMA --model joint-pairs cell-type conditional analysis.
#
# Run inside Docker via:
#   ./run_in_docker.sh bash pipeline/07_magma_joint_pairs.sh \
#       <genes_raw> <spec_covar> <out_prefix>
#
# Mirrors Duncan's duncan_repo/MAGMA/4.conditionalAnalysis.sh flags exactly:
#   --gene-results <step2.genes.raw>
#   --gene-covar   <sig-only specificity matrix>
#   --model        joint-pairs direction=greater
#   --out          <output prefix>

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <genes_raw> <spec_covar> <out_prefix>"
    exit 1
fi

GENES_RAW="$1"
SPEC_COVAR="$2"
OUT_PREFIX="$3"

echo "MAGMA joint-pairs conditional analysis"
echo "  --gene-results: ${GENES_RAW}"
echo "  --gene-covar:   ${SPEC_COVAR}"
echo "  --out:          ${OUT_PREFIX}"
echo ""

for f in "${GENES_RAW}" "${SPEC_COVAR}"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: missing input: $f"
        exit 1
    fi
done

# Matches Duncan's 4.conditionalAnalysis.sh exactly (same MAGMA flags)
/usr/local/bin/magma \
    --gene-results "${GENES_RAW}" \
    --gene-covar   "${SPEC_COVAR}" \
    --model        joint-pairs direction=greater \
    --out          "${OUT_PREFIX}"

echo ""
echo "Joint-pairs output written to: ${OUT_PREFIX}.gsa.out"
