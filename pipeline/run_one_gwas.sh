#!/usr/bin/env bash
# =============================================================================
# run_one_gwas.sh — End-to-end orchestrator for one GWAS
# =============================================================================
# Usage: bash pipeline/run_one_gwas.sh <gwas_key>
# Example: bash pipeline/run_one_gwas.sh af_nielsen2018
#
# Runs all 3 pipeline stages (prep + MAGMA + cell-type) for ONE GWAS using
# a per-GWAS env file. Multiple invocations can run in parallel terminals
# safely as long as gwas_key differs (each gets its own _env_<key>.env).
#
# All output logged to /tmp/run_<gwas_key>_<timestamp>.log
# =============================================================================

set -euo pipefail

GWAS_KEY="${1:?Usage: $0 <gwas_key>}"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPELINE_DIR="$PROJECT_ROOT/pipeline"
ENV_FILE="$PIPELINE_DIR/_env_${GWAS_KEY}.env"
TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="/tmp/run_${GWAS_KEY}_${TS}.log"

cd "$PROJECT_ROOT"

echo "================================================================"
echo "run_one_gwas.sh starting"
echo "  GWAS_KEY: $GWAS_KEY"
echo "  ENV_FILE: $ENV_FILE"
echo "  LOG:      $LOG_FILE"
echo "  Started:  $(date)"
echo "================================================================"

echo ""
echo "=== [1/4] Generating per-GWAS env file ==="
python3 "$PIPELINE_DIR/_emit_active_env.py" \
    --gwas-key "$GWAS_KEY" \
    --out-path "$ENV_FILE" 2>&1 | tee -a "$LOG_FILE"

echo ""
echo "=== [2/4] Stage 1: prep sumstats (host) ==="
cd "$PIPELINE_DIR"
python -u 01_prep_sumstats.py --gwas-key "$GWAS_KEY" 2>&1 | tee -a "$LOG_FILE"
cd "$PROJECT_ROOT"

echo ""
echo "=== [3/4] Stage 2: MAGMA gene-level (Docker, ~45 min) ==="
./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh "$ENV_FILE" 2>&1 | tee -a "$LOG_FILE"

echo ""
echo "=== [4/4] Stage 3: MAGMA cell-type (Docker, ~5 min) ==="
./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh "$ENV_FILE" 2>&1 | tee -a "$LOG_FILE"

OUTPUT_FOLDER=$(grep '^export OUTPUT_FOLDER=' "$ENV_FILE" | cut -d'"' -f2)

echo ""
echo "================================================================"
echo "run_one_gwas.sh DONE"
echo "  GWAS_KEY: $GWAS_KEY"
echo "  Finished: $(date)"
echo "  Log:      $LOG_FILE"
echo "  Outputs:  $PROJECT_ROOT/results/$OUTPUT_FOLDER/"
echo "================================================================"
