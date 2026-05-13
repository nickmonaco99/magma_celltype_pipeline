#!/usr/bin/env bash
# =============================================================================
# run_in_docker.sh — execute a command inside the MAGMA container
# =============================================================================
# Mounts the project root at /work, sets WORKDIR there, runs your command.
#
# Examples:
#   ./run_in_docker.sh /usr/local/bin/magma --version
#   ./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh
#   ./run_in_docker.sh bash             # interactive shell
#
# Assumes you've built the image:
#   docker build -t magma_celltype:latest .
# =============================================================================

set -euo pipefail

# Resolve project root (the directory containing this script)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Pick interactive vs non-interactive based on whether stdin is a TTY
TTY_FLAGS=""
if [[ -t 0 ]]; then
    TTY_FLAGS="-it"
fi

# Verify image exists; helpful error if not
if ! docker image inspect magma_celltype:latest >/dev/null 2>&1; then
    echo "ERROR: Docker image 'magma_celltype:latest' not found."
    echo "Build it first: cd $PROJECT_ROOT && docker build -t magma_celltype:latest ."
    exit 1
fi

# If no args, drop into interactive bash
if [[ $# -eq 0 ]]; then
    set -- "/bin/bash"
fi

docker run --rm $TTY_FLAGS \
    -v "$PROJECT_ROOT:/work" \
    -w /work \
    magma_celltype:latest \
    "$@"
