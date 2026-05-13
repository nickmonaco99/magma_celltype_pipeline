#!/usr/bin/env bash
# =============================================================================
# setup.sh — One-time bootstrap for the MAGMA cell-type pipeline (macOS + Docker)
# =============================================================================
# Usage:
#   cd /Users/nickmonaco/Desktop/"Science Side Projects"/"Allen Brain Seattle Trip"
#   bash magma_celltype_pipeline/setup.sh
# =============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT"

# -----------------------------------------------------------------------------
# 1. Create directory structure
# -----------------------------------------------------------------------------
echo ""
echo "=== Creating directory structure ==="
mkdir -p auxfiles gene-level gwas_sumstats results figures scratch magma_binary
mkdir -p gwas_sumstats/insomnia_watanabe2022_ukb
ls -d */ 2>/dev/null || true

# -----------------------------------------------------------------------------
# 2. Locate Duncan's repo
# -----------------------------------------------------------------------------
echo ""
echo "=== Locating Duncan's repo ==="
DUNCAN_REPO_SOURCE="$PROJECT_ROOT/../linking_cell_types_to_brain_phenotypes-main"
if [[ -d "$DUNCAN_REPO_SOURCE" ]]; then
    echo "Found at: $DUNCAN_REPO_SOURCE"
    if [[ ! -d "$PROJECT_ROOT/duncan_repo" ]]; then
        echo "Moving to: $PROJECT_ROOT/duncan_repo"
        mv "$DUNCAN_REPO_SOURCE" "$PROJECT_ROOT/duncan_repo"
    else
        echo "Already moved (duncan_repo/ exists)"
    fi
elif [[ -d "$PROJECT_ROOT/duncan_repo" ]]; then
    echo "Already in place (duncan_repo/ exists)"
else
    echo "WARNING: Duncan's repo not found at either location."
fi

# Symlink the gene.loc file so paths match Duncan's expected layout
if [[ -f "duncan_repo/Data/NCBI37.3.gene.loc" && ! -e "auxfiles/NCBI37.3.gene.loc" ]]; then
    echo "Symlinking gene.loc into auxfiles/"
    ln -s "../duncan_repo/Data/NCBI37.3.gene.loc" "auxfiles/NCBI37.3.gene.loc"
fi

# -----------------------------------------------------------------------------
# 3. Check Docker is available
# -----------------------------------------------------------------------------
echo ""
echo "=== Checking Docker ==="
ARCH=$(uname -m)
echo "Mac architecture: $ARCH"
if [[ "$ARCH" == "arm64" ]]; then
    echo "  → Apple Silicon. Using Docker (Linux containers) to run MAGMA."
elif [[ "$ARCH" == "x86_64" ]]; then
    echo "  → Intel Mac. Docker still recommended for parity with Duncan's Linux env."
fi

if ! command -v docker &> /dev/null; then
    cat <<EOF

ERROR: Docker is not installed.

Install Docker Desktop for Mac:
  https://www.docker.com/products/docker-desktop/

After install, open Docker Desktop once to grant permissions and let it start.
Then re-run this script.

EOF
    exit 1
fi

if ! docker info &> /dev/null; then
    echo "ERROR: Docker is installed but daemon is not running. Start Docker Desktop and re-run."
    exit 1
fi

echo "  ✓ Docker daemon running ($(docker --version))"

# -----------------------------------------------------------------------------
# 4. Build the MAGMA Docker image (only if magma_binary/magma is in place)
# -----------------------------------------------------------------------------
echo ""
echo "=== Docker image build ==="
if [[ -f "magma_binary/magma" ]]; then
    if docker image inspect magma_celltype:latest >/dev/null 2>&1; then
        echo "Image magma_celltype:latest already exists. Rebuild manually if needed:"
        echo "  docker build -t magma_celltype:latest ."
    else
        echo "Building magma_celltype:latest..."
        docker build -t magma_celltype:latest .
        echo ""
        echo "Smoke-testing the binary inside the container:"
        ./run_in_docker.sh /usr/local/bin/magma --version
    fi
else
    echo "SKIPPING: magma_binary/magma not yet downloaded."
    echo "Download MAGMA v1.10 Linux binary, place at magma_binary/magma, then run:"
    echo "  docker build -t magma_celltype:latest ."
fi

# -----------------------------------------------------------------------------
# 5. Manual downloads remaining
# -----------------------------------------------------------------------------
cat <<EOF

=============================================================================
MANUAL DOWNLOADS STILL REQUIRED
=============================================================================

1) MAGMA v1.10 Linux binary (NOT the Mac version)
   • URL: https://cncr.nl/research/magma/
   • Pick: 'Static binary for Linux' (~500 KB)
   • Unzip to: $PROJECT_ROOT/magma_binary/
   • So you have: magma_binary/magma  (single file, executable)
   • Then build the Docker image:
       cd $PROJECT_ROOT
       docker build -t magma_celltype:latest .

2) MAGMA auxiliary files (1000G EUR reference panel)
   • URL: https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe
   • Need: g1000_eur.bed, g1000_eur.bim, g1000_eur.fam, g1000_eur.synonyms
   • Save to: $PROJECT_ROOT/auxfiles/
   • Size: ~1 GB total

3) Siletti specificity matrix (precomputed)
   • URL: https://www.dropbox.com/scl/fi/3p5qoyfw5c3q8yf38s0di/conti_specificity_matrix.txt
   • Save as: $PROJECT_ROOT/gene-level/Siletti_l2_conti_specificity_matrix.txt
   • Size: ~100 MB

4) Watanabe insomnia UKB sumstats
   • URL: https://vu.data.surfsara.nl/index.php/s/06RsHECyWqlBRwq
   • Get the UKB-only file (large; NOT the meta with 23andMe top-10k)
   • Save into: $PROJECT_ROOT/gwas_sumstats/insomnia_watanabe2022_ukb/
   • Update 'raw_filename' in pipeline/00_config.yaml

=============================================================================
NEXT STEPS (after downloads)
=============================================================================

# Build the Docker image
docker build -t magma_celltype:latest .

# Open Cursor at the project root
cursor $PROJECT_ROOT

# In Cursor's integrated terminal:
#   pip install pandas numpy scipy matplotlib pyyaml ipykernel
# (these are for running the notebook on the host;
#  the container has its own copies for reproducibility)

# Open pipeline/01_prep_sumstats.py
# Run cells 1-4 to inspect the Watanabe file
# Update 00_config.yaml columns block, re-run from cell 1

EOF
