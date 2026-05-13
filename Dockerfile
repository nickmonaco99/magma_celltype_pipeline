# =============================================================================
# Dockerfile — minimal ubuntu:22.04 with MAGMA v1.10 for the cell-type pipeline
# =============================================================================
# Why this exists:
#   - Apple Silicon (M3 Pro) doesn't run MAGMA's Mac binary natively
#   - Rosetta works but has intermittent issues with MAGMA's mmap'd reads on the
#     large 1000G reference panel (.bed file is ~500MB)
#   - This container gives bit-identical Linux behaviour to Duncan's environment
#     (their Example_results logs show MAGMA v1.10 on Linux)
#
# Build:
#   docker build -t magma_celltype:latest .
#
# Run (the project root gets mounted at /work):
#   docker run --rm -v "$(pwd):/work" -w /work magma_celltype:latest \
#     /usr/local/bin/magma --version
#
# Or use the wrapper script:
#   ./run_in_docker.sh /usr/local/bin/magma --version
# =============================================================================

FROM ubuntu:22.04

# Avoid prompts during apt operations
ENV DEBIAN_FRONTEND=noninteractive

# Minimal runtime deps. MAGMA is statically linked but needs libc6 etc.
# We also install Python + a few stats packages so QC notebooks can run inside
# the container if you ever want full reproducibility (host Python differences
# can be a subtle reproducibility hazard).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libc6 \
        libgomp1 \
        wget \
        unzip \
        python3 \
        python3-pip \
        r-base \
    && rm -rf /var/lib/apt/lists/*

# Python deps for the QC notebooks
RUN pip3 install --no-cache-dir \
        pandas==2.2.* \
        numpy==1.26.* \
        scipy==1.13.* \
        matplotlib==3.8.* \
        pyyaml==6.* \
        jupyterlab==4.*

# R deps for Duncan's scripts 3 and 5
RUN R -e 'install.packages(c("readr", "dplyr", "tidyr"), repos="https://cloud.r-project.org/")'

# Install MAGMA v1.10 from the binary you downloaded into magma_binary/
# (we COPY rather than wget so the build is fully self-contained and doesn't
# depend on the MAGMA website staying up)
COPY magma_binary/magma /usr/local/bin/magma
RUN chmod +x /usr/local/bin/magma

# Sanity-check the binary works during build (fails build if you copied a
# Mac binary by mistake — Linux can't exec it)
RUN /usr/local/bin/magma --version 2>&1 | grep -q "v1.10" || \
    (echo "ERROR: MAGMA binary is not v1.10 or not a Linux binary"; exit 1)

# Project root will be mounted here
WORKDIR /work

# Default to bash so `docker run -it` drops into a shell; the wrapper script
# passes a specific command which overrides this.
CMD ["/bin/bash"]
