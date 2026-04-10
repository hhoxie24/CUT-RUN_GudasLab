#!/bin/zsh
set -eo pipefail

###############################################################################
# 00_setup_conda_envs.sh
#
# Creates two conda environments for the CUT&RUN pipeline.
# Idempotent: safe to re-run; skips creation if env already exists.
#
# ENVIRONMENT ARCHITECTURE:
#   CNR:       bowtie2, samtools, bedtools, picard, trim-galore, cutadapt,
#              macs3, multiqc, openjdk 17, python 3.12
#   DEEPTOOLS: deeptools, samtools, bedtools, python 3.12
#
# Two envs avoid dependency conflicts on arm64:
#   1. Picard (Java) + R: libgfortran @rpath resolution broken on arm64
#   2. deeptools bundles pysam/htslib: conflicts with Picard JDK requirements
#
# USAGE:
#   zsh 00_setup_conda_envs.sh 2>&1 | tee conda_setup.log
#
# PREREQUISITES:
#   miniforge: brew install --cask miniforge && conda init zsh
###############################################################################

echo "======================================="
echo " CUT&RUN Pipeline: Conda Environment Setup"
echo " Date: $(date)"
echo " Architecture: $(uname -m)"
echo "======================================="
echo

# Verify conda is available and is miniforge (not Anaconda)
CONDA_ROOT="$(conda info --base 2>/dev/null || echo '')"
if [[ -z "$CONDA_ROOT" ]]; then
    echo "ERROR: conda not found. Install miniforge first."
    exit 1
fi
if [[ "$CONDA_ROOT" == *"anaconda"* ]]; then
    echo "ERROR: Active conda is Anaconda ($CONDA_ROOT). Switch to miniforge."
    exit 1
fi
echo "Conda root: $CONDA_ROOT"

SOLVER="conda"
command -v mamba &>/dev/null && SOLVER="mamba"
echo "Solver: $SOLVER"
echo

setup_env() {
    local ENV_NAME="$1"; shift
    if conda env list 2>/dev/null | grep -qw "$ENV_NAME"; then
        echo "Environment '$ENV_NAME' already exists."
        conda list -n "$ENV_NAME" 2>/dev/null | grep -v "^#"
        return 0
    fi
    echo "Creating environment '$ENV_NAME'..."
    $SOLVER create -n "$ENV_NAME" -y "$@"
    echo "'$ENV_NAME' created."
}

echo "=== ENV 1: CNR ==="
setup_env CNR \
    -c conda-forge -c bioconda \
    python=3.12 bowtie2 samtools bedtools picard \
    trim-galore cutadapt macs3 multiqc openjdk=17

echo
echo "=== ENV 2: DEEPTOOLS ==="
setup_env DEEPTOOLS \
    -c conda-forge -c bioconda \
    python=3.12 deeptools samtools bedtools

echo
echo "=== Verifying installations ==="

eval "$(conda shell.zsh hook)"

echo "--- CNR env ---"
conda activate CNR
for TOOL in bowtie2 samtools bedtools picard trim_galore cutadapt macs3 python java; do
    command -v "$TOOL" &>/dev/null && echo "  $TOOL: OK" || echo "  $TOOL: NOT FOUND"
done
conda deactivate

echo "--- DEEPTOOLS env ---"
conda activate DEEPTOOLS
for TOOL in bamCoverage samtools bedtools python; do
    command -v "$TOOL" &>/dev/null && echo "  $TOOL: OK" || echo "  $TOOL: NOT FOUND"
done
conda deactivate

echo
echo "Setup complete. Each pipeline script self-activates its conda env."
