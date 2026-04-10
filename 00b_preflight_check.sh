#!/bin/zsh
set -eo pipefail

###############################################################################
# 00b_preflight_check.sh
#
# Validates: conda envs, tools, genome indices, FASTQs, HOMER, SEACR, disk.
# Run BEFORE starting the pipeline (script 02 onwards).
#
# USAGE:
#   zsh 00b_preflight_check.sh 2>&1 | tee preflight.log
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

# Load config before conda activation (parse_config uses python from PATH or system)
# We'll activate conda envs during tool checks below

PASS=0; FAIL=0; WARN=0
pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }
warn() { echo "  WARN: $1"; WARN=$((WARN + 1)); }

echo "======================================="
echo " CUT&RUN Pipeline: Pre-flight Check"
echo " Date: $(date)"
echo "======================================="
echo

# 1. Conda
echo "--- 1. Conda environment ---"
CONDA_ROOT="$(conda info --base 2>/dev/null || echo '')"
if [[ -z "$CONDA_ROOT" ]]; then
    fail "conda not found in PATH"
elif [[ "$CONDA_ROOT" == *"anaconda"* ]]; then
    fail "Active conda is Anaconda — should be miniforge"
else
    pass "Conda root: $CONDA_ROOT"
fi

eval "$(conda shell.zsh hook)"

# Need conda active to run python for config parsing
conda activate CNR 2>/dev/null || conda activate base
load_config "$SCRIPT_DIR/config.json"
conda deactivate 2>/dev/null

echo
echo " Project: $PROJECT"
echo

for ENV in CNR DEEPTOOLS; do
    conda env list 2>/dev/null | grep -qw "$ENV" && pass "Environment '$ENV' exists" || fail "Environment '$ENV' not found"
done
echo

# 2. CNR tools
echo "--- 2. CNR env tools ---"
conda activate CNR
for TOOL in bowtie2 samtools bedtools picard trim_galore cutadapt macs3 multiqc java python; do
    command -v "$TOOL" &>/dev/null && pass "$TOOL: found" || fail "$TOOL: not found"
done
conda deactivate
echo

# 3. DEEPTOOLS tools
echo "--- 3. DEEPTOOLS env tools ---"
conda activate DEEPTOOLS
for TOOL in bamCoverage samtools bedtools python; do
    command -v "$TOOL" &>/dev/null && pass "$TOOL: found" || fail "$TOOL: not found"
done
conda deactivate
echo

# 4. Genome indices
echo "--- 4. Reference genome indices ---"
for IDX_PREFIX in "$BT2_MM10" "$BT2_ECOLI"; do
    IDX_NAME=$(basename "$IDX_PREFIX")
    MISSING=0
    for EXT in .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2; do
        [[ -f "${IDX_PREFIX}${EXT}" ]] || [[ -f "${IDX_PREFIX}${EXT}l" ]] || MISSING=$((MISSING + 1))
    done
    [[ $MISSING -eq 0 ]] && pass "Bowtie2 index: $IDX_NAME" || fail "Bowtie2 index: $IDX_NAME ($MISSING files missing)"
done
echo

# 5. FASTQs
echo "--- 5. FASTQ inputs ---"
if [[ ! -d "$FASTQ_DIR" ]]; then
    fail "FASTQ directory not found: $FASTQ_DIR"
else
    R1_COUNT=$(ls "$FASTQ_DIR"/*_R1_001.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$R1_COUNT" -eq 0 ]]; then
        fail "No *_R1_001.fastq.gz files found"
    else
        pass "Found $R1_COUNT samples in $FASTQ_DIR"
        UNPAIRED=0
        for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
            SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
            [[ ! -f "$FASTQ_DIR/${SAMPLE}_R2_001.fastq.gz" ]] && { fail "Missing R2 for: $SAMPLE"; UNPAIRED=$((UNPAIRED + 1)); }
        done
        [[ $UNPAIRED -eq 0 ]] && pass "All samples have paired R1/R2 files"
    fi
fi
echo

# 6. HOMER
echo "--- 6. HOMER ---"
[[ -f "$HOMER_BIN/annotatePeaks.pl" ]] && pass "HOMER found at $HOMER_BIN" || warn "HOMER not found (needed for scripts 05-07)"
echo

# 7. SEACR
echo "--- 7. SEACR ---"
[[ -f "$SEACR_SCRIPT" ]] && pass "SEACR found: $SEACR_SCRIPT" || warn "SEACR not found (needed for script 04)"
echo

# 8. System resources
echo "--- 8. System resources ---"
detect_threads
pass "Threads: $THREADS"
if command -v sysctl &>/dev/null; then
    RAM_BYTES=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
    [[ $RAM_BYTES -gt 0 ]] && pass "RAM: $((RAM_BYTES / 1073741824)) GB"
fi
echo

# Summary
echo "==========================================="
echo " SUMMARY: PASS=$PASS  WARN=$WARN  FAIL=$FAIL"
echo "==========================================="
if [[ $FAIL -gt 0 ]]; then
    echo "  STATUS: BLOCKED — fix FAIL items before running pipeline."
    exit 1
elif [[ $WARN -gt 0 ]]; then
    echo "  STATUS: READY with warnings."
    exit 0
else
    echo "  STATUS: READY — all checks passed."
    exit 0
fi
