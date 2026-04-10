#!/bin/zsh
###############################################################################
# parse_config.sh
#
# Sourced by all pipeline scripts to load paths from config.json.
# Requires: python3 (available in both CNR and DEEPTOOLS envs)
#
# Usage (in pipeline scripts):
#   SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
#   source "$SCRIPT_DIR/scripts/parse_config.sh"
#   load_config "$SCRIPT_DIR/config.json"
#
# After load_config, the following variables are set:
#   PROJECT, FASTQ_DIR, BT2_MM10, BT2_ECOLI,
#   HOMER_BIN, HOMER_CPP, SEACR_SCRIPT,
#   ENV_ALIGN, ENV_DEEPTOOLS,
#   PARAM_MAX_INSERT, PARAM_MAPQ, PARAM_FRAG_MAX, PARAM_BW_BIN,
#   PARAM_TG_CORES, PARAM_SEACR_MODE, PARAM_SEACR_THRESH
###############################################################################

load_config() {
    local CONFIG_FILE="$1"

    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "ERROR: config.json not found at $CONFIG_FILE"
        echo "       Copy config.json.example to config.json and edit paths."
        exit 1
    fi

    # Parse JSON with python (available in conda envs; no jq dependency)
    eval "$(python3 - "$CONFIG_FILE" << 'PYEOF'
import json, sys

with open(sys.argv[1]) as f:
    cfg = json.load(f)

proj = cfg["project_dir"]
fq   = cfg.get("fastq_dir", "02_fastq")

print(f'PROJECT="{proj}"')
print(f'FASTQ_DIR="{proj}/{fq}"')
print(f'BT2_MM10="{cfg["genome"]["bowtie2_index_mm10"]}"')
print(f'BT2_ECOLI="{cfg["genome"]["bowtie2_index_ecoli"]}"')
print(f'HOMER_BIN="{cfg["tools"]["homer_bin"]}"')
print(f'HOMER_CPP="{cfg["tools"].get("homer_cpp", cfg["tools"]["homer_bin"].replace("/bin", "/cpp"))}"')
print(f'SEACR_SCRIPT="{cfg["tools"]["seacr_script"]}"')
print(f'ENV_ALIGN="{cfg["conda_envs"]["alignment"]}"')
print(f'ENV_DEEPTOOLS="{cfg["conda_envs"]["deeptools"]}"')

p = cfg.get("parameters", {})
print(f'PARAM_MAX_INSERT="{p.get("bowtie2_max_insert", 700)}"')
print(f'PARAM_MAPQ="{p.get("mapq_threshold", 30)}"')
print(f'PARAM_FRAG_MAX="{p.get("fragment_max_size", 1000)}"')
print(f'PARAM_BW_BIN="{p.get("bigwig_bin_size", 10)}"')
print(f'PARAM_TG_CORES="{p.get("trim_galore_cores", 4)}"')
print(f'PARAM_SEACR_MODE="{p.get("seacr_mode", "non")}"')
print(f'PARAM_SEACR_THRESH="{p.get("seacr_threshold", "stringent")}"')
PYEOF
)"

    # Derived output directories (relative to PROJECT)
    TRIM_DIR="$PROJECT/02b_trimmed"
    ALIGN_MM10="$PROJECT/03_align_mm10"
    ALIGN_ECOLI="$PROJECT/04_align_ecoli"
    SPIKE_DIR="$PROJECT/05_spikein"
    QC_DIR="$PROJECT/06_fragment_qc"
    BIGWIG_DIR="$PROJECT/07_bigwig"
    MACS3_DIR="$PROJECT/08_peaks_macs3"
    SEACR_DIR="$PROJECT/09_peaks_seacr"
    BEDGRAPH_DIR="$SEACR_DIR/bedgraph"

    LOG_MM10="$ALIGN_MM10/logs"
    LOG_ECOLI="$ALIGN_ECOLI/logs"
    LOG_TRIM="$TRIM_DIR/logs"

    NORM_TSV="$SPIKE_DIR/spikein_normalization.tsv"
    GENOME_SIZES="$ALIGN_MM10/mm10.chrom.sizes"

    # Validate critical paths exist
    local ERRORS=0
    for VAR_NAME VAR_VAL in \
        "project_dir" "$PROJECT" \
        "bowtie2_index_mm10" "$BT2_MM10" \
        "bowtie2_index_ecoli" "$BT2_ECOLI"; do
        # Check if the value still contains "path/to"
        if [[ "$VAR_VAL" == *"/path/to/"* ]]; then
            echo "ERROR: $VAR_NAME is not configured in config.json (still contains placeholder path)"
            ERRORS=$((ERRORS + 1))
        fi
    done

    if [[ $ERRORS -gt 0 ]]; then
        echo "Fix the above path(s) in config.json before running."
        exit 1
    fi
}

# Thread detection for Apple Silicon (P-cores only) or fallback
detect_threads() {
    if sysctl -n hw.perflevel0.logicalcpu &>/dev/null; then
        THREADS=$(sysctl -n hw.perflevel0.logicalcpu)
    elif command -v nproc &>/dev/null; then
        THREADS=$(nproc)
    else
        THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
    fi
}
