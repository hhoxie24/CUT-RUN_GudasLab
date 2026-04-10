#!/bin/zsh
set -eo pipefail

###############################################################################
# 04_peakcall_v2.sh
#
# SEACR peak calling on spike-in-normalized bedGraphs.
# MACS3 section included but commented out (enable as needed).
#
# Inputs:
#   - 09_peaks_seacr/bedgraph/*.mm10.spikeNorm.bedgraph (from 03_bw_bg_v2.sh)
#   - 03_align_mm10/*.mm10.mapq30.bam (for MACS3, if enabled)
#
# Requires: CNR conda env (auto-activated)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate CNR 2>/dev/null || conda activate "$ENV_ALIGN"

load_config "$SCRIPT_DIR/config.json"
detect_threads

mkdir -p "$BIGWIG_DIR" "$MACS3_DIR" "$SEACR_DIR" "$BEDGRAPH_DIR"

if [[ ! -f "$NORM_TSV" ]]; then
    echo "ERROR: $NORM_TSV not found. Run 02_align_spikein_v2.sh first."
    exit 1
fi

PIPELINE_LOG="$PROJECT/04_peakcall_v2_$(date +%Y%m%d_%H%M%S).log"

{
echo "======================================="
echo " 04_peakcall_v2.sh"
echo " Date: $(date)"
echo " SEACR script: $SEACR_SCRIPT"
echo " SEACR mode: $PARAM_SEACR_MODE"
echo " SEACR threshold: $PARAM_SEACR_THRESH"
echo " Using $THREADS threads"
echo "======================================="
echo

###############################################################################
# MACS3 peak calling (commented out — enable as needed)
###############################################################################
# echo "Running MACS3 peak calling..."
#
# for BAM in "$ALIGN_MM10"/*.mm10.mapq30.bam; do
#     BASENAME=$(basename "$BAM" .mm10.mapq30.bam)
#     TARGET=$(echo "$BASENAME" | cut -d'_' -f2)
#     REP=$(echo "$BASENAME" | cut -d'_' -f3)
#
#     if [[ "$TARGET" == "IgG" ]]; then continue; fi
#
#     CTRL_BAM=$(ls "$ALIGN_MM10"/WT_IgG_"${REP}"_*".mm10.mapq30.bam" 2>/dev/null | head -n1 || true)
#     if [[ -z "$CTRL_BAM" ]]; then
#         echo "WARNING: No IgG control for $BASENAME (rep $REP), skipping."
#         continue
#     fi
#
#     OUT_PREFIX="$MACS3_DIR/${BASENAME}"
#
#     if [[ "$TARGET" == "H3K27me3" ]]; then
#         macs3 callpeak -t "$BAM" -c "$CTRL_BAM" -f BAMPE -g mm \
#             --keep-dup all --broad --broad-cutoff 0.1 -n "$OUT_PREFIX"
#     else
#         macs3 callpeak -t "$BAM" -c "$CTRL_BAM" -f BAMPE -g mm \
#             --keep-dup all -n "$OUT_PREFIX"
#     fi
#     echo
# done

###############################################################################
# SEACR peak calling on spike-in-normalized bedGraphs
###############################################################################
echo "Running SEACR peak calling..."

for BG in "$BEDGRAPH_DIR"/*.mm10.spikeNorm.bedgraph; do
    BASENAME=$(basename "$BG" .mm10.spikeNorm.bedgraph)

    TARGET=$(echo "$BASENAME" | cut -d'_' -f2)
    REP=$(echo "$BASENAME" | cut -d'_' -f3)

    # Skip IgG as treatment
    if [[ "$TARGET" == "IgG" ]]; then continue; fi

    CTRL_BG=$(ls "$BEDGRAPH_DIR"/WT_IgG_"${REP}"_*".mm10.spikeNorm.bedgraph" 2>/dev/null | head -n1 || true)

    if [[ -z "$CTRL_BG" ]]; then
        echo "WARNING: No IgG bedGraph for $BASENAME (rep $REP), skipping SEACR."
        continue
    fi

    OUT_PREFIX="$SEACR_DIR/${BASENAME}_SEACR"

    echo "SEACR treatment: $BG"
    echo "      control:   $CTRL_BG"
    echo "      prefix:    $OUT_PREFIX"

    # "non" = no internal normalization (inputs are pre-normalized by spike-in)
    bash "$SEACR_SCRIPT" \
        "$BG" \
        "$CTRL_BG" \
        "$PARAM_SEACR_MODE" \
        "$PARAM_SEACR_THRESH" \
        "$OUT_PREFIX"

    echo
done

echo "SEACR peak calling complete."
echo "=== 04_peakcall_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
