#!/bin/zsh
set -eo pipefail

###############################################################################
# 05_merge_annotate_v2.sh
#
# Merge MACS3 peaks across replicates per target using HOMER mergePeaks,
# then annotate merged peaks with HOMER annotatePeaks.pl.
#
# Requires: CNR conda env + HOMER
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate CNR 2>/dev/null || conda activate "$ENV_ALIGN"

load_config "$SCRIPT_DIR/config.json"

export PATH="$HOMER_BIN:$HOMER_CPP:$PATH"

OUT_DIR="$MACS3_DIR/HomerMergePeaks"
ANO_DIR="$OUT_DIR/AnnoPeaks"

PIPELINE_LOG="$PROJECT/05_merge_annotate_v2_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR" "$ANO_DIR"

{
echo "========================================================="
echo " 05_merge_annotate_v2.sh"
echo " Date: $(date)"
echo " Input: $MACS3_DIR"
echo " HOMER: $HOMER_BIN"
echo "========================================================="

if ! command -v mergePeaks &> /dev/null; then
    echo "ERROR: HOMER tools not found in PATH."
    exit 1
fi

TARGETS=($(ls "$MACS3_DIR"/*_peaks.*Peak 2>/dev/null | xargs -n 1 basename | cut -d'_' -f2 | sort | uniq | grep -v "IgG"))

for TARGET in "${TARGETS[@]}"; do
    echo
    echo "Processing Target: $TARGET"
    echo "---------------------------------------------------------"

    REPLICATE_FILES=("$MACS3_DIR"/WT_${TARGET}_*_peaks.*Peak)

    if [ ${#REPLICATE_FILES[@]} -eq 0 ]; then
        echo "WARNING: No files found for $TARGET. Skipping."
        continue
    fi

    echo "  Found ${#REPLICATE_FILES[@]} replicates"

    MERGED_TXT="$OUT_DIR/WT_${TARGET}_merged.txt"
    mergePeaks -d given "${REPLICATE_FILES[@]}" > "$MERGED_TXT"

    ANNOTATED_TXT="$ANO_DIR/WT_${TARGET}_annotated.txt"
    perl "$HOMER_BIN/annotatePeaks.pl" "$MERGED_TXT" mm10 > "$ANNOTATED_TXT"

    FINAL_BED="$OUT_DIR/WT_${TARGET}_merged.bed"
    awk 'NR>1 {print $2, $3, $4, $1, "1000", "."}' OFS="\t" "$MERGED_TXT" > "$FINAL_BED"

    echo "  Generated: $FINAL_BED"
    echo "  Generated: $ANNOTATED_TXT"
done

echo
echo "=== 05_merge_annotate_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
