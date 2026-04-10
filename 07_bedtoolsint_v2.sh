#!/bin/zsh
set -eo pipefail

###############################################################################
# 07_bedtoolsint_v2.sh
#
# Pairwise bedtools intersect between RARg replicates,
# then HOMER annotation of the combined intersection file.
#
# Requires: CNR conda env + HOMER
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate CNR 2>/dev/null || conda activate "$ENV_ALIGN"

load_config "$SCRIPT_DIR/config.json"

export PATH="$HOMER_BIN:$PATH"

OUT_DIR="$MACS3_DIR/Intersect"

PIPELINE_LOG="$PROJECT/07_bedtoolsint_v2_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

{
echo "========================================================="
echo " 07_bedtoolsint_v2.sh"
echo " Date: $(date)"
echo " Input: $MACS3_DIR"
echo " Output: $OUT_DIR"
echo "========================================================="

# Pairwise intersections between RARg replicates
# NOTE: Adjust filenames below to match your MACS3 output naming
bedtools intersect \
    -a "$MACS3_DIR"/WT_RARg_2_*_peaks.narrowPeak \
    -b "$MACS3_DIR"/WT_RARg_1_*_peaks.narrowPeak \
    > "$OUT_DIR/RARg_2x1.txt"

bedtools intersect \
    -a "$MACS3_DIR"/WT_RARg_2_*_peaks.narrowPeak \
    -b "$MACS3_DIR"/WT_RARg_3_*_peaks.narrowPeak \
    > "$OUT_DIR/RARg_2x3.txt"

bedtools intersect \
    -a "$MACS3_DIR"/WT_RARg_1_*_peaks.narrowPeak \
    -b "$MACS3_DIR"/WT_RARg_3_*_peaks.narrowPeak \
    > "$OUT_DIR/RARg_1x3.txt"

echo "  RARg_2x1: $(wc -l < "$OUT_DIR/RARg_2x1.txt") peaks"
echo "  RARg_2x3: $(wc -l < "$OUT_DIR/RARg_2x3.txt") peaks"
echo "  RARg_1x3: $(wc -l < "$OUT_DIR/RARg_1x3.txt") peaks"

# Combine and annotate
cat "$OUT_DIR"/RARg_*.txt | sort -k1,1 -k2,2n | uniq > "$OUT_DIR/RARg_all.txt"

perl "$HOMER_BIN/annotatePeaks.pl" \
    "$OUT_DIR/RARg_all.txt" mm10 \
    > "$OUT_DIR/RARg_all_anno.txt"

echo "  Annotated: RARg_all_anno.txt"
echo
echo "=== 07_bedtoolsint_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
