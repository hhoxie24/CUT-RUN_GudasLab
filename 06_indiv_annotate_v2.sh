#!/bin/zsh
set -eo pipefail

###############################################################################
# 06_indiv_annotate_v2.sh
#
# Annotates every individual MACS3 and SEACR peak file against mm10.
# Useful for checking replicate consistency before merging.
#
# Requires: CNR conda env + HOMER
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate CNR 2>/dev/null || conda activate "$ENV_ALIGN"

load_config "$SCRIPT_DIR/config.json"

export PATH="$HOMER_BIN:$HOMER_CPP:$PATH"

OUT_DIR_M="$MACS3_DIR/HomerAnno"
OUT_DIR_S="$SEACR_DIR/HomerAnno"

PIPELINE_LOG="$PROJECT/06_indiv_annotate_v2_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR_M" "$OUT_DIR_S"

{
echo "========================================================="
echo " 06_indiv_annotate_v2.sh"
echo " Date: $(date)"
echo "========================================================="

# MACS3 peak annotations
echo " Annotating Individual MACS3 Peaks"
echo " Input:  $MACS3_DIR"
echo " Output: $OUT_DIR_M"
echo "========================================================="

for PEAK_FILE in "$MACS3_DIR"/*_peaks.*Peak; do
    if [[ ! -f "$PEAK_FILE" ]]; then continue; fi
    FILENAME=$(basename "$PEAK_FILE")
    BASENAME="${FILENAME%.*}"
    echo "Processing: $FILENAME"
    annotatePeaks.pl "$PEAK_FILE" mm10 > "$OUT_DIR_M/${BASENAME}_annotated.txt"
done

echo
echo "Done. Individual MACS3 annotations in $OUT_DIR_M"

# SEACR peak annotations
echo
echo "========================================================="
echo " Annotating Individual SEACR Peaks"
echo " Input:  $SEACR_DIR"
echo " Output: $OUT_DIR_S"
echo "========================================================="

for PEAK_FILE in "$SEACR_DIR"/*_SEACR.stringent.bed; do
    if [[ ! -f "$PEAK_FILE" ]]; then continue; fi
    FILENAME=$(basename "$PEAK_FILE")
    BASENAME="${FILENAME%.*}"
    BASENAME="${BASENAME%.*}"

    echo "Processing: $FILENAME"

    CLEAN_BED="$OUT_DIR_S/${BASENAME}.clean.bed"
    awk -v pfx="$BASENAME" '{print $1, $2, $3, pfx"_"NR, $5, "."}' OFS="\t" "$PEAK_FILE" > "$CLEAN_BED"

    annotatePeaks.pl "$CLEAN_BED" mm10 > "$OUT_DIR_S/${BASENAME}_annotated.txt"
done

echo
echo "Done. Individual SEACR annotations in $OUT_DIR_S"
echo
echo "=== 06_indiv_annotate_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
