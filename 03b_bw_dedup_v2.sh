#!/bin/zsh
set -eo pipefail

###############################################################################
# 03b_bw_dedup_v2.sh
#
# Optional: generates deduplicated BigWigs + bedGraphs.
# Reads existing mapq30.dedup.bam files and spikein_normalization.tsv.
# No redundant alignment or duplicate-marking steps.
#
# Outputs:
#   - 07_bigwig/dedup/<SAMPLE>.mm10.dedup.spikeNorm.bw
#   - 09_peaks_seacr/bedgraph/dedup/<SAMPLE>.mm10.dedup.spikeNorm.bedgraph
#   - 07_bigwig_cpm/dedup/<SAMPLE>.dedup.cpm.bw
#
# Requires: DEEPTOOLS conda env (auto-activated)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate DEEPTOOLS 2>/dev/null || conda activate "$ENV_DEEPTOOLS"

load_config "$SCRIPT_DIR/config.json"
detect_threads

DEDUP_BW_DIR="$PROJECT/07_bigwig/dedup"
DEDUP_BG_DIR="$SEACR_DIR/bedgraph/dedup"
DEDUP_CPM_DIR="$PROJECT/07_bigwig_cpm/dedup"

mkdir -p "$DEDUP_BW_DIR" "$DEDUP_BG_DIR" "$DEDUP_CPM_DIR"

if [[ ! -f "$NORM_TSV" ]]; then
    echo "ERROR: $NORM_TSV not found. Run 02_align_spikein_v2.sh first."
    exit 1
fi

PIPELINE_LOG="$PROJECT/03b_bw_dedup_v2_$(date +%Y%m%d_%H%M%S).log"

{
echo "======================================="
echo " 03b_bw_dedup_v2.sh"
echo " Date: $(date)"
echo " BAM source: $ALIGN_MM10 (mapq30.dedup)"
echo " Using $THREADS threads"
echo "======================================="
echo

# Genome sizes file
if [[ ! -f "$GENOME_SIZES" ]]; then
    FIRST_BAM=$(ls "$ALIGN_MM10"/*.mm10.mapq30.dedup.bam | head -1)
    samtools view -H "$FIRST_BAM" \
        | grep "^@SQ" \
        | awk -F'\t' '{gsub(/SN:/, "", $2); gsub(/LN:/, "", $3); print $2"\t"$3}' \
        > "$GENOME_SIZES"
fi

# Spike-in BigWigs (dedup)
echo "Generating spike-in-normalized BigWigs (deduplicated BAMs)..."
tail -n +2 "$NORM_TSV" | while read SAMPLE MM10_UNIQ ECOLI_UNIQ SPIKE_FRAC SCALE; do
    BAM="$ALIGN_MM10/${SAMPLE}.mm10.mapq30.dedup.bam"
    [[ ! -f "$BAM" ]] && { echo "WARNING: $BAM not found, skipping."; continue; }

    bamCoverage -b "$BAM" -o "$DEDUP_BW_DIR/${SAMPLE}.mm10.dedup.spikeNorm.bw" \
        --normalizeUsing None --scaleFactor "$SCALE" \
        --binSize "$PARAM_BW_BIN" --extendReads -p $THREADS
    echo "  $SAMPLE (sf=$SCALE)"
done
echo

# SEACR bedGraphs (dedup)
echo "Generating SEACR fragment bedGraphs (deduplicated)..."
tail -n +2 "$NORM_TSV" | while read SAMPLE MM10_UNIQ ECOLI_UNIQ SPIKE_FRAC SCALE; do
    BAM="$ALIGN_MM10/${SAMPLE}.mm10.mapq30.dedup.bam"
    [[ ! -f "$BAM" ]] && continue

    samtools sort -n -@ $THREADS "$BAM" \
    | bedtools bamtobed -bedpe -i stdin \
    | awk -v maxfrag="$PARAM_FRAG_MAX" '$1==$4 && $6-$2 < maxfrag {print $0}' \
    | cut -f 1,2,6 \
    | sort -k1,1 -k2,2n -k3,3n \
    > "$DEDUP_BG_DIR/${SAMPLE}.dedup.fragments.bed"

    bedtools genomecov -bg \
        -i "$DEDUP_BG_DIR/${SAMPLE}.dedup.fragments.bed" \
        -g "$GENOME_SIZES" \
    | awk -v sf="$SCALE" 'BEGIN{OFS="\t"} {$4 = $4 * sf; print}' \
    > "$DEDUP_BG_DIR/${SAMPLE}.mm10.dedup.spikeNorm.bedgraph"

    rm -f "$DEDUP_BG_DIR/${SAMPLE}.dedup.fragments.bed"
    echo "  $SAMPLE"
done
echo

# CPM BigWigs (dedup)
echo "Generating CPM BigWigs (deduplicated)..."
for BAM in "$ALIGN_MM10"/*.mm10.mapq30.dedup.bam; do
    BASENAME=$(basename "$BAM" .mm10.mapq30.dedup.bam)
    bamCoverage --bam "$BAM" --outFileName "$DEDUP_CPM_DIR/${BASENAME}.dedup.cpm.bw" \
        --binSize "$PARAM_BW_BIN" --normalizeUsing CPM \
        --ignoreForNormalization chrX chrY chrM --extendReads -p $THREADS
    echo "  $BASENAME"
done

echo
echo "=== 03b_bw_dedup_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
