#!/bin/zsh
set -eo pipefail

###############################################################################
# 03_bw_bg_v2.sh
#
# Generates spike-in-normalized BigWigs and SEACR-compatible bedGraphs
# from NON-DEDUPLICATED MAPQ30 BAMs (duplicates retained for target samples).
# Also generates CPM-normalized BigWigs.
#
# Inputs (from 02_align_spikein_v2.sh):
#   - 03_align_mm10/*.mm10.mapq30.bam
#   - 05_spikein/spikein_normalization.tsv
#
# Outputs:
#   - 07_bigwig/<SAMPLE>.mm10.spikeNorm.bw
#   - 09_peaks_seacr/bedgraph/<SAMPLE>.mm10.spikeNorm.bedgraph
#   - 07_bigwig_cpm/<SAMPLE>.cpm.bw
#
# Requires: DEEPTOOLS conda env (auto-activated)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate DEEPTOOLS 2>/dev/null || conda activate "$ENV_DEEPTOOLS"

load_config "$SCRIPT_DIR/config.json"
detect_threads

mkdir -p "$BIGWIG_DIR" "$SEACR_DIR" "$BEDGRAPH_DIR"

if [[ ! -f "$NORM_TSV" ]]; then
    echo "ERROR: $NORM_TSV not found. Run 02_align_spikein_v2.sh first."
    exit 1
fi

PIPELINE_LOG="$PROJECT/03_bw_bg_v2_$(date +%Y%m%d_%H%M%S).log"

{
echo "======================================="
echo " 03_bw_bg_v2.sh"
echo " Date: $(date)"
echo " Spike-in table: $NORM_TSV"
echo " BAM source: $ALIGN_MM10 (mapq30, dups retained)"
echo " Using $THREADS threads"
echo "======================================="
echo

# Generate genome sizes file from BAM headers
if [[ ! -f "$GENOME_SIZES" ]]; then
    echo "Generating chromosome sizes from first BAM header..."
    FIRST_BAM=$(ls "$ALIGN_MM10"/*.mm10.mapq30.bam | head -1)
    samtools view -H "$FIRST_BAM" \
        | grep "^@SQ" \
        | awk -F'\t' '{gsub(/SN:/, "", $2); gsub(/LN:/, "", $3); print $2"\t"$3}' \
        > "$GENOME_SIZES"
    echo "Wrote: $GENOME_SIZES ($(wc -l < "$GENOME_SIZES") contigs)"
fi

###############################################################################
# Spike-in-normalized BigWigs (deeptools bamCoverage)
###############################################################################
echo "Generating spike-in-normalized BigWigs (non-deduplicated BAMs)..."

tail -n +2 "$NORM_TSV" | while read SAMPLE MM10_UNIQ ECOLI_UNIQ SPIKE_FRAC SCALE; do
    BAM="$ALIGN_MM10/${SAMPLE}.mm10.mapq30.bam"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping BigWig."
        continue
    fi

    BW_OUT="$BIGWIG_DIR/${SAMPLE}.mm10.spikeNorm.bw"

    echo "Sample: $SAMPLE  scaleFactor: $SCALE"

    bamCoverage \
        -b "$BAM" \
        -o "$BW_OUT" \
        --normalizeUsing None \
        --scaleFactor "$SCALE" \
        --binSize "$PARAM_BW_BIN" \
        --extendReads \
        -p $THREADS

    echo
done

echo "Spike-in BigWigs generated."
echo

###############################################################################
# SEACR-compatible bedGraphs (bedtools fragment method)
###############################################################################
echo "Generating SEACR-compatible fragment bedGraphs..."

tail -n +2 "$NORM_TSV" | while read SAMPLE MM10_UNIQ ECOLI_UNIQ SPIKE_FRAC SCALE; do
    BAM="$ALIGN_MM10/${SAMPLE}.mm10.mapq30.bam"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping bedGraph."
        continue
    fi

    BG_OUT="$BEDGRAPH_DIR/${SAMPLE}.mm10.spikeNorm.bedgraph"

    echo "Sample: $SAMPLE  scaleFactor: $SCALE"

    # Name-sort → fragment BED → coverage bedGraph → spike-in scaling
    samtools sort -n -@ $THREADS "$BAM" \
    | bedtools bamtobed -bedpe -i stdin \
    | awk -v maxfrag="$PARAM_FRAG_MAX" '$1==$4 && $6-$2 < maxfrag {print $0}' \
    | cut -f 1,2,6 \
    | sort -k1,1 -k2,2n -k3,3n \
    > "$BEDGRAPH_DIR/${SAMPLE}.fragments.bed"

    bedtools genomecov \
        -bg \
        -i "$BEDGRAPH_DIR/${SAMPLE}.fragments.bed" \
        -g "$GENOME_SIZES" \
    | awk -v sf="$SCALE" 'BEGIN{OFS="\t"} {$4 = $4 * sf; print}' \
    > "$BG_OUT"

    rm -f "$BEDGRAPH_DIR/${SAMPLE}.fragments.bed"
    echo
done

echo "SEACR bedGraphs generated."
echo

###############################################################################
# CPM-normalized BigWigs
###############################################################################
CPM_DIR="$PROJECT/07_bigwig_cpm"
mkdir -p "$CPM_DIR"

echo "Generating CPM-normalized BigWigs..."

for BAM in "$ALIGN_MM10"/*.mm10.mapq30.bam; do
    BASENAME=$(basename "$BAM" .mm10.mapq30.bam)
    BW="$CPM_DIR/${BASENAME}.cpm.bw"

    echo "Processing: $BASENAME"

    bamCoverage \
        --bam "$BAM" \
        --outFileName "$BW" \
        --binSize "$PARAM_BW_BIN" \
        --normalizeUsing CPM \
        --ignoreForNormalization chrX chrY chrM \
        --extendReads \
        -p $THREADS
done

echo "CPM tracks generated in $CPM_DIR"
echo
echo "=== 03_bw_bg_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
