#!/bin/zsh
set -eo pipefail

###############################################################################
# 02_align_spikein_v2.sh
#
# CUT&RUN pipeline: Trim → mm10 alignment → MarkDup → QC → MAPQ filter → E. coli spike-in
#
# Produces (per sample):
#   - markdup BAM:      03_align_mm10/<SAMPLE>.mm10.markdup.bam
#   - MAPQ30 w/ dups:   03_align_mm10/<SAMPLE>.mm10.mapq30.bam
#   - MAPQ30 deduped:   03_align_mm10/<SAMPLE>.mm10.mapq30.dedup.bam
#   - E. coli BAM:      04_align_ecoli/<SAMPLE>.ecoli.sorted.bam
#   - Spike-in table:   05_spikein/spikein_normalization.tsv
#
# Requires: CNR conda env (auto-activated)
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/scripts/parse_config.sh"

eval "$(conda shell.zsh hook)"
conda activate CNR 2>/dev/null || conda activate "$ENV_ALIGN"

load_config "$SCRIPT_DIR/config.json"
detect_threads

mkdir -p "$ALIGN_MM10" "$ALIGN_ECOLI" "$LOG_MM10" "$LOG_ECOLI" \
         "$SPIKE_DIR" "$QC_DIR" "$TRIM_DIR" "$LOG_TRIM"

PIPELINE_LOG="$PROJECT/02_align_spikein_v2_$(date +%Y%m%d_%H%M%S).log"

{
echo "======================================="
echo " 02_align_spikein_v2.sh"
echo " Date: $(date)"
echo " Project: $PROJECT"
echo " FASTQ dir: $FASTQ_DIR"
echo " Trimmed dir: $TRIM_DIR"
echo " mm10 index: $BT2_MM10"
echo " E. coli index: $BT2_ECOLI"
echo " Using $THREADS threads"
echo "======================================="
echo

###############################################################################
# STEP 1: Trim, Align, Mark Duplicates, QC, Filter
###############################################################################
cd "$FASTQ_DIR"

for R1 in *_R1_001.fastq.gz; do
    SAMPLE=${R1%_R1_001.fastq.gz}
    R2="${SAMPLE}_R2_001.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: R2 missing for $SAMPLE, skipping."
        continue
    fi

    echo "==============================================="
    echo "Processing sample: $SAMPLE"
    echo "-----------------------------------------------"

    # Adapter trimming with Trim Galore
    TRIM_R1="$TRIM_DIR/${SAMPLE}_R1_001_val_1.fq.gz"
    TRIM_R2="$TRIM_DIR/${SAMPLE}_R2_001_val_2.fq.gz"

    if [[ -f "$TRIM_R1" && -f "$TRIM_R2" ]]; then
        echo "[Trim] Trimmed files exist for $SAMPLE, skipping trim step."
    else
        echo "[Trim] Running Trim Galore on $SAMPLE ..."
        trim_galore \
            --paired \
            --fastqc \
            --cores "$PARAM_TG_CORES" \
            --output_dir "$TRIM_DIR" \
            "$R1" "$R2" \
            2>&1 | tee "$LOG_TRIM/${SAMPLE}.trim_galore.log"
    fi

    # Align to mm10 (end-to-end, dovetail for CUT&RUN short fragments)
    echo "[mm10] Aligning $SAMPLE ..."

    bowtie2 \
        --end-to-end --very-sensitive --dovetail \
        --no-mixed --no-discordant --phred33 -I 10 -X "$PARAM_MAX_INSERT" \
        -x "$BT2_MM10" \
        -1 "$TRIM_R1" \
        -2 "$TRIM_R2" \
        -p $THREADS \
        2> "$LOG_MM10/${SAMPLE}.mm10.bowtie2.log" \
    | samtools view -bS - \
    | samtools sort -@ $THREADS -o "$ALIGN_MM10/${SAMPLE}.mm10.sorted.bam" -

    # Mark duplicates (flag only, do not remove)
    echo "[mm10] Marking duplicates for $SAMPLE ..."

    samtools sort -n -@ $THREADS \
        -o "$ALIGN_MM10/${SAMPLE}.mm10.namesort.bam" \
        "$ALIGN_MM10/${SAMPLE}.mm10.sorted.bam"

    samtools fixmate -m -@ $THREADS \
        "$ALIGN_MM10/${SAMPLE}.mm10.namesort.bam" \
        "$ALIGN_MM10/${SAMPLE}.mm10.fixmate.bam"

    samtools sort -@ $THREADS \
        -o "$ALIGN_MM10/${SAMPLE}.mm10.coordsort.bam" \
        "$ALIGN_MM10/${SAMPLE}.mm10.fixmate.bam"

    MARKDUP="$ALIGN_MM10/${SAMPLE}.mm10.markdup.bam"
    samtools markdup -@ $THREADS \
        "$ALIGN_MM10/${SAMPLE}.mm10.coordsort.bam" "$MARKDUP"
    samtools index "$MARKDUP"

    # Cleanup intermediate BAMs
    rm -f "$ALIGN_MM10/${SAMPLE}.mm10.namesort.bam" \
          "$ALIGN_MM10/${SAMPLE}.mm10.fixmate.bam" \
          "$ALIGN_MM10/${SAMPLE}.mm10.coordsort.bam" \
          "$ALIGN_MM10/${SAMPLE}.mm10.sorted.bam"

    # QC: samtools stats + Picard insert size metrics
    samtools stats "$MARKDUP" > "$LOG_MM10/${SAMPLE}.mm10.stats"

    picard CollectInsertSizeMetrics \
          I="$MARKDUP" \
          O="$LOG_MM10/${SAMPLE}.insert_size_metrics.txt" \
          H="$LOG_MM10/${SAMPLE}.insert_size_histogram.pdf" \
          M=0.5 \
    || {
        echo "WARNING: Picard histogram PDF failed for $SAMPLE (likely R/libgfortran issue on arm64)."
        echo "         Metrics .txt was still written."
        if [[ ! -s "$LOG_MM10/${SAMPLE}.insert_size_metrics.txt" ]]; then
            picard CollectInsertSizeMetrics \
                  I="$MARKDUP" \
                  O="$LOG_MM10/${SAMPLE}.insert_size_metrics.txt" \
                  H=/dev/null \
                  M=0.5 \
            || echo "ERROR: Picard metrics collection also failed for $SAMPLE."
        fi
    }

    # MAPQ-filtered output BAMs
    echo "[mm10] Creating MAPQ>=$PARAM_MAPQ BAMs for $SAMPLE ..."

    MAPQ30_BAM="$ALIGN_MM10/${SAMPLE}.mm10.mapq30.bam"
    MAPQ30_DEDUP_BAM="$ALIGN_MM10/${SAMPLE}.mm10.mapq30.dedup.bam"

    # Retain duplicates (target sample peak calling)
    samtools view -b -q "$PARAM_MAPQ" -f 2 "$MARKDUP" > "$MAPQ30_BAM"
    samtools index "$MAPQ30_BAM"

    # Remove duplicates (IgG control peak calling, visualization)
    samtools view -b -q "$PARAM_MAPQ" -f 2 -F 1024 "$MARKDUP" > "$MAPQ30_DEDUP_BAM"
    samtools index "$MAPQ30_DEDUP_BAM"

    # E. coli spike-in alignment
    echo "[Ecoli] Aligning $SAMPLE ..."

    bowtie2 \
        --end-to-end --very-sensitive --dovetail \
        --no-mixed --no-discordant --phred33 -I 10 -X "$PARAM_MAX_INSERT" \
        -x "$BT2_ECOLI" \
        -1 "$TRIM_R1" \
        -2 "$TRIM_R2" \
        -p $THREADS \
        2> "$LOG_ECOLI/${SAMPLE}.ecoli.bowtie2.log" \
    | samtools view -bS - \
    | samtools sort -@ $THREADS -o "$ALIGN_ECOLI/${SAMPLE}.ecoli.sorted.bam" -

    samtools index "$ALIGN_ECOLI/${SAMPLE}.ecoli.sorted.bam"

done

echo "All samples aligned, duplicates marked, MAPQ>=$PARAM_MAPQ BAMs created."
echo

###############################################################################
# STEP 2: Spike-in normalization table
###############################################################################
RAW_TSV="$SPIKE_DIR/spikein_counts_raw.tsv"

echo "sample	mm10_unique	ecoli_unique" > "$RAW_TSV"

for MMLOG in "$LOG_MM10"/*.mm10.bowtie2.log; do
    SAMPLE_BASENAME=$(basename "$MMLOG" .mm10.bowtie2.log)
    ECLOG="$LOG_ECOLI/${SAMPLE_BASENAME}.ecoli.bowtie2.log"

    if [[ ! -f "$ECLOG" ]]; then
        echo "WARNING: No E. coli log for $SAMPLE_BASENAME, skipping."
        continue
    fi

    MM10_UNIQ=$(grep "aligned concordantly exactly 1 time" "$MMLOG" | awk '{print $1}')
    ECOLI_UNIQ=$(grep "aligned concordantly exactly 1 time" "$ECLOG" | awk '{print $1}')

    if [[ -z "$MM10_UNIQ" || -z "$ECOLI_UNIQ" ]]; then
        echo "WARNING: Could not parse unique counts for $SAMPLE_BASENAME, skipping."
        continue
    fi

    echo -e "${SAMPLE_BASENAME}\t${MM10_UNIQ}\t${ECOLI_UNIQ}" >> "$RAW_TSV"
done

echo "Raw spike-in counts written to: $RAW_TSV"
echo "Computing spike fractions and scaleFactors..."

python3 - "$RAW_TSV" "$NORM_TSV" << 'EOF'
import sys

raw_path, out_path = sys.argv[1], sys.argv[2]

rows = []
with open(raw_path) as f:
    header = f.readline()
    for line in f:
        if not line.strip():
            continue
        sample, mm10, ecoli = line.strip().split("\t")
        mm10_u = int(mm10)
        ecoli_u = int(ecoli)
        spike_fraction = float(ecoli_u) / mm10_u if mm10_u > 0 else 0.0
        rows.append((sample, mm10_u, ecoli_u, spike_fraction))

if not rows:
    print("No rows found in spikein_counts_raw.tsv", file=sys.stderr)
    sys.exit(1)

non_zero = [r[3] for r in rows if r[3] > 0]
if not non_zero:
    print("All spike fractions are zero; cannot compute scale factors.", file=sys.stderr)
    sys.exit(1)

min_frac = min(non_zero)

with open(out_path, "w") as out:
    out.write("sample\tmm10_unique\tecoli_unique\tspike_fraction\tscaleFactor\n")
    for sample, mm10_u, ecoli_u, frac in rows:
        if frac > 0:
            sf = min_frac / frac
        else:
            sf = 1.0
        out.write(f"{sample}\t{mm10_u}\t{ecoli_u}\t{frac:.6g}\t{sf:.6g}\n")
EOF

echo "Spike-in normalization table written to: $NORM_TSV"
echo

cd "$PROJECT"
multiqc .

echo "=== 02_align_spikein_v2.sh COMPLETE ==="

} 2>&1 | tee "$PIPELINE_LOG"
