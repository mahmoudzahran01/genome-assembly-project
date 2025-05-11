#!/bin/bash
#SBATCH --job-name=improve_assembly_subsampled
#SBATCH --output=improve_assembly_%j.out
#SBATCH --error=improve_assembly_%j.err
#SBATCH --time=24:00:00  # Reduced since we're using subsamples
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G  # Reduced for subsamples

echo "==================================================================================="
echo "Scincus mitranus Genome Assembly Improvement (Task 2.3) - SUBSAMPLED DATA ONLY"
echo "Running on $(hostname)"
echo "Start time: $(date)"
echo "==================================================================================="

# Set up directories (using existing subsamples)
SCRATCH_DIR="/ibex/scratch/$USER/genome_assembly"
SUBSAMPLE_DIR="$SCRATCH_DIR/data/subsampled"  # Existing subsamples
ASSEMBLY_DIR="$SCRATCH_DIR/assembly"
QC_DIR="$SCRATCH_DIR/qc"
IMPROVED_DIR="$SCRATCH_DIR/improved"
MERQURY_DIR="$SCRATCH_DIR/merqury"
REPORT_DIR="$SCRATCH_DIR/report"
TOOLS_DIR="$SCRATCH_DIR/tools"
TEMP_DIR="$SCRATCH_DIR/temp"

# Create necessary directories
mkdir -p $QC_DIR/{fastqc,seqstats,kraken,trimmed,corrected}
mkdir -p $IMPROVED_DIR $MERQURY_DIR/{original,improved} $REPORT_DIR $TOOLS_DIR $TEMP_DIR

# Load modules
module purge
module load fastqc trimmomatic porechop kraken2 java samtools seqtk jellyfish merqury/1.3 blast minimap2 hifiasm/0.19.5

# Define subsampled datasets (from your listed files)
HIC_R1="$SUBSAMPLE_DIR/hic_r1_10p.fastq"
HIC_R2="$SUBSAMPLE_DIR/hic_r2_10p.fastq"
NANOPORE="$SUBSAMPLE_DIR/nanopore_10p.fastq"
PACBIO_SEQ="$SUBSAMPLE_DIR/pacbio_seq_10p.fastq"
PACBIO_REV="$SUBSAMPLE_DIR/pacbio_rev_10p.fastq"

# Merge PacBio reads (if both SEQ and REV exist)
if [ -f "$PACBIO_SEQ" ] && [ -f "$PACBIO_REV" ]; then
  cat $PACBIO_SEQ $PACBIO_REV > $TEMP_DIR/pacbio_merged.fastq
  PACBIO_FOR_ASSEMBLY="$TEMP_DIR/pacbio_merged.fastq"
elif [ -f "$PACBIO_SEQ" ]; then
  PACBIO_FOR_ASSEMBLY="$PACBIO_SEQ"
else
  echo "ERROR: No PacBio subsamples found!" >&2
  exit 1
fi

#######################################################
### PART 1: Quality Assessment (Subsampled Data)
#######################################################
echo "Quality assessment of subsampled reads..."
run_fastqc() {
  fastqc -t 8 -o $QC_DIR/fastqc $1
}
get_long_read_stats() {
  seqtk fqchk $1 | head -n 3 > $QC_DIR/seqstats/$(basename $1 .fastq)_stats.txt
}

# Hi-C
run_fastqc $HIC_R1
run_fastqc $HIC_R2

# Nanopore
get_long_read_stats $NANOPORE

# PacBio
get_long_read_stats $PACBIO_FOR_ASSEMBLY

#######################################################
### PART 2: Contamination Screening
#######################################################
echo "Contamination screening with Kraken2..."
run_kraken() {
  kraken2 --db /ibex/reference/KSL/kraken2/standard \
          --threads 16 \
          --output $1.out \
          --report $1.report \
          $2
}

run_kraken "$QC_DIR/kraken/hic" "$HIC_R1"
run_kraken "$QC_DIR/kraken/ont" "$NANOPORE"
run_kraken "$QC_DIR/kraken/pacbio" "$PACBIO_FOR_ASSEMBLY"

#######################################################
### PART 3: Read Filtering
#######################################################
echo "Filtering contaminated reads..."
filter_contamination() {
  seqtk subseq $1 <(grep -E "unclassified|Reptilia" $2.out | cut -f2) > $3
}

# Hi-C (trim first)
trimmomatic PE -threads 16 $HIC_R1 $HIC_R2 \
  $QC_DIR/trimmed/hic_r1_paired.fastq $QC_DIR/trimmed/hic_r1_unpaired.fastq \
  $QC_DIR/trimmed/hic_r2_paired.fastq $QC_DIR/trimmed/hic_r2_unpaired.fastq \
  ILLUMINACLIP:/ibex/reference/KSL/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

filter_contamination "$QC_DIR/trimmed/hic_r1_paired.fastq" "$QC_DIR/kraken/hic" "$QC_DIR/trimmed/hic_r1_clean.fastq"
filter_contamination "$QC_DIR/trimmed/hic_r2_paired.fastq" "$QC_DIR/kraken/hic" "$QC_DIR/trimmed/hic_r2_clean.fastq"

# Nanopore
porechop -i $NANOPORE -o $QC_DIR/trimmed/ont_trimmed.fastq -t 16
filter_contamination "$QC_DIR/trimmed/ont_trimmed.fastq" "$QC_DIR/kraken/ont" "$QC_DIR/trimmed/ont_clean.fastq"

# PacBio (no trimming, just filter)
filter_contamination "$PACBIO_FOR_ASSEMBLY" "$QC_DIR/kraken/pacbio" "$QC_DIR/trimmed/pacbio_clean.fastq"

#######################################################
### PART 4: Error Correction (ONT only)
#######################################################
echo "Error correction for Nanopore reads..."
run_herro() {
  $TOOLS_DIR/herro/build/bin/herro -t 16 -r $1 -i $2 -o $3
}
run_dechat() {
  minimap2 -t 16 -cx map-ont $1 $2 | awk '{if($12>0.8) print}' > $3/dechat_aln.paf
  seqtk subseq $2 <(cut -f1 $3/dechat_aln.paf | sort -u) > $3/dechat_corrected.fastq
}

mkdir -p $QC_DIR/corrected/{herro,dechat}

# Method 1: Herro
if [ ! -d "$TOOLS_DIR/herro" ]; then
  git clone https://github.com/lbcb-sci/herro.git $TOOLS_DIR/herro
  cd $TOOLS_DIR/herro && mkdir build && cd build
  cmake -DCMAKE_BUILD_TYPE=Release .. && make
fi
run_herro "$QC_DIR/trimmed/pacbio_clean.fastq" "$QC_DIR/trimmed/ont_clean.fastq" "$QC_DIR/corrected/herro"

# Method 2: DeChat (simplified)
run_dechat "$QC_DIR/trimmed/pacbio_clean.fastq" "$QC_DIR/trimmed/ont_clean.fastq" "$QC_DIR/corrected/dechat"

# Combine corrected reads (take union of both methods)
cat $QC_DIR/corrected/herro/corrected_reads.fasta $QC_DIR/corrected/dechat/dechat_corrected.fastq | \
  seqtk seq -F '#' > $QC_DIR/corrected/ont_combined.fastq

#######################################################
### PART 5: Reassembly with Hifiasm
#######################################################
echo "Reassembling with improved subsampled reads..."
hifiasm -o $IMPROVED_DIR/scincus_mitranus_subsampled \
        -t 32 \
        --ul $QC_DIR/corrected/ont_combined.fastq \
        --h1 $QC_DIR/trimmed/hic_r1_clean.fastq \
        --h2 $QC_DIR/trimmed/hic_r2_clean.fastq \
        $QC_DIR/trimmed/pacbio_clean.fastq

# Convert GFA to FASTA
awk '/^S/{print ">"$2"\n"$3}' $IMPROVED_DIR/*.p_ctg.gfa > $IMPROVED_DIR/improved_assembly.fasta

#######################################################
### PART 6: QV Evaluation
#######################################################
echo "Evaluating QV scores with Merqury..."
export MERQURY=$TOOLS_DIR/merqury
[ ! -d "$MERQURY" ] && git clone https://github.com/marbl/merqury.git $MERQURY

# Build k-mer DB
jellyfish count -C -m 21 -s 1G -t 16 $QC_DIR/trimmed/pacbio_clean.fastq -o $MERQURY_DIR/reads.jf
$MERQURY/meryl count k=21 $QC_DIR/trimmed/pacbio_clean.fastq output $MERQURY_DIR/reads.meryl

# Evaluate assemblies
$MERQURY/merqury.sh $MERQURY_DIR/reads.meryl $ASSEMBLY_DIR/*.fasta $MERQURY_DIR/original/original
$MERQURY/merqury.sh $MERQURY_DIR/reads.meryl $IMPROVED_DIR/improved_assembly.fasta $MERQURY_DIR/improved/improved

#######################################################
### PART 7: Final Report
#######################################################
echo "Generating final report..."
cat > $REPORT_DIR/improvement_report.txt << EOF
Scincus mitranus Assembly Improvement (Subsampled Data)
=====================================================
Date: $(date)
Job ID: $SLURM_JOB_ID

=== INPUT DATA ===
- Hi-C subsamples: $(basename $HIC_R1), $(basename $HIC_R2)
- Nanopore subsample: $(basename $NANOPORE)
- PacBio subsamples: $(basename $PACBIO_SEQ), $(basename $PACBIO_REV)

=== KEY METRICS ===
$(paste <(echo -e "Metric\nContigs\nAssembly Size (bp)\nN50\nQV Score") \
      <(echo -e "Original\n$(grep -c "^>" $ASSEMBLY_DIR/*.fasta)\n$(grep -v "^>" $ASSEMBLY_DIR/*.fasta | wc -c)\n$(calculate_n50 $ASSEMBLY_DIR/*.fasta)\n$(cat $MERQURY_DIR/original/original.qv 2>/dev/null | awk '{print $2}')") \
      <(echo -e "Improved\n$(grep -c "^>" $IMPROVED_DIR/improved_assembly.fasta)\n$(grep -v "^>" $IMPROVED_DIR/improved_assembly.fasta | wc -c)\n$(calculate_n50 $IMPROVED_DIR/improved_assembly.fasta)\n$(cat $MERQURY_DIR/improved/improved.qv 2>/dev/null | awk '{print $2}')") | column -t -s $'\t')

=== ERROR CORRECTION COMPARISON ===
- Herro corrected: $(grep -c "^@" $QC_DIR/corrected/herro/corrected_reads.fasta)
- DeChat corrected: $(grep -c "^@" $QC_DIR/corrected/dechat/dechat_corrected.fastq)
- Combined ONT reads: $(grep -c "^@" $QC_DIR/corrected/ont_combined.fastq)

=== RECOMMENDATIONS ===
1. For full data: Increase --mem to 192G and --time to 48:00:00
2. To further improve QV: Try different k-mer sizes in Merqury
3. For annotation: Use RNA-Seq data from /ibex/reference/course/cs249/lizard/annotation/
EOF

echo "==================================================================================="
echo "Pipeline completed at $(date)"
echo "Final report: $REPORT_DIR/improvement_report.txt"
echo "==================================================================================="
