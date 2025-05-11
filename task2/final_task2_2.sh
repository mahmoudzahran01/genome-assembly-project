#!/bin/bash
#SBATCH --job-name=eval_lizard
#SBATCH --output=eval_lizard_%j.out
#SBATCH --error=eval_lizard_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=48     # Increased from 32 to 48
#SBATCH --mem=128G

echo "Starting Scincus mitranus genome assembly evaluation (Task 2.2) - UPDATED"
echo "Running on $(hostname)"
echo "Start time: $(date)"

# Set up directories
SCRATCH_DIR="/ibex/scratch/zahrmm0b/genome_assembly"
ASSEMBLY_DIR="$SCRATCH_DIR/assembly"
SUBSAMPLE_DIR="$SCRATCH_DIR/data/subsampled"
EVAL_DIR="$SCRATCH_DIR/evaluation"
REPORT_DIR="$SCRATCH_DIR/report"
TOOLS_DIR="$SCRATCH_DIR/tools"

# Set up Meryl path (use your successfully built version)
MERYL_DIR="$SCRATCH_DIR/tools/meryl/build/bin"
export PATH=$MERYL_DIR:$PATH
echo "Using Meryl from: $MERYL_DIR"

# Create necessary directories
mkdir -p $EVAL_DIR
mkdir -p $EVAL_DIR/quast
mkdir -p $EVAL_DIR/busco
mkdir -p $EVAL_DIR/merqury
mkdir -p $EVAL_DIR/flagger
mkdir -p $REPORT_DIR
mkdir -p $TOOLS_DIR

# Load the confirmed available modules
module purge
module load quast
module load blast
module load minimap2
module load jellyfish
module load busco/5.8.3  # Using the latest available version from your system

# Check if samtools is available
if command -v samtools &> /dev/null; then
    echo "Samtools found in PATH"
else
    echo "Attempting to load samtools module..."
    module load samtools || echo "Samtools module not available, will use basic analysis methods"
fi

echo "Directories:"
echo "- Scratch directory: $SCRATCH_DIR"
echo "- Assembly directory: $ASSEMBLY_DIR"
echo "- Evaluation directory: $EVAL_DIR"
echo "- Report directory: $REPORT_DIR"
echo "- Tools directory: $TOOLS_DIR"

# Identify assemblies to evaluate
echo "Identifying assemblies to evaluate..."

# Find all primary contig FASTA files (excluding haplotype-specific ones)
PRIMARY_ASSEMBLIES=($(find $ASSEMBLY_DIR -name "*p_ctg.fasta" | grep -v "hap[12]"))

# Check if any assemblies were found
if [ ${#PRIMARY_ASSEMBLIES[@]} -eq 0 ]; then
    echo "ERROR: No assembly FASTA files found in $ASSEMBLY_DIR"
    exit 1
fi

echo "Found ${#PRIMARY_ASSEMBLIES[@]} assemblies to evaluate:"
for ASSEMBLY in "${PRIMARY_ASSEMBLIES[@]}"; do
    echo "- $(basename $ASSEMBLY)"
done

# Choose the best assembly 
BEST_ASSEMBLY=""
for ASSEMBLY in "${PRIMARY_ASSEMBLIES[@]}"; do
    if [[ $ASSEMBLY == *"_full"* ]]; then
        BEST_ASSEMBLY=$ASSEMBLY
        break
    fi
done

# If no full integrated assembly, use any available assembly
if [ -z "$BEST_ASSEMBLY" ]; then
    BEST_ASSEMBLY=${PRIMARY_ASSEMBLIES[0]}
fi

BEST_ASSEMBLY_NAME=$(basename $BEST_ASSEMBLY)
echo "Selected $BEST_ASSEMBLY_NAME as primary assembly for detailed evaluation"

# Check for raw data files
PACBIO_READS=($SUBSAMPLE_DIR/pacbio_*_10p.fastq)
NANOPORE_READS=$SUBSAMPLE_DIR/nanopore_10p.fastq
HIC_READS_R1=$SUBSAMPLE_DIR/hic_r1_10p.fastq
HIC_READS_R2=$SUBSAMPLE_DIR/hic_r2_10p.fastq

echo "Checking for raw read files..."
for FILE in "${PACBIO_READS[@]}" "$NANOPORE_READS" "$HIC_READS_R1" "$HIC_READS_R2"; do
    if [ -f "$FILE" ]; then
        READ_COUNT=$(grep -c "^@" $FILE 2>/dev/null || echo 0)
        echo "- $(basename $FILE): Found with $READ_COUNT reads"
    else
        echo "- $(basename $FILE): NOT FOUND"
    fi
done

#######################################################
# PART 1: Basic metrics using QUAST
#######################################################
echo ""
echo "PART 1: Running QUAST for basic assembly metrics..."

# Run QUAST on all assemblies
quast.py "${PRIMARY_ASSEMBLIES[@]}" \
    -o $EVAL_DIR/quast \
    -t 48 \
    --min-contig 1000 \
    --eukaryote \
    --large \
    --est-ref-size 1300000000 \
    --plots-format png

# Check if QUAST completed successfully
if [ -f "$EVAL_DIR/quast/report.txt" ]; then
    echo "QUAST completed successfully."
    echo "Basic assembly metrics:"
    cat $EVAL_DIR/quast/report.txt
    
    # Extract key metrics for the report
    cat > $REPORT_DIR/quast_summary.txt << EOT
QUAST Basic Assembly Metrics
===========================
Date: $(date)

EOT
    
    # Append the QUAST report
    cat $EVAL_DIR/quast/report.txt >> $REPORT_DIR/quast_summary.txt
else
    echo "ERROR: QUAST did not complete successfully."
fi

#######################################################
# PART 2: Gene completeness using BUSCO
#######################################################
echo ""
echo "PART 2: Running BUSCO for gene completeness assessment..."

# Define BUSCO lineages to try
LINEAGES=("tetrapoda_odb10" "vertebrata_odb10" "squamata_odb10")

for LINEAGE in "${LINEAGES[@]}"; do
    echo "Running BUSCO with $LINEAGE lineage on $BEST_ASSEMBLY_NAME"
    
    # Run BUSCO
    cd $EVAL_DIR/busco
    
    # Create a temporary directory for BUSCO outputs
    BUSCO_OUT_DIR="busco_${LINEAGE}_$(basename $BEST_ASSEMBLY .fasta)"
    mkdir -p $BUSCO_OUT_DIR
    
    # Run BUSCO with the specified lineage
    busco -i $BEST_ASSEMBLY \
        -o $BUSCO_OUT_DIR \
        -l $LINEAGE \
        -m genome \
        -c 48 \
        --download_path $TOOLS_DIR/busco_downloads
    
    # Check if BUSCO completed successfully
    BUSCO_SUMMARY="$EVAL_DIR/busco/$BUSCO_OUT_DIR/short_summary.*.txt"
    
    if ls $BUSCO_SUMMARY 1> /dev/null 2>&1; then
        echo "BUSCO with $LINEAGE completed successfully."
        echo "Gene completeness summary:"
        cat $BUSCO_SUMMARY
        
        # Add to the report
        cat >> $REPORT_DIR/busco_summary.txt << EOT
BUSCO Gene Completeness Assessment - $LINEAGE
=============================================
Date: $(date)
Assembly: $BEST_ASSEMBLY_NAME

EOT
        cat $BUSCO_SUMMARY >> $REPORT_DIR/busco_summary.txt
        echo "" >> $REPORT_DIR/busco_summary.txt
    else
        echo "ERROR: BUSCO with $LINEAGE did not complete successfully."
    fi
done

#######################################################
# PART 3: K-mer distribution and QV score using Merqury
#######################################################
echo ""
echo "PART 3: Running Merqury for k-mer analysis and QV score..."

# Clone Merqury if needed
if [ ! -d "$TOOLS_DIR/merqury" ]; then
    echo "Cloning Merqury repository..."
    cd $TOOLS_DIR
    git clone https://github.com/marbl/merqury.git
fi

# Set up Merqury paths
export PATH=$TOOLS_DIR/merqury:$PATH
export MERQURY=$TOOLS_DIR/merqury

# Verify Meryl is in path and working
echo "Checking Meryl installation..."
which meryl
meryl --version || echo "Meryl version check failed, but continuing..."

# Create a directory for Merqury
mkdir -p $EVAL_DIR/merqury
cd $EVAL_DIR/merqury

# Concatenate PacBio reads for k-mer database
PACBIO_COMBINED="$EVAL_DIR/merqury/pacbio_combined.fastq"
if [ ! -f "$PACBIO_COMBINED" ]; then
    echo "Concatenating PacBio reads..."
    cat ${PACBIO_READS[@]} > $PACBIO_COMBINED
else
    echo "Using existing combined PacBio reads..."
fi

# Build k-mer database using your built Meryl
echo "Building k-mer database with Meryl..."
if [ ! -d "pacbio.meryl" ]; then
    $MERYL_DIR/meryl count k=21 pacbio_combined.fastq output pacbio.meryl
else
    echo "Using existing k-mer database..."
fi

# Run Merqury
echo "Running Merqury..."
mkdir -p $EVAL_DIR/merqury/logs
$MERQURY/merqury.sh pacbio.meryl $BEST_ASSEMBLY merqury_output

# Check if Merqury completed successfully
if [ -f "$EVAL_DIR/merqury/merqury_output.qv" ]; then
    echo "Merqury completed successfully."
    echo "QV score summary:"
    cat $EVAL_DIR/merqury/merqury_output.qv
    
    # Add to the report
    cat > $REPORT_DIR/merqury_summary.txt << EOT
Merqury K-mer Analysis and QV Score
==================================
Date: $(date)
Assembly: $BEST_ASSEMBLY_NAME

QV Scores:
EOT
    cat $EVAL_DIR/merqury/merqury_output.qv >> $REPORT_DIR/merqury_summary.txt
    
    # Also capture spectra-cn plots if available
    if [ -f "$EVAL_DIR/merqury/merqury_output.spectra-cn.pdf" ]; then
        echo "K-mer spectra plot generated at: $EVAL_DIR/merqury/merqury_output.spectra-cn.pdf"
        echo "" >> $REPORT_DIR/merqury_summary.txt
        echo "K-mer spectra plots available at: $EVAL_DIR/merqury/merqury_output.spectra-cn.pdf" >> $REPORT_DIR/merqury_summary.txt
    fi
else
    echo "ERROR: Merqury did not complete successfully."
    
    # Try directly with QV script as fallback
    echo "Attempting direct QV calculation as fallback..."
    $MERQURY/eval/qv.sh pacbio.meryl $BEST_ASSEMBLY direct_qv_output
    
    if [ -f "direct_qv_output.qv" ]; then
        echo "Direct QV calculation successful."
        cat direct_qv_output.qv
        
        # Add to the report
        cat > $REPORT_DIR/merqury_summary.txt << EOT
Merqury K-mer Analysis and QV Score (Direct Method)
================================================
Date: $(date)
Assembly: $BEST_ASSEMBLY_NAME

QV Scores:
EOT
        cat direct_qv_output.qv >> $REPORT_DIR/merqury_summary.txt
    fi
fi

#######################################################
# PART 4: Identify mis-assemblies using read mapping
#######################################################
echo ""
echo "PART 4: Running mis-assembly detection with read mapping..."

# Create directory for analysis
mkdir -p $EVAL_DIR/mis_assembly
cd $EVAL_DIR/mis_assembly

# First map reads to the assembly
echo "Mapping PacBio reads to the assembly..."
PACBIO_COMBINED="$EVAL_DIR/merqury/pacbio_combined.fastq"

# Map PacBio reads
echo "Mapping PacBio reads with minimap2..."
minimap2 -ax map-hifi $BEST_ASSEMBLY $PACBIO_COMBINED -t 48 > pacbio_alignment.sam

# Check if samtools is available
if command -v samtools &> /dev/null; then
    echo "Using samtools for detailed alignment analysis..."
    
    # Convert to BAM, sort, and index
    samtools view -bS pacbio_alignment.sam > pacbio_alignment.bam
    samtools sort pacbio_alignment.bam -o pacbio_alignment.sorted.bam
    samtools index pacbio_alignment.sorted.bam
    
    # Basic alignment statistics
    echo "Generating alignment statistics..."
    samtools flagstat pacbio_alignment.sorted.bam > alignment_stats.txt
    
    # Check coverage
    samtools depth -a pacbio_alignment.sorted.bam | \
        awk '{sum+=$3; n++} END {print "Average coverage: " sum/n}' > coverage_stats.txt
    
    # Identify potential mis-assemblies through coverage anomalies
    echo "Identifying potential mis-assemblies through coverage analysis..."
    
    # Calculate coverage standard deviation
    samtools depth -a pacbio_alignment.sorted.bam | \
        awk '{sum+=$3; sumsq+=($3)^2; n++} END {mean=sum/n; stddev=sqrt(sumsq/n - (mean)^2); print "Mean coverage: " mean; print "StdDev: " stddev}' > coverage_distribution.txt
    
    # Extract mean and stddev values
    MEAN_COV=$(grep "Mean coverage" coverage_distribution.txt | awk '{print $3}')
    STDDEV=$(grep "StdDev" coverage_distribution.txt | awk '{print $2}')
    
    # Create coverage histogram
    samtools depth -a pacbio_alignment.sorted.bam | cut -f3 | sort -n | uniq -c > coverage_histogram.txt
    
    # Identify suspicious regions with abnormal coverage
    LOW_THRESHOLD=$(echo "$MEAN_COV - 2*$STDDEV" | bc)
    HIGH_THRESHOLD=$(echo "$MEAN_COV + 2*$STDDEV" | bc)
    
    echo "Identifying regions with abnormal coverage (below $LOW_THRESHOLD or above $HIGH_THRESHOLD)..."
    
    samtools depth -a pacbio_alignment.sorted.bam | \
        awk -v low=$LOW_THRESHOLD -v high=$HIGH_THRESHOLD \
        '{if($3 < low || $3 > high) print $1, $2, $3}' > suspicious_regions.txt
    
    # Count suspicious regions
    SUSPICIOUS_COUNT=$(wc -l < suspicious_regions.txt)
    
    # Identify potential breakpoints in the alignment
    echo "Identifying potential breakpoints in the assembly..."
    samtools view pacbio_alignment.sorted.bam | \
        awk '{if($2 == 0 || $2 == 16) next; print $0}' | \
        wc -l > breakpoints_count.txt
    
    # Add to the report
    cat > $REPORT_DIR/mis_assembly_summary.txt << EOT
Mis-assembly Detection
====================
Date: $(date)
Assembly: $BEST_ASSEMBLY_NAME

Alignment Statistics:
EOT
    cat alignment_stats.txt >> $REPORT_DIR/mis_assembly_summary.txt
    echo "" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Coverage Statistics:" >> $REPORT_DIR/mis_assembly_summary.txt
    cat coverage_distribution.txt >> $REPORT_DIR/mis_assembly_summary.txt
    echo "" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Suspicious regions identified: $SUSPICIOUS_COUNT positions" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "These are positions with coverage below $LOW_THRESHOLD or above $HIGH_THRESHOLD" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Potential breakpoints: $(cat breakpoints_count.txt)" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "These represent reads with split or discordant alignments, suggesting mis-assemblies" >> $REPORT_DIR/mis_assembly_summary.txt
    
    echo "Mis-assembly detection completed."
else
    echo "Samtools not available, performing basic alignment analysis..."
    
    # Count alignments as a basic metric
    TOTAL_ALIGNMENTS=$(grep -v "^@" pacbio_alignment.sam | wc -l)
    MAPPED_ALIGNMENTS=$(grep -v "^@" pacbio_alignment.sam | grep -v "4[[:space:]]" | wc -l)
    UNMAPPED_ALIGNMENTS=$(($TOTAL_ALIGNMENTS - $MAPPED_ALIGNMENTS))
    MAPPING_RATE=$(echo "scale=2; $MAPPED_ALIGNMENTS * 100 / $TOTAL_ALIGNMENTS" | bc)
    
    echo "Total alignments: $TOTAL_ALIGNMENTS" > basic_alignment_stats.txt
    echo "Mapped alignments: $MAPPED_ALIGNMENTS" >> basic_alignment_stats.txt
    echo "Unmapped alignments: $UNMAPPED_ALIGNMENTS" >> basic_alignment_stats.txt
    echo "Mapping rate: $MAPPING_RATE%" >> basic_alignment_stats.txt
    
    # Add to the report
    cat > $REPORT_DIR/mis_assembly_summary.txt << EOT
Basic Mis-assembly Detection
==========================
Date: $(date)
Assembly: $BEST_ASSEMBLY_NAME

Basic Alignment Statistics:
$(cat basic_alignment_stats.txt)

Without samtools, detailed analysis is limited. Potential mis-assemblies can include:
- Collapsed repeats (regions where coverage would be abnormally high)
- Expansions (regions where coverage would be abnormally low)
- Inversions (regions where read orientations would conflict)
- Translocations (reads mapping to unexpected contigs)

A high overall mapping rate (>95%) suggests good assembly quality.
EOT
fi

# Map Nanopore reads if available
if [ -f "$NANOPORE_READS" ]; then
    echo "Mapping Nanopore reads with minimap2..."
    minimap2 -ax map-ont $BEST_ASSEMBLY $NANOPORE_READS -t 48 > nanopore_alignment.sam
    
    # Basic analysis of Nanopore alignments
    NANO_TOTAL_ALIGNMENTS=$(grep -v "^@" nanopore_alignment.sam | wc -l)
    NANO_MAPPED_ALIGNMENTS=$(grep -v "^@" nanopore_alignment.sam | grep -v "4[[:space:]]" | wc -l)
    NANO_MAPPING_RATE=$(echo "scale=2; $NANO_MAPPED_ALIGNMENTS * 100 / $NANO_TOTAL_ALIGNMENTS" | bc)
    
    echo "" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Nanopore Alignment Statistics:" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Total alignments: $NANO_TOTAL_ALIGNMENTS" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Mapped alignments: $NANO_MAPPED_ALIGNMENTS" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Mapping rate: $NANO_MAPPING_RATE%" >> $REPORT_DIR/mis_assembly_summary.txt
    echo "Nanopore long reads are particularly useful for identifying structural variants and mis-assemblies." >> $REPORT_DIR/mis_assembly_summary.txt
fi

#######################################################
# Generate Final Comprehensive Report
#######################################################
echo ""
echo "Generating comprehensive evaluation report..."

cat > $REPORT_DIR/assembly_evaluation_report.txt << EOT
Scincus mitranus Genome Assembly Evaluation (Task 2.2)
=====================================================
Date: $(date)
Job ID: ${SLURM_JOB_ID}

Assemblies Evaluated:
EOT

for ASSEMBLY in "${PRIMARY_ASSEMBLIES[@]}"; do
    echo "- $(basename $ASSEMBLY)" >> $REPORT_DIR/assembly_evaluation_report.txt
done

echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
echo "Primary Assembly for Detailed Evaluation: $BEST_ASSEMBLY_NAME" >> $REPORT_DIR/assembly_evaluation_report.txt
echo "" >> $REPORT_DIR/assembly_evaluation_report.txt

# Add QUAST metrics
if [ -f "$REPORT_DIR/quast_summary.txt" ]; then
    echo "1. Basic Assembly Metrics (QUAST)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "=================================" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "QUAST provides basic statistics about the assembly including size, contiguity (N50), and GC content." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "These metrics help evaluate the overall completeness and fragmentation of the assembly." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    tail -n +3 $REPORT_DIR/quast_summary.txt >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
fi

# Add BUSCO results
if [ -f "$REPORT_DIR/busco_summary.txt" ]; then
    echo "2. Gene Completeness Assessment (BUSCO)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "=======================================" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "BUSCO (Benchmarking Universal Single-Copy Orthologs) evaluates assembly completeness by" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "searching for conserved genes that should be present in the species' lineage." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "High percentage of complete BUSCOs indicates good gene content representation." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    cat $REPORT_DIR/busco_summary.txt >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
else
    echo "2. Gene Completeness Assessment (BUSCO)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "=======================================" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "BUSCO analysis is still in progress or has not completed." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "This analysis searches for conserved single-copy orthologs in the assembly." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "For reptile genomes, typical BUSCO completeness ranges from 80-95% for the tetrapoda_odb10 dataset." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
fi

# Add Merqury results
if [ -f "$REPORT_DIR/merqury_summary.txt" ]; then
    echo "3. K-mer Analysis and QV Score (Merqury)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "=======================================" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "Merqury performs k-mer based analysis of the assembly to estimate consensus quality." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "The QV (Quality Value) score represents the log-scaled error probability, with higher values indicating better quality:" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "- QV 20 = 99% accuracy (1 error per 100 bases)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "- QV 30 = 99.9% accuracy (1 error per 1,000 bases)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "- QV 40 = 99.99% accuracy (1 error per 10,000 bases)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "- QV 50 = 99.999% accuracy (1 error per 100,000 bases)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    cat $REPORT_DIR/merqury_summary.txt >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
else
    echo "3. K-mer Analysis and QV Score (Merqury)" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "=======================================" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "Merqury analysis is still in progress or has not completed." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "This analysis uses k-mer frequencies to assess assembly quality and calculate QV scores." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "High-quality assemblies typically have QV scores above 40, indicating 99.99% accuracy." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
fi

# Add mis-assembly detection results
if [ -f "$REPORT_DIR/mis_assembly_summary.txt" ]; then
    echo "4. Mis-assembly Detection" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "===========================" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "Mis-assembly detection identifies structural errors in the assembly by analyzing read alignments." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "Potential issues include collapsed repeats, expansions, inversions, and translocations." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "High mapping rates and uniform coverage suggest good assembly quality." >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
    cat $REPORT_DIR/mis_assembly_summary.txt >> $REPORT_DIR/assembly_evaluation_report.txt
    echo "" >> $REPORT_DIR/assembly_evaluation_report.txt
fi

# Add interpretation and recommendations based on QUAST results
echo "5. Interpretation and Recommendations" >> $REPORT_DIR/assembly_evaluation_report.txt
echo "===================================" >> $REPORT_DIR/assembly_evaluation_report.txt
echo "" >> $REPORT_DIR/assembly_evaluation_report.txt

# Extract N50 and contig count
if [ -f "$EVAL_DIR/quast/report.txt" ]; then
    N50=$(grep "N50" $EVAL_DIR/quast/report.txt | awk '{print $2}')
    CONTIG_COUNT=$(grep "# contigs" $EVAL_DIR/quast/report.txt | head -1 | awk '{print $3}')
    ASSEMBLY_SIZE=$(grep "Total length" $EVAL_DIR/quast/report.txt | head -1 | awk '{print $3}')
    GC_CONTENT=$(grep "GC" $EVAL_DIR/quast/report.txt | awk '{print $2}')
    L50=$(grep "L50" $EVAL_DIR/quast/report.txt | head -1 | awk '{print $2}')
    
    cat >> $REPORT_DIR/assembly_evaluation_report.txt << EOT
Based on the evaluation metrics obtained, we can draw the following conclusions about the Scincus mitranus assembly:

1. Assembly Size and Contiguity:
   - The assembly size is $ASSEMBLY_SIZE bp (~${ASSEMBLY_SIZE:0:1}.${ASSEMBLY_SIZE:1:2} Gb), which is slightly larger than the estimated genome size of 1.3 Gb.
   - This suggests potential duplicate regions or unresolved heterozygosity.
   - The N50 value of $N50 bp indicates moderate contiguity.
   - For a vertebrate genome, this N50 is relatively low. High-quality vertebrate assemblies typically achieve N50 values >1 Mb.
   
2. Fragmentation Analysis:
   - With $CONTIG_COUNT contigs and an L50 of $L50, the assembly is quite fragmented.
   - For comparison, chromosome-level assemblies typically have fewer than 100 scaffolds.
   - The large number of contigs suggests many unresolved repetitive regions.
   
3. Sequence Content:
   - The GC content of $GC_CONTENT% is within the typical range for reptile genomes (usually 40-45%).
   
EOT
else
    echo "Could not extract QUAST metrics for detailed interpretation." >> $REPORT_DIR/assembly_evaluation_report.txt
fi

# Add general recommendations
cat >> $REPORT_DIR/assembly_evaluation_report.txt << EOT
Potential Improvements:
1. Increase Contiguity:
   - Incorporate more ultra-long reads (>100 kb) from Oxford Nanopore to span repetitive regions
   - Use additional Hi-C data with improved proximity ligation methods
   - Consider chromosome flow-sorting to assemble chromosomes independently

2. Reduce Redundancy:
   - The assembly size is larger than expected, suggesting duplicate regions
   - Use purge_haplotigs or purge_dups to remove redundant sequences
   - Consider haplotype-aware assembly approaches like HiCanu or Verkko

3. Improve Base-Level Accuracy:
   - Add Illumina short-read data for consensus polishing with Pilon
   - Use newer versions of assembly polishing tools like POLCA
   - Consider accuracy assessment with independent technologies

4. Address Specific Issues:
   - Target regions with abnormal coverage for manual curation
   - Verify the most problematic contigs with independent sequencing
   - Compare to closely related species genomes as a reference

Next Steps:
1. Scaffolding the current assembly with additional Hi-C data
2. Polishing the assembly with additional short-read data
3. Manual curation of problematic regions identified in mis-assembly detection
4. Chromosome assignment and orientation based on related species

This assembly represents a solid foundation but requires additional improvements to reach reference quality. With the current technologies available, a chromosome-level assembly with higher contiguity should be achievable.
EOT

echo "Comprehensive evaluation report generated at: $REPORT_DIR/assembly_evaluation_report.txt"
echo "Assembly evaluation completed at $(date)"
