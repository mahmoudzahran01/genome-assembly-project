#!/bin/bash
#SBATCH --job-name=lizard
#SBATCH --output=lizard_%j.out
#SBATCH --error=lizard_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

echo "Starting Scincus mitranus genome assembly (Task 2.1)"
echo "Running on $(hostname)"
echo "Start time: $(date)"

# Set up directories with exact paths
SCRATCH_DIR="/ibex/scratch/zahrmm0b/genome_assembly"
SUBSAMPLE_DIR="$SCRATCH_DIR/data/subsampled"
ASSEMBLY_DIR="$SCRATCH_DIR/assembly"
REPORT_DIR="$SCRATCH_DIR/report"

# Make sure output directories exist
mkdir -p $ASSEMBLY_DIR
mkdir -p $REPORT_DIR

echo "Working with the following directories:"
echo "- Scratch directory: $SCRATCH_DIR"
echo "- Subsampled data: $SUBSAMPLE_DIR"
echo "- Assembly output: $ASSEMBLY_DIR"
echo "- Report directory: $REPORT_DIR"

# Define the datasets
PACBIO_SEQ="/ibex/reference/course/cs249/lizard/input/ncbi_upload/lizard_liver_seq.fastq.gz"
PACBIO_REV="/ibex/reference/course/cs249/lizard/input/ncbi_upload/lizard_liver_rev.fastq.gz"
NANOPORE="/ibex/reference/course/cs249/lizard/input/ncbi_upload/lizard_ont.fastq.gz"
HIC_R1="/ibex/reference/course/cs249/lizard/input/ncbi_upload/lizard_hic_R1.fastq.gz"
HIC_R2="/ibex/reference/course/cs249/lizard/input/ncbi_upload/lizard_hic_R2.fastq.gz"

# Define paths for subsampled files
SUB_SEQ="$SUBSAMPLE_DIR/pacbio_seq_10p.fastq"
SUB_REV="$SUBSAMPLE_DIR/pacbio_rev_10p.fastq"
SUB_ONT="$SUBSAMPLE_DIR/nanopore_10p.fastq"
SUB_HIC_R1="$SUBSAMPLE_DIR/hic_r1_10p.fastq"
SUB_HIC_R2="$SUBSAMPLE_DIR/hic_r2_10p.fastq"

# Load necessary modules
echo "Loading modules..."
module purge
module load seqtk
module load hifiasm/0.19.5

# Log input files
echo "Input files:"
echo "- PacBio SEQ: $PACBIO_SEQ"
echo "- PacBio REV: $PACBIO_REV"
echo "- Nanopore: $NANOPORE"
echo "- Hi-C R1: $HIC_R1"
echo "- Hi-C R2: $HIC_R2"

echo "Subsampled files location: $SUBSAMPLE_DIR"
echo "Checking subsampled files..."

# Function to check if a file exists and has content
file_has_content() {
    local file=$1
    if [ -f "$file" ] && [ -s "$file" ]; then
        local read_count=$(grep -c "^@" "$file")
        if [ "$read_count" -gt 0 ]; then
            return 0  # File exists and has content
        fi
    fi
    return 1  # File doesn't exist or has no content
}

# Function to subsample a file if needed
subsample_if_needed() {
    local source_file=$1
    local target_file=$2
    local temp_dir=$3
    
    if file_has_content "$target_file"; then
        echo "$(basename $target_file) exists and has reads, skipping subsampling."
        return 0
    else
        echo "Subsampling $(basename $source_file) to create $(basename $target_file)..."
        gunzip -c $source_file > $temp_dir/temp.fastq
        seqtk sample -s42 $temp_dir/temp.fastq 0.1 > $target_file
        rm $temp_dir/temp.fastq
        if file_has_content "$target_file"; then
            echo "Subsampling successful."
            return 0
        else
            echo "ERROR: Subsampling failed for $(basename $source_file)."
            return 1
        fi
    fi
}

# Check if subsampled files exist and have content
for FILE in "$SUB_SEQ" "$SUB_REV" "$SUB_ONT" "$SUB_HIC_R1" "$SUB_HIC_R2"; do
    if file_has_content "$FILE"; then
        READS=$(grep -c "^@" $FILE)
        echo "$(basename $FILE) exists and has $READS reads"
    else
        echo "WARNING: $(basename $FILE) does not exist or is empty"
    fi
done

# Get read counts for the report
READS_SEQ=$(grep -c "^@" $SUB_SEQ 2>/dev/null || echo 0)
READS_REV=$(grep -c "^@" $SUB_REV 2>/dev/null || echo 0)
READS_ONT=$(grep -c "^@" $SUB_ONT 2>/dev/null || echo 0)
READS_HIC_R1=$(grep -c "^@" $SUB_HIC_R1 2>/dev/null || echo 0)
READS_HIC_R2=$(grep -c "^@" $SUB_HIC_R2 2>/dev/null || echo 0)

echo "Current subsampled read counts:"
echo "- PacBio SEQ: $READS_SEQ reads"
echo "- PacBio REV: $READS_REV reads"
echo "- Nanopore: $READS_ONT reads"
echo "- Hi-C R1: $READS_HIC_R1 reads"
echo "- Hi-C R2: $READS_HIC_R2 reads"

# Create a temp directory for processing
TEMP_DIR="$SCRATCH_DIR/temp"
mkdir -p $TEMP_DIR

# Process any missing or empty subsampled files
echo "Processing any missing or empty subsampled files..."

# Process PacBio SEQ
if ! file_has_content "$SUB_SEQ"; then
    subsample_if_needed "$PACBIO_SEQ" "$SUB_SEQ" "$TEMP_DIR"
fi

# Process PacBio REV
if ! file_has_content "$SUB_REV"; then
    subsample_if_needed "$PACBIO_REV" "$SUB_REV" "$TEMP_DIR"
fi

# Process Nanopore
if ! file_has_content "$SUB_ONT"; then
    subsample_if_needed "$NANOPORE" "$SUB_ONT" "$TEMP_DIR"
fi

# Process Hi-C R1
if ! file_has_content "$SUB_HIC_R1"; then
    subsample_if_needed "$HIC_R1" "$SUB_HIC_R1" "$TEMP_DIR"
fi

# Process Hi-C R2
if ! file_has_content "$SUB_HIC_R2"; then
    subsample_if_needed "$HIC_R2" "$SUB_HIC_R2" "$TEMP_DIR"
fi

# Get updated read counts
READS_SEQ=$(grep -c "^@" $SUB_SEQ 2>/dev/null || echo 0)
READS_REV=$(grep -c "^@" $SUB_REV 2>/dev/null || echo 0)
READS_ONT=$(grep -c "^@" $SUB_ONT 2>/dev/null || echo 0)
READS_HIC_R1=$(grep -c "^@" $SUB_HIC_R1 2>/dev/null || echo 0)
READS_HIC_R2=$(grep -c "^@" $SUB_HIC_R2 2>/dev/null || echo 0)

echo "Updated subsampled read counts:"
echo "- PacBio SEQ: $READS_SEQ reads"
echo "- PacBio REV: $READS_REV reads"
echo "- Nanopore: $READS_ONT reads"
echo "- Hi-C R1: $READS_HIC_R1 reads"
echo "- Hi-C R2: $READS_HIC_R2 reads"

# Check for existing assembly files
echo "Checking for existing assembly files..."
cd $ASSEMBLY_DIR

# Find the most recent primary contig GFA file (excluding haplotype specific files)
P_CTG_GFA=$(find $ASSEMBLY_DIR -name "scincus_mitranus*.p_ctg.gfa" | grep -v "hap[12]" | sort -r | head -1)

# PHASE 1: PacBio Assembly (if needed)
if [ -n "$P_CTG_GFA" ]; then
    echo "Found existing assembly file: $P_CTG_GFA"
    # Extract the base name for output files
    BASE_NAME=$(basename $P_CTG_GFA .p_ctg.gfa)
    echo "Using base name: $BASE_NAME"
    
    # Check if the corresponding FASTA file exists
    FASTA_OUTPUT="$ASSEMBLY_DIR/${BASE_NAME}.p_ctg.fasta"
    
    if [ ! -f "$FASTA_OUTPUT" ]; then
        echo "FASTA file doesn't exist. Converting GFA to FASTA..."
        awk '/^S/{print ">"$2"\n"$3}' $P_CTG_GFA > $FASTA_OUTPUT
    else
        echo "FASTA file already exists: $FASTA_OUTPUT"
    fi
else
    echo "No existing assembly found. Need to run Hifiasm."
    
    # Check if PacBio files have content before proceeding
    if [ "$READS_SEQ" -eq 0 ] || [ "$READS_REV" -eq 0 ]; then
        echo "ERROR: PacBio files have no reads. Cannot proceed with assembly."
        echo "Please check that the files exist at: $SUBSAMPLE_DIR"
        exit 1
    fi
    
    # Run with PacBio data
    echo "Running HiFi assembly with PacBio data"
    hifiasm -o scincus_mitranus -t 32 $SUB_SEQ $SUB_REV
    
    # Check for any p_ctg.gfa file pattern
    P_CTG_GFA=$(find $ASSEMBLY_DIR -name "scincus_mitranus*.p_ctg.gfa" | grep -v "hap[12]" | sort -r | head -1)
    
    if [ -z "$P_CTG_GFA" ]; then
        echo "ERROR: Assembly did not produce any primary contig GFA file."
        exit 1
    fi
    
    # Extract the base name for output files
    BASE_NAME=$(basename $P_CTG_GFA .p_ctg.gfa)
    echo "Assembly successful. Using base name: $BASE_NAME"
    
    # Convert GFA to FASTA
    echo "Converting GFA to FASTA..."
    FASTA_OUTPUT="$ASSEMBLY_DIR/${BASE_NAME}.p_ctg.fasta"
    awk '/^S/{print ">"$2"\n"$3}' $P_CTG_GFA > $FASTA_OUTPUT
fi

# PHASE 2: Process the Nanopore data
echo "Processing Nanopore data..."
if [ "$READS_ONT" -eq 0 ]; then
    echo "WARNING: Nanopore data has no reads. Skipping Nanopore processing."
else
    echo "Nanopore data available with $READS_ONT reads."
    
    # Run Hifiasm with PacBio and Nanopore data
    echo "Running assembly with PacBio + Nanopore data"
    NANOPORE_OUTPUT="scincus_mitranus_nanopore"
    
    hifiasm -o $NANOPORE_OUTPUT -t 32 --ul $SUB_ONT $SUB_SEQ $SUB_REV
    
    # Check if the Nanopore-assisted assembly was successful
    NANO_P_CTG_GFA=$(find $ASSEMBLY_DIR -name "${NANOPORE_OUTPUT}*.p_ctg.gfa" | grep -v "hap[12]" | sort -r | head -1)
    
    if [ -n "$NANO_P_CTG_GFA" ]; then
        echo "Nanopore-assisted assembly successful: $NANO_P_CTG_GFA"
        # Convert to FASTA
        NANO_BASE_NAME=$(basename $NANO_P_CTG_GFA .p_ctg.gfa)
        NANO_FASTA_OUTPUT="$ASSEMBLY_DIR/${NANO_BASE_NAME}.p_ctg.fasta"
        awk '/^S/{print ">"$2"\n"$3}' $NANO_P_CTG_GFA > $NANO_FASTA_OUTPUT
        
        # Calculate statistics
        NANO_CONTIG_COUNT=$(grep -c "^>" $NANO_FASTA_OUTPUT)
        NANO_ASSEMBLY_SIZE=$(grep -v "^>" $NANO_FASTA_OUTPUT | wc -c)
        
        echo "Nanopore-assisted assembly statistics:"
        echo "- Number of contigs: $NANO_CONTIG_COUNT"
        echo "- Total assembly size: $NANO_ASSEMBLY_SIZE bp"
    else
        echo "WARNING: Nanopore-assisted assembly did not produce a primary contig file."
    fi
fi

# PHASE 3: Process the Hi-C data
echo "Processing Hi-C data..."
if [ "$READS_HIC_R1" -eq 0 ] || [ "$READS_HIC_R2" -eq 0 ]; then
    echo "WARNING: Hi-C data has missing or empty files. Skipping Hi-C processing."
else
    echo "Hi-C data available with $READS_HIC_R1 and $READS_HIC_R2 reads."
    
    # Run Hifiasm with PacBio and Hi-C data
    echo "Running assembly with PacBio + Hi-C data"
    HIC_OUTPUT="scincus_mitranus_hic"
    
    hifiasm -o $HIC_OUTPUT -t 32 --h1 $SUB_HIC_R1 --h2 $SUB_HIC_R2 $SUB_SEQ $SUB_REV
    
    # Check if the Hi-C-assisted assembly was successful
    HIC_P_CTG_GFA=$(find $ASSEMBLY_DIR -name "${HIC_OUTPUT}*.p_ctg.gfa" | grep -v "hap[12]" | sort -r | head -1)
    
    if [ -n "$HIC_P_CTG_GFA" ]; then
        echo "Hi-C-assisted assembly successful: $HIC_P_CTG_GFA"
        # Convert to FASTA
        HIC_BASE_NAME=$(basename $HIC_P_CTG_GFA .p_ctg.gfa)
        HIC_FASTA_OUTPUT="$ASSEMBLY_DIR/${HIC_BASE_NAME}.p_ctg.fasta"
        awk '/^S/{print ">"$2"\n"$3}' $HIC_P_CTG_GFA > $HIC_FASTA_OUTPUT
        
        # Calculate statistics
        HIC_CONTIG_COUNT=$(grep -c "^>" $HIC_FASTA_OUTPUT)
        HIC_ASSEMBLY_SIZE=$(grep -v "^>" $HIC_FASTA_OUTPUT | wc -c)
        
        echo "Hi-C-assisted assembly statistics:"
        echo "- Number of contigs: $HIC_CONTIG_COUNT"
        echo "- Total assembly size: $HIC_ASSEMBLY_SIZE bp"
    else
        echo "WARNING: Hi-C-assisted assembly did not produce a primary contig file."
    fi
fi

# PHASE 4: Full integrated assembly with all datasets
echo "Processing full integrated assembly with all datasets..."
if [ "$READS_SEQ" -gt 0 ] && [ "$READS_REV" -gt 0 ] && [ "$READS_ONT" -gt 0 ] && [ "$READS_HIC_R1" -gt 0 ] && [ "$READS_HIC_R2" -gt 0 ]; then
    echo "All datasets available. Running full integrated assembly."
    
    FULL_OUTPUT="scincus_mitranus_full"
    
    hifiasm -o $FULL_OUTPUT -t 32 --ul $SUB_ONT --h1 $SUB_HIC_R1 --h2 $SUB_HIC_R2 $SUB_SEQ $SUB_REV
    
    # Check if the full integrated assembly was successful
    FULL_P_CTG_GFA=$(find $ASSEMBLY_DIR -name "${FULL_OUTPUT}*.p_ctg.gfa" | grep -v "hap[12]" | sort -r | head -1)
    
    if [ -n "$FULL_P_CTG_GFA" ]; then
        echo "Full integrated assembly successful: $FULL_P_CTG_GFA"
        # Convert to FASTA
        FULL_BASE_NAME=$(basename $FULL_P_CTG_GFA .p_ctg.gfa)
        FULL_FASTA_OUTPUT="$ASSEMBLY_DIR/${FULL_BASE_NAME}.p_ctg.fasta"
        awk '/^S/{print ">"$2"\n"$3}' $FULL_P_CTG_GFA > $FULL_FASTA_OUTPUT
        
        # Calculate statistics
        FULL_CONTIG_COUNT=$(grep -c "^>" $FULL_FASTA_OUTPUT)
        FULL_ASSEMBLY_SIZE=$(grep -v "^>" $FULL_FASTA_OUTPUT | wc -c)
        
        echo "Full integrated assembly statistics:"
        echo "- Number of contigs: $FULL_CONTIG_COUNT"
        echo "- Total assembly size: $FULL_ASSEMBLY_SIZE bp"
    else
        echo "WARNING: Full integrated assembly did not produce a primary contig file."
    fi
else
    echo "Not all datasets are available. Skipping full integrated assembly."
fi

# Calculate assembly statistics for the original assembly
CONTIG_COUNT=$(grep -c "^>" $FASTA_OUTPUT)
ASSEMBLY_SIZE=$(grep -v "^>" $FASTA_OUTPUT | wc -c)

echo "Original assembly statistics:"
echo "- Number of contigs: $CONTIG_COUNT"
echo "- Total assembly size: $ASSEMBLY_SIZE bp"

# Create a report with assembly information
cat > $REPORT_DIR/assembly_report.txt << EOT
Scincus mitranus Genome Assembly (Task 2.1)
===========================================
Date: $(date)
Job ID: ${SLURM_JOB_ID}
Assembly method: Hifiasm
Subsampling ratio: 10%

Files used:
- PacBio SEQ: $PACBIO_SEQ (subsampled to 10%, $READS_SEQ reads)
- PacBio REV: $PACBIO_REV (subsampled to 10%, $READS_REV reads)
- Nanopore: $NANOPORE (subsampled to 10%, $READS_ONT reads)
- Hi-C R1: $HIC_R1 (subsampled to 10%, $READS_HIC_R1 reads)
- Hi-C R2: $HIC_R2 (subsampled to 10%, $READS_HIC_R2 reads)

Original Assembly Statistics:
- File: $P_CTG_GFA
- FASTA: $FASTA_OUTPUT
- Number of contigs: $CONTIG_COUNT
- Total assembly size: $ASSEMBLY_SIZE bp
EOT

# Add Nanopore assembly stats to report if available
if [ -n "${NANO_FASTA_OUTPUT:-}" ]; then
    cat >> $REPORT_DIR/assembly_report.txt << EOT

Nanopore-assisted Assembly:
- File: $NANO_P_CTG_GFA
- FASTA: $NANO_FASTA_OUTPUT
- Number of contigs: $NANO_CONTIG_COUNT
- Total assembly size: $NANO_ASSEMBLY_SIZE bp
EOT
fi

# Add Hi-C assembly stats to report if available
if [ -n "${HIC_FASTA_OUTPUT:-}" ]; then
    cat >> $REPORT_DIR/assembly_report.txt << EOT

Hi-C-assisted Assembly:
- File: $HIC_P_CTG_GFA
- FASTA: $HIC_FASTA_OUTPUT
- Number of contigs: $HIC_CONTIG_COUNT
- Total assembly size: $HIC_ASSEMBLY_SIZE bp
EOT
fi

# Add full integrated assembly stats to report if available
if [ -n "${FULL_FASTA_OUTPUT:-}" ]; then
    cat >> $REPORT_DIR/assembly_report.txt << EOT

Full Integrated Assembly:
- File: $FULL_P_CTG_GFA
- FASTA: $FULL_FASTA_OUTPUT
- Number of contigs: $FULL_CONTIG_COUNT
- Total assembly size: $FULL_ASSEMBLY_SIZE bp
EOT
fi

echo "Assembly report saved to: $REPORT_DIR/assembly_report.txt"
echo "All results kept in scratch space."
echo "Original assembly: $FASTA_OUTPUT"

# Display all generated assemblies
echo "All generated assemblies:"
if [ -n "${FASTA_OUTPUT:-}" ]; then echo "- Original PacBio: $FASTA_OUTPUT"; fi
if [ -n "${NANO_FASTA_OUTPUT:-}" ]; then echo "- PacBio + Nanopore: $NANO_FASTA_OUTPUT"; fi
if [ -n "${HIC_FASTA_OUTPUT:-}" ]; then echo "- PacBio + Hi-C: $HIC_FASTA_OUTPUT"; fi
if [ -n "${FULL_FASTA_OUTPUT:-}" ]; then echo "- Full integrated: $FULL_FASTA_OUTPUT"; fi

# Clean up
rm -rf $TEMP_DIR
