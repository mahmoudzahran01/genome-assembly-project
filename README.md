I'll update the README to include Task 2.4 for genome annotation. Here's the complete updated README with the new section:

# Genome Assembly Project

This repository contains implementations of genome assembly algorithms and comprehensive analysis pipelines for both simulated datasets and the *Scincus mitranus* (sandfish lizard) genome.

## Project Overview

This project addresses two main components of genome assembly:

1. **Task 1:** Implementation and evaluation of fundamental genome assembly algorithms
   - De Bruijn Graph (DBG) assembly
   - Overlap-Layout-Consensus (OLC) assembly
   - Application to synthetic datasets and MERS-CoV genome

2. **Task 2:** Assembly of the *Scincus mitranus* genome
   - Assembly using Hifiasm with multiple sequencing technologies
   - Comprehensive evaluation of assembly quality
   - Assembly improvement through read preprocessing and error correction
   - Genome annotation using RNA-Seq data and computational prediction

## Repository Structure

```
.
├── README.md                   # Project documentation
├── dbg_assembler.py            # De Bruijn Graph assembler implementation
├── olc_assembler.py            # Overlap-Layout-Consensus assembler implementation
├── make_gfa.py                 # GFA format generator for visualizations
├── run_task1.py                # Main script for Task 1 execution
├── task2/                      # Task 2 scripts
│   ├── process_all_datasets.sh # Task 2.1: Data processing and assembly
│   ├── final_task2_2.sh        # Task 2.2: Assembly evaluation
│   ├── task2_3_assembly_improvement_subsampled.sh # Task 2.3: Assembly improvement
│   └── lizard_annotation_pipeline.sh # Task 2.4: Genome annotation
└── [additional directories for data, results, etc.]
```

## Requirements

### Core Requirements
- Python 3.8+
- BioPython
- NetworkX
- NumPy
- Matplotlib

### External Tools
- **Assembly Tools:** 
  - Hifiasm (v0.19.5+)
  - SPAdes
- **Quality Assessment:**
  - QUAST
  - BUSCO (v5.8.3+)
  - Merqury (and Meryl)
- **Read Processing:**
  - FastQC
  - Trimmomatic 
  - Porechop
  - Kraken2
  - seqtk
- **Visualization:**
  - Bandage
- **Annotation Tools:**
  - BRAKER
  - CD-HIT
  - InterProScan
  - BLAST
  - seqkit

## Task 1: Fundamental Assembly Algorithms

### 1.1 De Bruijn Graph (DBG) Assembly
Our DBG implementation (`dbg_assembler.py`) constructs a genome assembly by:
- Breaking reads into k-mers of user-defined length
- Building a de Bruijn graph where nodes are (k-1)-mers and edges represent k-mers
- Finding Eulerian paths to identify contigs
- Outputting contigs in FASTA format

### 1.2 Overlap-Layout-Consensus (OLC) Assembly
Our OLC implementation (`olc_assembler.py`) assembles genomes by:
- Computing all-vs-all read overlaps using user-defined minimum overlap length
- Constructing an overlap graph where nodes are reads and edges represent overlaps
- Identifying non-branching paths to form layout
- Calculating consensus sequences from aligned reads
- Outputting contigs in FASTA format

### Running Task 1

The `run_task1.py` script orchestrates the entire Task 1 workflow:

```bash
python run_task1.py --output_dir <path/to/output> [options]
```

**Key Parameters:**
- `--output_dir`: Directory for storing all results (required)
- `--data_dir`: Directory containing input data (default: "task1/data")
- `--skip_dbg`: Skip DBG assembly
- `--skip_olc`: Skip OLC assembly
- `--skip_spades`: Skip SPAdes assembly
- `--skip_evaluation`: Skip evaluation with QUAST

**What the script does:**
1. Creates necessary directory structure
2. Executes DBG assembly with different k-mer sizes on all datasets
   - Constructs assembly graph for reads_b.fastq with k=40 (Task 1.3.1)
   - Runs DBG on reads_r.fastq with k=35 and k=45 (Task 1.3.2)
   - Applies DBG to MERS datasets with appropriate k-mer settings (Task 1.3.3)
3. Executes OLC assembly on all datasets with suitable overlap parameters
4. Runs professional SPAdes assembler for comparison (Task 1.3.4)
5. Evaluates all assemblies with QUAST against reference genomes

**Example usage for specific cases:**
```bash
# Run DBG assembly with k=40 on a specific dataset
python dbg_assembler.py --fastq data/task1/synthetic_1/reads_b.fastq --k 40 --output results/dbg_test

# Run OLC assembly with minimum overlap of 30bp
python olc_assembler.py --fastq data/task1/synthetic_2/no_error_ont_hq_50x.fastq --min-overlap 30 --output results/olc_test
```

## Task 2: Scincus mitranus Genome Assembly

### Task 2.1: Genome Assembly (process_all_datasets.sh)

This script performs the actual genome assembly of *Scincus mitranus* using Hifiasm:

```bash
sbatch task2/process_all_datasets.sh
```

**What the script does:**
1. Sets up directory structure for assembly outputs
2. Subsamples the original datasets to 10% for computational efficiency:
   - PacBio Sequel II reads (SEQ)
   - PacBio Revio reads (REV)
   - Oxford Nanopore reads
   - Hi-C paired-end reads (R1 and R2)
3. Implements a progressive, four-phase assembly approach:
   - **Phase 1:** PacBio-only assembly (base assembly)
   - **Phase 2:** PacBio + Nanopore integration using `--ul` flag
   - **Phase 3:** PacBio + Hi-C integration using `--h1`/`--h2` flags
   - **Phase 4:** Full integrated assembly with all data types
4. Converts the GFA outputs to FASTA format
5. Generates a comprehensive assembly report with statistics

**Resource requirements:** 
- 32 CPUs, 128GB memory, 48-hour runtime
- Runs on Ibex cluster with appropriate modules

### Task 2.2: Assembly Evaluation (final_task2_2.sh)

This script performs a comprehensive evaluation of the *Scincus mitranus* genome assembly:

```bash
sbatch task2/final_task2_2.sh
```

**What the script does:**
1. Identifies all available assemblies from Task 2.1 and selects the best one for detailed evaluation
2. Implements four required evaluation approaches:
   - **Basic Metrics (QUAST):** Assembly size, contiguity (N50), GC content
   - **Gene Completeness (BUSCO):** Tests with tetrapoda, vertebrata, and squamata lineages
   - **K-mer Analysis (Merqury):** K-mer distribution and Quality Value (QV) scores
   - **Mis-assembly Detection:** Read mapping with statistical analysis of coverage
3. Generates a unified report with detailed interpretation of all metrics
4. Provides specific recommendations for assembly improvement

**Resource requirements:**
- 48 CPUs, 128GB memory, 24-hour runtime
- Requires QUAST, BUSCO, Merqury, minimap2, and samtools

### Task 2.3: Assembly Improvement (task2_3_assembly_improvement_subsampled.sh)

This optional script implements quality control and improvements for the assembly:

```bash
sbatch task2/task2_3_assembly_improvement_subsampled.sh
```

**What the script does:**
1. **Quality Assessment:** Evaluates all read datasets with FastQC and seqtk
2. **Contamination Screening:** Uses Kraken2 to identify and filter non-reptilian sequences
3. **Read Filtering and Preprocessing:**
   - Trims adapters from Hi-C reads with Trimmomatic
   - Removes ONT adapters with Porechop
   - Filters contaminated reads from all datasets
4. **Error Correction for ONT Reads:**
   - Implements two methods (Herro and DeChat)
   - Combines results from both approaches
5. **Reassembly with Improved Reads:** Uses Hifiasm with all cleaned datasets
6. **QV Evaluation:** Compares original and improved assemblies
7. **Report Generation:** Creates detailed comparison of assembly quality before and after improvement

**Resource requirements:**
- 32 CPUs, 64GB memory, 24-hour runtime
- Requires multiple bioinformatics tools (FastQC, Kraken2, Trimmomatic, etc.)

### Task 2.4: Genome Annotation (lizard_annotation_pipeline.sh)

This script performs comprehensive genome annotation of the *Scincus mitranus* assembly using RNA-Seq data:

```bash
sbatch task2/lizard_annotation_pipeline.sh
```

**What the script does:**
1. **IsoSeq Data Processing and Quality Control:**
   - Analyzes RNA-Seq data from eye and liver tissues using seqkit
   - Detects potential contamination with Kraken2
   - Checks for transcript completeness (poly-A tails)
   - Clusters isoforms using CD-HIT
   - Aligns transcripts to the genome with minimap2
   - Identifies potential chimeric alignments

2. **Gene Prediction with BRAKER:**
   - Merges RNA-Seq alignments from both tissues
   - Performs soft-masking of the genome for repeat regions
   - Applies BRAKER pipeline for gene prediction using IsoSeq evidence
   - Processes and extracts predicted gene models and protein sequences

3. **Functional Annotation:**
   - Performs BLASTP against UniProt SwissProt database for homology-based annotation
   - Applies InterProScan for domain and functional motif prediction
   - Assigns Gene Ontology (GO) terms based on InterProScan results

4. **Report Generation:**
   - Creates comprehensive reports for each annotation step
   - Generates a final annotation report with integrated results
   - Organizes output files in a structured directory format

**Resource requirements:**
- 48 CPUs, 256GB memory, 72-hour runtime
- Requires specialized annotation tools (BRAKER, BLAST, InterProScan, etc.)

**Key outputs:**
- Gene predictions in GFF3 format
- Protein sequences in FASTA format
- Functional annotations with GO terms and domain information
- Comprehensive reports on tissue-specific transcriptomes

## Key Outputs

### Task 1 Outputs
- `dbg/task_1_3_1/`: DBG assembly and graph for Task 1.3.1
- `dbg/task_1_3_2/k_35/` and `dbg/task_1_3_2/k_45/`: Assemblies with different k values
- `dbg/mers/` and `olc/mers/`: Assemblies of MERS datasets
- `spades/mers/`: SPAdes assembly results
- `evaluations/`: QUAST evaluation results for all assemblies

### Task 2 Outputs
- `assembly/*.p_ctg.fasta`: Primary contig assemblies from different approaches
- `evaluation/quast/`: Basic assembly metrics
- `evaluation/busco/`: Gene completeness assessment
- `evaluation/merqury/`: K-mer analysis and QV scores
- `evaluation/mis_assembly/`: Analysis of potential misassemblies
- `improved/improved_assembly.fasta`: Final improved assembly
- `report/assembly_evaluation_report.txt`: Comprehensive evaluation report
- `annotation/scincus_mitranus_genes.gff3`: Gene predictions
- `annotation/scincus_mitranus_proteins.faa`: Protein sequences
- `annotation/functional/`: Functional annotation results
- `report/annotation/final_annotation_report.txt`: Genome annotation report

## Interpretation of Key Results

### Assembly Metrics
- **N50**: Length where 50% of the assembly is in contigs of this size or larger
- **QV Score**: Base-level accuracy measure
  - QV 40 = 99.99% accuracy (1 error per 10,000 bases)
  - QV 50 = 99.999% accuracy (1 error per 100,000 bases)
- **BUSCO Completeness**: Percentage of conserved genes found in the assembly

### Performance Comparison
- The DBG assembler shows superior performance on short, error-free reads
- The OLC assembler handles long reads with errors more effectively
- The integrated approach (PacBio+Nanopore+Hi-C) produces the most contiguous assembly
- Assembly improvement techniques increase QV scores by reducing error rates
- Annotation with IsoSeq RNA data improves gene model accuracy, especially for splice variants

## Advanced Usage Notes

### K-mer Size Selection
- **Small k-mer size** (e.g., k=21-31): Better handles error-prone reads but may collapse repeats
- **Large k-mer size** (e.g., k=45-55): Better resolves repeats but may fragment error-containing regions
- For error-containing datasets, k-mer frequencies can be filtered with `--min-freq` parameter

### Memory Considerations
- For full-sized datasets, increase memory allocation (192GB+ recommended)
- When processing the complete *Scincus mitranus* dataset, increase runtime to 48+ hours
- Annotation with BRAKER and InterProScan requires significant memory (256GB recommended)

### Output Visualization
- GFA files can be directly loaded into Bandage for visualization
- Use the `make_gfa.py` script to convert contig layouts to GFA format
- Gene annotations can be visualized with genome browsers like IGV or JBrowse

## Acknowledgments

This work utilized Claude 3.7 Sonnet (Anthropic) for implementation guidance, debugging assistance, and report drafting. All final code implementation, analysis, and interpretations were validated through systematic evaluation of assembly outputs.

## License

MIT License