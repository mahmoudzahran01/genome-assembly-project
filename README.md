# Genome Assembly Project

This repository contains implementations of genome assembly algorithms and analysis of assembly results for both model datasets and the *Scincus mitranus* genome.

## Project Overview

This project consists of two main tasks:

1. Implementation and evaluation of De Bruijn Graph (DBG) and Overlap-Layout-Consensus (OLC) assembly algorithms on simulated datasets
2. Assembly of the *Scincus mitranus* genome using multiple sequencing technologies (PacBio HiFi, Oxford Nanopore, and Hi-C)

## Repository Structure

```
.
├── README.md                   # Project documentation
├── assembly_runner.py          # Assembly execution wrapper
├── dbg_assembler.py            # De Bruijn Graph assembler implementation
├── make_gfa.py                 # GFA format generator for visualizations
├── olc_assembler.py            # Overlap-Layout-Consensus assembler implementation
├── run_task1.py                # Main script for Task 1 assembly
├── data/                       # Data directory
│   └── task1/                  # Task 1 simulated datasets
│       ├── synthetic_1/        # Read and reference files for synthetic dataset 1
│       │   ├── reads_*.fastq   # Raw sequencing reads
│       │   └── reference_*.fasta # Reference sequences
│       └── synthetic_2/        # Read and reference files for synthetic dataset 2
│           ├── GCF_000901155.1_ViralProj183710_genomic.fna # Reference genome
│           ├── no_error_*.fastq # Error-free simulated reads 
│           └── *_hq_*.fastq    # Reads with simulated errors
├── evaluations/                # Evaluation results
│   ├── mers/                   # Mer-based evaluations
│   │   ├── dbg/                # DBG assembly evaluations
│   │   │   ├── dbg_*_no_error/ # Error-free dataset evaluations
│   │   │   ├── dbg_*_with_error/ # Error-containing dataset evaluations
│   │   │   └── spades_*/       # SPAdes comparisons
│   │   └── olc/                # OLC assembly evaluations
│   └── synthetic_1/            # Evaluations on synthetic dataset 1
│       └── dbg/                # DBG assembly evaluations and comparisons
│           ├── comparison.*    # K-mer size comparison files
│           ├── k_35/           # K=35 assembly results
│           └── k_45/           # K=45 assembly results
├── results/                    # Output results directory
│   ├── task1/                  # Task 1 results
│   │   └── spades/             # SPAdes assembly results
│   └── task2/                  # Task 2 results
│       ├── assembly/           # Assembled genomes
│       │   ├── scincus_mitranus.bp.p_ctg.fasta           # PacBio-only assembly
│       │   ├── scincus_mitranus_nanopore.bp.p_ctg.fasta  # PacBio+Nanopore assembly
│       │   ├── scincus_mitranus_hic.hic.p_ctg.fasta      # PacBio+Hi-C assembly
│       │   └── scincus_mitranus_full.hic.p_ctg.fasta     # Integrated assembly
│       ├── evaluation/         # Quality assessment
│       │   ├── busco/          # Completeness assessment
│       │   ├── merqury/        # K-mer based quality evaluation
│       │   └── quast/          # Assembly statistics
│       ├── improvement/        # Improved assemblies
│       │   ├── correction_logs/  # Read correction logs
│       │   └── improved_assembly.fasta  # Final improved assembly
│       └── reports/            # Analysis reports
│           └── *.txt           # Various summary and analysis files
├── scripts/                    # Helper scripts
│   └── utility/                # Utility scripts
├── task2/                      # Task 2 scripts
│   ├── final_task2_2.sh        # Task 2.2 assembly script
│   ├── process_all_datasets.sh # Process and assemble datasets for Task 1
│   └── task2_3_assembly_improvement_subsampled.sh # Assembly improvement pipeline
└── visualizations/             # Visualizations and figures
    ├── *.png                   # Assembly graph visualizations
    ├── task1/                  # Task 1 visualizations
    │   ├── busco/              # BUSCO visualization
    │   ├── merqury/            # Merqury quality assessment
    │   ├── mis_assembly/       # Misassembly analysis
    │   ├── plots/              # Various summary plots
    │   │   └── *.pdf           # Plot files (quality, coverage, etc.)
    │   └── quast/              # QUAST results
    └── task2/                  # Task 2 visualizations
        └── */                  # Various visualization categories
```

## Requirements

- Python 3.8+
- BioPython
- Matplotlib
- NetworkX
- NumPy
- QUAST
- Bandage
- Hifiasm (v0.19.5+)
- Kraken2
- Trimmomatic
- Porechop

## Task 1: Assembly Algorithms

### Implementation

Two fundamentally different genome assembly approaches were implemented:

- **De Bruijn Graph (DBG) Assembly**: Breaks reads into k-mers and constructs a directed graph where nodes represent (k-1)-mers and edges represent consecutive overlaps
- **Overlap-Layout-Consensus (OLC) Assembly**: Uses whole reads as nodes and computes pairwise overlaps to construct an overlap graph

### Running the Assembly

To run the assembly on all datasets for Task 1:

```bash
./task2/process_all_datasets.sh
```

This script is the primary workflow for Task 1, processing all error-free and error-containing datasets using both DBG and OLC assembly methods. It:
- Runs the DBG assembler on HiSeq and ONT datasets (with and without errors)
- Runs the OLC assembler on all datasets with varying overlap parameters
- Executes comparative analysis against SPAdes assemblies
- Generates evaluation metrics for all assemblies

For a specific analysis of one dataset:

```bash
python run_task1.py --input data/task1/synthetic_2/reads_hiseq_5k.fastq --output results/task1/ --algorithm dbg --k 40
```

Parameters:

- `--input`: Input FASTQ file
- `--output`: Output directory
- `--algorithm`: Assembly algorithm (dbg or olc)
- `--k`: k-mer size (for DBG) or minimum overlap (for OLC)
- `--min_freq`: Minimum k-mer frequency (default: 2)
- `--min_overlap`: Minimum overlap length for OLC (default: 30)

### Analyzing Results

Assembly graphs can be visualized using Bandage with the GFA files generated by `make_gfa.py`. Example visualizations are available in the `visualizations/` directory.

The impact of k-mer size on assembly quality can be explored in the `evaluations/synthetic_1/dbg/` directory, which includes comparisons between k=35 and k=45 assemblies.

## Task 2: Scincus mitranus Genome Assembly

### Assembly Strategies

The *Scincus mitranus* genome was assembled using a multi-phase approach:

1. **Task 2.1: Initial Read Processing and Quality Control**
   - Quality assessment of raw PacBio HiFi, Oxford Nanopore, and Hi-C reads
   - Adapter removal and read filtering
   - Contamination screening using Kraken2
   - Read subsampling for optimal coverage depth

2. **Task 2.2: Genome Assembly with Different Technologies**
   Four assembly strategies were evaluated:
   - PacBio-only assembly
   - PacBio + Nanopore assembly
   - PacBio + Hi-C assembly
   - Full integrated assembly (combining all three data types)

3. **Task 2.3: Assembly Improvement**
   - Error correction for long reads
   - Quality assessment and comparison of different assembly versions
   - Final polishing and gap filling

### Running the Assembly

To run the main assembly task:

```bash
./task2/final_task2_2.sh
```

This script orchestrates the PacBio HiFi assembly with Hifiasm, testing different combinations of input data.

### Assembly Improvement

To run the assembly improvement pipeline:

```bash
./task2/task2_3_assembly_improvement_subsampled.sh
```

This script implements a multi-step improvement process:

1. Quality control and filtering
2. Error correction for long reads
3. Quality assessment and evaluation
4. Visualization of assembly graphs

## Key Findings

- The DBG approach performed excellently on error-free data but struggled with error-containing ONT reads
- The OLC approach showed greater resilience to errors in long-read data
- For *Scincus mitranus* assembly, the PacBio+Nanopore integration provided the best results, with 47% fewer contigs and 28% higher N50 compared to PacBio-only assembly
- Hi-C data alone provided limited benefits at the contig level without additional scaffolding

Assembly statistics and detailed analyses are available in the `results/task2/reports/` directory.

## Visualizations

Key visualizations include:

- `bandage_1_3_1.png`: De Bruijn graph visualization for task 1.3.1
- `k35_assembly_graph.png`: Assembly graph with k=35
- `k45_assembly_graph.png`: Assembly graph with k=45
- Additional plots and quality assessment visualizations in the `visualizations/` directory

## Acknowledgments

This work utilized Claude 3.7 Sonnet (Anthropic) for implementation guidance, debugging assistance, and report drafting. All final code implementation, analysis, and interpretations were validated through systematic evaluation of assembly outputs against reference sequences and thorough cross-referencing of metrics against actual output files.

## License

MIT License