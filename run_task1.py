#!/usr/bin/env python3
"""
Bioinformatics Algorithms - Assignment 2
Task 1 Runner Script: Genome Assembly and Evaluation

This script runs all components of Task 1:
1. DBG and OLC assemblers on synthetic datasets
2. Comparison of assembly results with QUAST
3. Application to MERS-CoV data
4. Comparison with professional assemblers

Usage:
    python run_task1.py --output_dir path/to/output
"""

import os
import argparse
import subprocess
import logging
import sys
import shutil
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("task1_runner.log")
    ]
)
logger = logging.getLogger(__name__)

def run_command(cmd, desc=None):
    """Run a shell command and log its output."""
    if desc:
        logger.info(f"Running: {desc}")
    logger.debug(f"Command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True)
        logger.debug(f"Output: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}")
        logger.error(f"stderr: {e.stderr}")
        logger.error(f"stdout: {e.stdout}")
        return False

def ensure_dir(directory):
    """Ensure a directory exists, create if it doesn't."""
    if not os.path.exists(directory):
        os.makedirs(directory)
        logger.info(f"Created directory: {directory}")

def main():
    parser = argparse.ArgumentParser(description="Run Task 1 of Assignment 2")
    parser.add_argument("--output_dir", required=True, help="Output directory for all results")
    parser.add_argument("--data_dir", default="task1/data", help="Directory containing input data")
    parser.add_argument("--skip_dbg", action="store_true", help="Skip DBG assembly")
    parser.add_argument("--skip_olc", action="store_true", help="Skip OLC assembly")
    parser.add_argument("--skip_spades", action="store_true", help="Skip SPAdes assembly")
    parser.add_argument("--skip_evaluation", action="store_true", help="Skip evaluation with QUAST")
    
    args = parser.parse_args()
    
    # Create necessary directories
    output_dir = args.output_dir
    ensure_dir(output_dir)
    
    dbg_dir = os.path.join(output_dir, "dbg")
    olc_dir = os.path.join(output_dir, "olc")
    spades_dir = os.path.join(output_dir, "spades")
    eval_dir = os.path.join(output_dir, "evaluations")
    
    data_dir = args.data_dir
    synthetic_1_dir = os.path.join(data_dir, "synthetic_1")
    synthetic_2_dir = os.path.join(data_dir, "synthetic_2")
    
    ensure_dir(dbg_dir)
    ensure_dir(olc_dir)
    ensure_dir(spades_dir)
    ensure_dir(eval_dir)
    
    # Define data paths
    datasets = {
        "synthetic_1": {
            "reads_r": os.path.join(synthetic_1_dir, "reads_r.fastq"),
            "reads_b": os.path.join(synthetic_1_dir, "reads_b.fastq"),
            "reference_r": os.path.join(synthetic_1_dir, "reference_r.fasta"),
            "reference_b": os.path.join(synthetic_1_dir, "reference_b.fasta")
        },
        "mers": {
            "no_error_hiseq": os.path.join(synthetic_2_dir, "no_error_reads_hiseq_5k.fastq"),
            "with_error_hiseq": os.path.join(synthetic_2_dir, "reads_hiseq_5k.fastq"),
            "no_error_ont": os.path.join(synthetic_2_dir, "no_error_ont_hq_50x.fastq"),
            "with_error_ont": os.path.join(synthetic_2_dir, "ont_hq_50x.fastq"),
            "reference": os.path.join(synthetic_2_dir, "GCF_000901155.1_ViralProj183710_genomic.fna")
        }
    }
    
    # Task 1.1: Run DBG Assembler
    if not args.skip_dbg:
        logger.info("=== Running Task 1.1: DBG Assembly ===")
        
        # Task 1.3.1: Construct assembly graph for reads_b.fastq with k=40
        dbg_task_1_3_1_dir = os.path.join(dbg_dir, "task_1_3_1")
        ensure_dir(dbg_task_1_3_1_dir)
        run_command(
            f"python dbg_assembler.py --fastq {datasets['synthetic_1']['reads_b']} --k 40 "
            f"--output {dbg_task_1_3_1_dir}",
            "DBG Assembly for Task 1.3.1"
        )
        
        # Task 1.3.2: DBG on reads_r.fastq with k=35 and k=45
        dbg_task_1_3_2_dir = os.path.join(dbg_dir, "task_1_3_2")
        ensure_dir(dbg_task_1_3_2_dir)
        
        for k in [35, 45]:
            k_dir = os.path.join(dbg_task_1_3_2_dir, f"k_{k}")
            ensure_dir(k_dir)
            run_command(
                f"python dbg_assembler.py --fastq {datasets['synthetic_1']['reads_r']} --k {k} "
                f"--output {k_dir}",
                f"DBG Assembly for Task 1.3.2 with k={k}"
            )
        
        # Task 1.3.3: DBG on MERS datasets
        dbg_mers_dir = os.path.join(dbg_dir, "mers")
        ensure_dir(dbg_mers_dir)
        
        mers_datasets = [
            ("hiseq_no_error", datasets["mers"]["no_error_hiseq"]),
            ("hiseq_with_error", datasets["mers"]["with_error_hiseq"]),
            ("ont_no_error", datasets["mers"]["no_error_ont"]),
            ("ont_with_error", datasets["mers"]["with_error_ont"])
        ]
        
        for name, fastq in mers_datasets:
            dataset_dir = os.path.join(dbg_mers_dir, name)
            ensure_dir(dataset_dir)
            
            # Use different k values for error-containing reads
            if "with_error" in name:
                k_values = [21, 25, 31]
            else:
                k_values = [31]
            
            for k in k_values:
                run_command(
                    f"python dbg_assembler.py --fastq {fastq} --k {k} --min-freq 3 "
                    f"--output {dataset_dir}",
                    f"DBG Assembly for MERS {name} with k={k}"
                )
    
    # Task 1.2: Run OLC Assembler
    if not args.skip_olc:
        logger.info("=== Running Task 1.2: OLC Assembly ===")
        
        # Task 1.3.3: OLC on MERS datasets
        olc_mers_dir = os.path.join(olc_dir, "mers")
        ensure_dir(olc_mers_dir)
        
        mers_datasets = [
            ("hiseq_no_error", datasets["mers"]["no_error_hiseq"]),
            ("hiseq_with_error", datasets["mers"]["with_error_hiseq"]),
            ("ont_no_error", datasets["mers"]["no_error_ont"]),
            ("ont_with_error", datasets["mers"]["with_error_ont"])
        ]
        
        for name, fastq in mers_datasets:
            dataset_dir = os.path.join(olc_mers_dir, name)
            ensure_dir(dataset_dir)
            
            # Different min_overlap values for different datasets
            if "ont" in name:
                min_overlap = 50
            else:
                min_overlap = 20
            
            run_command(
                f"python olc_assembler.py --fastq {fastq} --min-overlap {min_overlap} "
                f"--output {dataset_dir}",
                f"OLC Assembly for MERS {name} with min_overlap={min_overlap}"
            )
    
    # Task 1.3.4: Run professional assemblers (SPAdes)
    if not args.skip_spades:
        logger.info("=== Running Task 1.3.4: SPAdes Assembly ===")
        
        spades_mers_dir = os.path.join(spades_dir, "mers")
        ensure_dir(spades_mers_dir)
        
        # Run SPAdes on HiSeq datasets
        spades_hiseq_no_error_dir = os.path.join(spades_mers_dir, "hiseq_no_error")
        ensure_dir(spades_hiseq_no_error_dir)
        run_command(
            f"spades.py -s {datasets['mers']['no_error_hiseq']} "
            f"-o {spades_hiseq_no_error_dir} --isolate",
            "SPAdes Assembly for MERS HiSeq no error"
        )
        
        spades_hiseq_with_error_dir = os.path.join(spades_mers_dir, "hiseq_with_error")
        ensure_dir(spades_hiseq_with_error_dir)
        run_command(
            f"spades.py -s {datasets['mers']['with_error_hiseq']} "
            f"-o {spades_hiseq_with_error_dir} --isolate",
            "SPAdes Assembly for MERS HiSeq with error"
        )
        
        # Try SPAdes for ONT data (may not work well)
        spades_ont_no_error_dir = os.path.join(spades_mers_dir, "ont_no_error")
        ensure_dir(spades_ont_no_error_dir)
        run_command(
            f"spades.py --only-assembler -s {datasets['mers']['no_error_ont']} "
            f"-o {spades_ont_no_error_dir} --isolate",
            "SPAdes Assembly for MERS ONT no error"
        )
        
        spades_ont_with_error_dir = os.path.join(spades_mers_dir, "ont_with_error")
        ensure_dir(spades_ont_with_error_dir)
        run_command(
            f"spades.py --only-assembler -s {datasets['mers']['with_error_ont']} "
            f"-o {spades_ont_with_error_dir} --isolate",
            "SPAdes Assembly for MERS ONT with error"
        )
    
    # Evaluate assemblies with QUAST
    if not args.skip_evaluation:
        logger.info("=== Running Assembly Evaluation with QUAST ===")
        
        # Task 1.3.2: Evaluate DBG on reads_r.fastq with k=35 and k=45
        eval_task_1_3_2_dir = os.path.join(eval_dir, "task_1_3_2")
        ensure_dir(eval_task_1_3_2_dir)
        
        for k in [35, 45]:
            k_eval_dir = os.path.join(eval_task_1_3_2_dir, f"k_{k}")
            ensure_dir(k_eval_dir)
            
            contigs_file = os.path.join(dbg_dir, "task_1_3_2", f"k_{k}", "contigs_k{k}.fasta")
            reference_file = datasets["synthetic_1"]["reference_r"]
            
            run_command(
                f"quast.py -o {k_eval_dir} -r {reference_file} {contigs_file}",
                f"QUAST evaluation for DBG assembly of reads_r.fastq with k={k}"
            )
        
        # Task 1.3.3: Evaluate DBG and OLC on MERS datasets
        eval_mers_dir = os.path.join(eval_dir, "mers")
        ensure_dir(eval_mers_dir)
        
        mers_datasets = [
            ("hiseq_no_error", "no_error_hiseq"),
            ("hiseq_with_error", "with_error_hiseq"),
            ("ont_no_error", "no_error_ont"),
            ("ont_with_error", "with_error_ont")
        ]
        
        for name, data_key in mers_datasets:
            # Evaluate DBG assemblies
            dbg_eval_dir = os.path.join(eval_mers_dir, f"dbg_{name}")
            ensure_dir(dbg_eval_dir)
            
            dbg_contigs_file = os.path.join(dbg_mers_dir, name, "contigs_k31.fasta")
            if os.path.exists(dbg_contigs_file):
                run_command(
                    f"quast.py -o {dbg_eval_dir} -r {datasets['mers']['reference']} {dbg_contigs_file} --min-contig 0",
                    f"QUAST evaluation for DBG assembly of MERS {name}"
                )
            
            # Evaluate OLC assemblies
            olc_eval_dir = os.path.join(eval_mers_dir, f"olc_{name}")
            ensure_dir(olc_eval_dir)
            
            min_overlap = 50 if "ont" in name else 20
            olc_contigs_file = os.path.join(olc_mers_dir, name, f"contigs_olc_min{min_overlap}.fasta")
            if os.path.exists(olc_contigs_file):
                run_command(
                    f"quast.py -o {olc_eval_dir} -r {datasets['mers']['reference']} {olc_contigs_file}",
                    f"QUAST evaluation for OLC assembly of MERS {name}"
                )
            
            # Evaluate SPAdes assemblies
            spades_eval_dir = os.path.join(eval_mers_dir, f"spades_{name}")
            ensure_dir(spades_eval_dir)
            
            spades_scaffolds_file = os.path.join(spades_mers_dir, name, "scaffolds.fasta")
            if os.path.exists(spades_scaffolds_file):
                run_command(
                    f"quast.py -o {spades_eval_dir} -r {datasets['mers']['reference']} {spades_scaffolds_file}",
                    f"QUAST evaluation for SPAdes assembly of MERS {name}"
                )
    
    logger.info("Task 1 completed successfully!")

if __name__ == "__main__":
    main()