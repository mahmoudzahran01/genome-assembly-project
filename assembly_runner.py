#!/usr/bin/env python3
"""
Assembly Runner

This module provides a unified interface for running different assembly algorithms
and evaluating the results.
"""

import os
import argparse
import logging
import subprocess
import json
from typing import List, Dict

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Import our assemblers
from dbg_assembler import run_dbg_assembly

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def json_serializable(obj):
    """
    Convert numpy types to Python native types for JSON serialization.
    
    Args:
        obj: Object to convert
        
    Returns:
        JSON serializable object
    """
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def run_quast(contigs_file, reference_file, output_dir):
    """
    Run QUAST for assembly evaluation.
    
    Args:
        contigs_file: Path to contigs FASTA file
        reference_file: Path to reference genome FASTA file
        output_dir: Directory for QUAST output
        
    Returns:
        Dictionary of evaluation metrics
    """
    quast_cmd = [
        "quast.py",
        "-o", output_dir,
        "-r", reference_file,
        contigs_file
    ]
    
    try:
        logger.info("Running QUAST: {}".format(' '.join(quast_cmd)))
        subprocess.run(quast_cmd, check=True, capture_output=True, text=True)
        
        # Parse QUAST report
        report_file = os.path.join(output_dir, "report.tsv")
        if os.path.exists(report_file):
            metrics = {}
            with open(report_file, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:  # Skip header
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        metric_name = parts[0]
                        value = parts[1]
                        try:
                            # Convert numeric values
                            metrics[metric_name] = float(value) if '.' in value else int(value)
                        except ValueError:
                            metrics[metric_name] = value
            
            return metrics
        else:
            logger.error("QUAST report not found: {}".format(report_file))
            return {"error": "QUAST report not found"}
    
    except subprocess.CalledProcessError as e:
        logger.error("QUAST failed: {}".format(e.stderr))
        return {"error": e.stderr}
    except Exception as e:
        logger.error("Error running QUAST: {}".format(e))
        return {"error": str(e)}


def compare_assemblies(results_list, output_file):
    """
    Compare multiple assembly results and generate visualizations.
    
    Args:
        results_list: List of assembly result dictionaries
        output_file: Path to output file (PDF)
    """
    # Extract key metrics for comparison
    comparison_data = []
    
    for results in results_list:
        k = results.get("k", "unknown")
        metrics = results.get("assembly", {})
        quast = results.get("quast", {})
        
        data_point = {
            "k": k,
            "num_contigs": metrics.get("num_contigs", 0),
            "total_length": metrics.get("total_length", 0),
            "n50": metrics.get("n50", 0),
            "n90": metrics.get("n90", 0),
            "largest_contig": metrics.get("max_length", 0),
            "gc_content": metrics.get("gc_content", 0),
        }
        
        # Add QUAST metrics if available
        if quast:
            data_point.update({
                "genome_fraction": quast.get("Genome fraction (%)", 0),
                "misassemblies": quast.get("# misassemblies", 0),
                "mismatches": quast.get("# mismatches per 100 kbp", 0),
                "indels": quast.get("# indels per 100 kbp", 0)
            })
        
        comparison_data.append(data_point)
    
    if not comparison_data:
        logger.error("No data for comparison")
        return
    
    # Create DataFrame
    df = pd.DataFrame(comparison_data)
    df = df.sort_values("k")
    
    # Create visualizations
    plt.figure(figsize=(12, 10))
    
    # Plot 1: Number of contigs vs k
    plt.subplot(2, 2, 1)
    plt.plot(df["k"], df["num_contigs"], 'o-')
    plt.xlabel("k-mer size")
    plt.ylabel("Number of contigs")
    plt.title("Number of contigs vs k-mer size")
    plt.grid(True)
    
    # Plot 2: N50 vs k
    plt.subplot(2, 2, 2)
    plt.plot(df["k"], df["n50"], 'o-')
    plt.xlabel("k-mer size")
    plt.ylabel("N50")
    plt.title("N50 vs k-mer size")
    plt.grid(True)
    
    # Plot 3: Total assembly length vs k
    plt.subplot(2, 2, 3)
    plt.plot(df["k"], df["total_length"], 'o-')
    plt.xlabel("k-mer size")
    plt.ylabel("Total assembly length")
    plt.title("Assembly length vs k-mer size")
    plt.grid(True)
    
    # Plot 4: Largest contig vs k
    plt.subplot(2, 2, 4)
    plt.plot(df["k"], df["largest_contig"], 'o-')
    plt.xlabel("k-mer size")
    plt.ylabel("Largest contig length")
    plt.title("Largest contig vs k-mer size")
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(output_file)
    logger.info("Comparison plots saved to {}".format(output_file))
    
    # Save tabular data
    csv_file = output_file.replace(".pdf", ".csv")
    df.to_csv(csv_file, index=False)
    logger.info("Comparison data saved to {}".format(csv_file))


def run_task1_3(input_dir, output_dir, reference_file):
    """
    Run analysis for Task 1.3 (compare different k values).
    
    Args:
        input_dir: Directory containing input FASTQ files
        output_dir: Directory for output files
        reference_file: Path to reference genome FASTA file
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # k-mer sizes to test
    k_values = [35, 45]  # As specified in the assignment
    
    # Input FASTQ file
    fastq_file = os.path.join(input_dir, "reads_r.fastq")
    
    # Run assembly for each k value
    results_list = []
    
    for k in k_values:
        logger.info("Running assembly with k={}".format(k))
        
        # Create k-specific output directory
        k_output_dir = os.path.join(output_dir, "k_{}".format(k))
        os.makedirs(k_output_dir, exist_ok=True)
        
        # Run DBG assembly
        assembly_results = run_dbg_assembly(
            fastq_file=fastq_file,
            k=k,
            min_frequency=2,
            output_dir=k_output_dir,
            reference_file=reference_file
        )
        
        # Run QUAST evaluation
        contigs_file = assembly_results["files"]["contigs"]
        quast_dir = os.path.join(k_output_dir, "quast")
        quast_results = run_quast(contigs_file, reference_file, quast_dir)
        
        # Add to results list
        assembly_results["k"] = k
        assembly_results["quast"] = quast_results
        results_list.append(assembly_results)
        
        # Save individual results
        results_file = os.path.join(k_output_dir, "complete_results.json")
        with open(results_file, 'w') as f:
            json.dump(assembly_results, f, indent=2, default=json_serializable)
    
    # Compare assemblies
    comparison_file = os.path.join(output_dir, "comparison.pdf")
    compare_assemblies(results_list, comparison_file)


def main():
    """Main entry point for the program."""
    parser = argparse.ArgumentParser(description="Genome Assembly Runner")
    
    # Main options
    parser.add_argument("--task", choices=["dbg", "olc", "task1.3"], required=True,
                        help="Task to run")
    parser.add_argument("--input", required=True, help="Input directory or file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--reference", help="Reference genome FASTA file")
    
    # DBG options
    parser.add_argument("--k", type=int, default=31, help="k-mer size for DBG (default: 31)")
    parser.add_argument("--min-freq", type=int, default=1, 
                        help="Minimum k-mer frequency for DBG (default: 1)")
    
    # OLC options
    parser.add_argument("--min-overlap", type=int, default=20, 
                        help="Minimum overlap length for OLC (default: 20)")
    
    # Other options
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    try:
        if args.task == "dbg":
            # Run DBG assembly on a single file
            run_dbg_assembly(
                fastq_file=args.input,
                k=args.k,
                min_frequency=args.min_freq,
                output_dir=args.output,
                reference_file=args.reference
            )
            
        elif args.task == "task1.3":
            # Run Task 1.3 analysis
            run_task1_3(
                input_dir=args.input,
                output_dir=args.output,
                reference_file=args.reference
            )
            
        elif args.task == "olc":
            # Placeholder for OLC assembly
            logger.warning("OLC assembly not yet implemented")
            
        logger.info("Task {} completed successfully".format(args.task))
        
    except Exception as e:
        logger.error("Task failed: {}".format(e))
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    main()