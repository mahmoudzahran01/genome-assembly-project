#!/usr/bin/env python3
"""
Overlap-Layout-Consensus (OLC) Genome Assembler

This module implements an OLC approach to genome assembly,
taking FASTQ files as input and generating contigs as FASTA output.
"""

import os
import json
import argparse
import logging
from collections import defaultdict
import heapq

import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio import Align
from Bio.Align import PairwiseAligner

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


def read_fastq_sequences(fastq_file):
    """
    Read sequences from a single FASTQ file.
    
    Args:
        fastq_file: Path to FASTQ file
        
    Returns:
        Dictionary mapping read IDs to sequences
    """
    sequences = {}
    try:
        for record in SeqIO.parse(fastq_file, "fastq"):
            sequences[record.id] = str(record.seq)
    except Exception as e:
        logger.error("Error reading {}: {}".format(fastq_file, e))
        
    logger.info("Read {} sequences from {}".format(len(sequences), fastq_file))
    return sequences


def compute_overlap(seq1, seq2, min_overlap):
    """
    Compute the overlap between the suffix of seq1 and prefix of seq2.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        min_overlap: Minimum required overlap length
        
    Returns:
        Tuple of (overlap_length, score) or (0, 0) if no sufficient overlap
    """
    # Start with the minimum overlap and try increasingly larger overlaps
    max_overlap = min(len(seq1), len(seq2))
    
    # Only check overlaps >= min_overlap
    for overlap_len in range(min_overlap, max_overlap + 1):
        # Check if suffix of seq1 matches prefix of seq2
        if seq1[-overlap_len:] == seq2[:overlap_len]:
            # Calculate a score based on overlap length and quality
            score = overlap_len
            return (overlap_len, score)
    
    return (0, 0)


def compute_pairwise_overlaps(sequences, min_overlap, max_reads=None):
    """
    Compute all pairwise overlaps between sequences.
    
    Args:
        sequences: Dictionary mapping read IDs to sequences
        min_overlap: Minimum required overlap length
        max_reads: Maximum number of reads to process (for testing)
        
    Returns:
        List of (read_id1, read_id2, overlap_length, score) tuples
    """
    overlaps = []
    read_ids = list(sequences.keys())
    
    # Limit number of reads if specified (for testing)
    if max_reads and max_reads < len(read_ids):
        read_ids = read_ids[:max_reads]
        logger.info("Limited to {} reads for overlap computation".format(max_reads))
    
    total_pairs = len(read_ids) * (len(read_ids) - 1)
    logger.info("Computing overlaps for {} sequence pairs...".format(total_pairs))
    
    # Compute overlaps for all pairs
    for i, id1 in enumerate(read_ids):
        if i % 100 == 0:
            logger.info("Progress: {}/{} reads processed".format(i, len(read_ids)))
        
        seq1 = sequences[id1]
        
        for id2 in read_ids:
            if id1 != id2:  # Don't compare a sequence to itself
                seq2 = sequences[id2]
                
                # Compute overlap
                overlap_len, score = compute_overlap(seq1, seq2, min_overlap)
                
                if overlap_len >= min_overlap:
                    overlaps.append((id1, id2, overlap_len, score))
    
    logger.info("Found {} overlaps of length >= {}".format(len(overlaps), min_overlap))
    return overlaps


def build_overlap_graph(sequences, overlaps):
    """
    Build an overlap graph from the computed overlaps.
    
    Args:
        sequences: Dictionary mapping read IDs to sequences
        overlaps: List of (read_id1, read_id2, overlap_length, score) tuples
        
    Returns:
        NetworkX directed graph
    """
    graph = nx.DiGraph()
    
    # Add nodes
    for read_id in sequences.keys():
        graph.add_node(read_id, sequence=sequences[read_id])
    
    # Add edges
    for id1, id2, overlap_len, score in overlaps:
        graph.add_edge(id1, id2, overlap=overlap_len, score=score)
    
    logger.info("Built overlap graph with {} nodes and {} edges".format(
        graph.number_of_nodes(), graph.number_of_edges()))
    
    return graph


def find_contigs_from_paths(graph):
    """
    Find contigs by identifying non-branching paths in the graph.
    
    Args:
        graph: NetworkX directed graph
        
    Returns:
        List of paths, where each path is a list of read IDs
    """
    contigs = []
    visited = set()
    
    # Find all source nodes (nodes with no incoming edges)
    source_nodes = [node for node in graph.nodes() if graph.in_degree(node) == 0]
    
    # If no source nodes, find nodes with more outgoing than incoming edges
    if not source_nodes:
        source_nodes = [node for node in graph.nodes() 
                       if graph.out_degree(node) > graph.in_degree(node)]
    
    # If still no suitable start nodes, use any node with outgoing edges
    if not source_nodes:
        source_nodes = [node for node in graph.nodes() if graph.out_degree(node) > 0]
    
    # Start from each source node and find paths
    for start_node in source_nodes:
        if start_node in visited:
            continue
        
        # Start a new path
        path = [start_node]
        visited.add(start_node)
        
        # Follow the path as long as there is exactly one outgoing edge
        current = start_node
        while True:
            successors = list(graph.successors(current))
            if len(successors) == 1 and successors[0] not in path:
                # Continue the path
                next_node = successors[0]
                path.append(next_node)
                visited.add(next_node)
                current = next_node
            else:
                # End of non-branching path
                break
        
        # Only keep paths with multiple nodes
        if len(path) > 1:
            contigs.append(path)
    
    # Look for additional paths starting from nodes not yet visited
    for node in graph.nodes():
        if node not in visited and graph.out_degree(node) > 0:
            path = [node]
            visited.add(node)
            
            current = node
            while True:
                successors = list(graph.successors(current))
                if len(successors) == 1 and successors[0] not in path:
                    next_node = successors[0]
                    path.append(next_node)
                    visited.add(next_node)
                    current = next_node
                else:
                    break
            
            if len(path) > 1:
                contigs.append(path)
    
    # Find potential cycles
    cycle_basis = list(nx.simple_cycles(graph))
    for cycle in cycle_basis:
        if len(cycle) > 1 and all(node not in visited for node in cycle):
            for node in cycle:
                visited.add(node)
            contigs.append(cycle)
    
    logger.info("Found {} contigs (non-branching paths)".format(len(contigs)))
    return contigs


def generate_consensus_sequence(path, graph, sequences):
    """
    Generate a consensus sequence from a path of reads.
    
    Args:
        path: List of read IDs representing a contig path
        graph: NetworkX directed graph with overlap information
        sequences: Dictionary mapping read IDs to sequences
        
    Returns:
        Consensus sequence for the contig
    """
    if not path:
        return ""
    
    # Start with the first read
    consensus = sequences[path[0]]
    
    # Add subsequent reads, accounting for overlaps
    for i in range(1, len(path)):
        prev_id = path[i-1]
        curr_id = path[i]
        
        # Get the overlap length from the graph
        overlap_len = graph[prev_id][curr_id]["overlap"]
        
        # Add the non-overlapping part of the current read
        curr_seq = sequences[curr_id]
        consensus += curr_seq[overlap_len:]
    
    return consensus


def calculate_assembly_metrics(contigs):
    """
    Calculate basic assembly metrics.
    
    Args:
        contigs: List of assembled contig sequences
        
    Returns:
        Dictionary of metrics
    """
    if not contigs:
        return {
            "num_contigs": 0,
            "total_length": 0,
            "min_length": 0,
            "max_length": 0,
            "mean_length": 0,
            "n50": 0,
            "n90": 0,
            "gc_content": 0
        }
    
    # Calculate contig lengths
    lengths = [len(contig) for contig in contigs]
    total_length = sum(lengths)
    
    # Calculate N50 and N90
    sorted_lengths = sorted(lengths, reverse=True)
    cumulative_length = 0
    n50, n90 = 0, 0
    
    for length in sorted_lengths:
        cumulative_length += length
        if not n50 and cumulative_length >= total_length * 0.5:
            n50 = length
        if not n90 and cumulative_length >= total_length * 0.9:
            n90 = length
            break
    
    # Calculate GC content
    gc_count = sum(contig.count('G') + contig.count('C') for contig in contigs)
    gc_content = (gc_count / total_length) * 100 if total_length > 0 else 0
    
    metrics = {
        "num_contigs": len(contigs),
        "total_length": total_length,
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "mean_length": total_length / len(contigs) if contigs else 0,
        "n50": n50,
        "n90": n90,
        "gc_content": gc_content
    }
    
    return metrics


def export_graph_to_gfa(graph, output_file):
    """
    Export the overlap graph to GFA format for visualization.
    
    Args:
        graph: NetworkX directed graph
        output_file: Path to output GFA file
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write nodes (segments)
        for node in graph.nodes():
            sequence = graph.nodes[node].get('sequence', '')
            f.write("S\t{}\t{}\n".format(node, sequence))
        
        # Write edges (links)
        for source, target, data in graph.edges(data=True):
            overlap = data.get('overlap', 0)
            # In GFA, the orientation is + or -
            f.write("L\t{}\t+\t{}\t+\t{}M\n".format(source, target, overlap))
    
    logger.info("Exported overlap graph to GFA: {}".format(output_file))


def write_contigs_to_fasta(contigs, output_file):
    """
    Write contigs to a FASTA file.
    
    Args:
        contigs: List of contig sequences
        output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as f:
        for i, contig in enumerate(contigs):
            f.write(">contig_{} length={}\n".format(i+1, len(contig)))
            
            # Write sequence with 60 characters per line
            for j in range(0, len(contig), 60):
                f.write(contig[j:j+60] + "\n")
    
    logger.info("Wrote {} contigs to {}".format(len(contigs), output_file))


def align_contigs_to_reference(contigs, reference_file):
    """
    Align contigs to a reference genome and calculate alignment metrics.
    
    Args:
        contigs: List of contig sequences
        reference_file: Path to reference genome FASTA file
        
    Returns:
        Dictionary of alignment metrics
    """
    # Read reference genome
    reference_seq = ""
    try:
        for record in SeqIO.parse(reference_file, "fasta"):
            reference_seq = str(record.seq)
            break  # Just use the first sequence
    except Exception as e:
        logger.error("Error reading reference genome: {}".format(e))
        return {"error": str(e)}
    
    if not reference_seq:
        logger.error("No reference sequence found")
        return {"error": "No reference sequence found"}
    
    # Set up aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    # Align each contig and collect scores
    scores = []
    alignments = []
    
    for contig in contigs:
        if not contig:
            scores.append(0)
            alignments.append([])
            continue
        
        try:
            alignment = aligner.align(reference_seq, contig)[0]
            scores.append(alignment.score)
            alignments.append(alignment.path)
        except Exception as e:
            logger.error("Alignment error: {}".format(e))
            scores.append(0)
            alignments.append([])
    
    # Find best contig
    if scores:
        best_idx = np.argmax(scores)
        best_score = scores[best_idx]
        best_alignment = alignments[best_idx]
    else:
        best_idx = -1
        best_score = 0
        best_alignment = []
    
    return {
        "scores": scores,
        "best_score": best_score,
        "best_alignment": best_alignment,
        "best_contig_idx": best_idx
    }


def run_olc_assembly(fastq_file, min_overlap, output_dir, reference_file=None, max_reads=None):
    """
    Run the complete OLC assembly pipeline.
    
    Args:
        fastq_file: Path to input FASTQ file
        min_overlap: Minimum required overlap length
        output_dir: Directory for output files
        reference_file: Optional path to reference genome
        max_reads: Optional limit on number of reads to process
        
    Returns:
        Dictionary of results and metrics
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Read sequences
    sequences = read_fastq_sequences(fastq_file)
    
    # Limit sequences if specified (for testing)
    if max_reads and max_reads < len(sequences):
        seq_ids = list(sequences.keys())[:max_reads]
        sequences = {id: sequences[id] for id in seq_ids}
        logger.info("Limited to {} sequences for testing".format(len(sequences)))
    
    # 2. Compute pairwise overlaps
    overlaps = compute_pairwise_overlaps(sequences, min_overlap, max_reads)
    
    # 3. Build overlap graph
    graph = build_overlap_graph(sequences, overlaps)
    
    # 4. Find contig paths
    contig_paths = find_contigs_from_paths(graph)
    
    # 5. Generate consensus sequences
    contig_sequences = [
        generate_consensus_sequence(path, graph, sequences) 
        for path in contig_paths
    ]
    
    # 6. Calculate assembly metrics
    metrics = calculate_assembly_metrics(contig_sequences)
    
    # 7. Write outputs
    contigs_file = os.path.join(output_dir, "contigs_olc_min{}.fasta".format(min_overlap))
    write_contigs_to_fasta(contig_sequences, contigs_file)
    
    graph_file = os.path.join(output_dir, "assembly_graph_olc_min{}.gfa".format(min_overlap))
    export_graph_to_gfa(graph, graph_file)
    
    # 8. Align to reference if provided
    alignment_results = {}
    if reference_file:
        alignment_results = align_contigs_to_reference(contig_sequences, reference_file)
    
    # 9. Compile and return results
    results = {
        "sequences": {
            "total": len(sequences)
        },
        "overlaps": {
            "total": len(overlaps),
            "min_overlap": min_overlap
        },
        "graph": {
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges()
        },
        "assembly": metrics,
        "files": {
            "contigs": contigs_file,
            "graph": graph_file
        }
    }
    
    if alignment_results:
        results["alignment"] = alignment_results
    
    # Save results as JSON
    results_file = os.path.join(output_dir, "assembly_results_olc_min{}.json".format(min_overlap))
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=json_serializable)
    
    logger.info("Assembly complete. Results saved to {}".format(results_file))
    return results


def main():
    """Main entry point for the program."""
    parser = argparse.ArgumentParser(description="OLC Genome Assembler")
    
    # Required arguments
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output directory")
    
    # Optional arguments
    parser.add_argument("--min-overlap", type=int, default=20, 
                        help="Minimum overlap length (default: 20)")
    parser.add_argument("--reference", help="Reference genome FASTA file (optional)")
    parser.add_argument("--max-reads", type=int, help="Maximum number of reads to process (for testing)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Run the assembly pipeline
    try:
        run_olc_assembly(
            fastq_file=args.fastq,
            min_overlap=args.min_overlap,
            output_dir=args.output,
            reference_file=args.reference,
            max_reads=args.max_reads
        )
        logger.info("Assembly completed successfully")
    except Exception as e:
        logger.error("Assembly failed: {}".format(e))
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    main()