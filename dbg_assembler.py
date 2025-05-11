#!/usr/bin/env python3
"""
De Bruijn Graph (DBG) Genome Assembler

This module implements a De Bruijn Graph approach to genome assembly,
taking FASTQ files as input and generating contigs as FASTA output.
"""

import os
import json
import argparse
import logging
from collections import defaultdict

import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio import Align

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


def read_fastq_to_kmers(fastq_file, k):
    """
    Read sequences from a single FASTQ file and convert to k-mers.
    
    Args:
        fastq_file: Path to FASTQ file
        k: Size of k-mers to generate
        
    Returns:
        List of k-mers extracted from all sequences
    """
    kmers = []
    try:
        for record in SeqIO.parse(fastq_file, "fastq"):
            sequence = str(record.seq)
            # Extract all k-mers from this sequence
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if len(kmer) == k:  # Ensure full-length k-mers
                    kmers.append(kmer)
    except Exception as e:
        logger.error("Error reading {}: {}".format(fastq_file, e))
        
    logger.info("Extracted {} k-mers from {}".format(len(kmers), fastq_file))
    return kmers


def filter_kmers_by_frequency(kmers, min_frequency):
    """
    Filter k-mers based on their frequency in the dataset.
    
    Args:
        kmers: List of k-mers
        min_frequency: Minimum frequency threshold to keep a k-mer
        
    Returns:
        Filtered list of k-mers that meet the threshold
    """
    # Count k-mer frequencies
    kmer_counts = defaultdict(int)
    for kmer in kmers:
        kmer_counts[kmer] += 1
    
    # Filter k-mers below threshold
    filtered_kmers = [kmer for kmer in kmers if kmer_counts[kmer] >= min_frequency]
    
    # Report statistics
    total_unique = len(kmer_counts)
    kept_unique = len({k for k in filtered_kmers})
    logger.info("Filtering: {} unique k-mers before, {} after (threshold={})".format(
        total_unique, kept_unique, min_frequency))
    
    return filtered_kmers


def build_debruijn_graph(kmers):
    """
    Construct a De Bruijn graph from k-mers.
    
    Args:
        kmers: List of k-mers
        
    Returns:
        Tuple containing:
        - Set of nodes (k-1 mers)
        - List of edges (connections between k-1 mers)
        - String representing a start node
    """
    # Initialize data structures
    nodes = set()
    edges = []
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    # Build the graph
    for kmer in kmers:
        # Create nodes from (k-1)-mers
        prefix = kmer[:-1]  # First k-1 characters
        suffix = kmer[1:]   # Last k-1 characters
        
        # Add nodes to set
        nodes.add(prefix)
        nodes.add(suffix)
        
        # Add edge
        edges.append((prefix, suffix))
        
        # Track node degrees for determining start nodes
        out_degree[prefix] += 1
        in_degree[suffix] += 1
    
    # Find potential start nodes (more outgoing than incoming edges)
    start_candidates = []
    for node in nodes:
        degree_diff = out_degree[node] - in_degree[node]
        if degree_diff > 0:
            start_candidates.append((degree_diff, node))
    
    # Select start node
    start_node = None
    if start_candidates:
        # Use node with largest out-in degree difference
        start_candidates.sort(reverse=True)
        start_node = start_candidates[0][1]
    else:
        # If no obvious start, pick any node with outgoing edges
        for node in nodes:
            if out_degree[node] > 0:
                start_node = node
                break
    
    logger.info("Built De Bruijn graph with {} nodes and {} edges".format(len(nodes), len(edges)))
    
    return nodes, edges, start_node


def create_adjacency_list(edges):
    """
    Convert list of edges to an adjacency list representation.
    
    Args:
        edges: List of (source, target) edges
        
    Returns:
        Dictionary mapping source nodes to lists of target nodes
    """
    adjacency = defaultdict(list)
    for source, target in edges:
        adjacency[source].append(target)
    return dict(adjacency)


def find_start_node(edges):
    """
    Find a suitable start node based on degree differences.
    
    Args:
        edges: List of edges in the graph
        
    Returns:
        A node suitable for starting graph traversal
    """
    # Calculate in and out degrees
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    for source, target in edges:
        out_degree[source] += 1
        in_degree[target] += 1
    
    # Find nodes with out_degree > in_degree (potential starts)
    start_candidates = []
    for node in set(out_degree.keys()) | set(in_degree.keys()):
        degree_diff = out_degree[node] - in_degree[node]
        if degree_diff > 0:
            start_candidates.append((degree_diff, node))
    
    # Select the best start node
    if start_candidates:
        start_candidates.sort(reverse=True)
        return start_candidates[0][1]
    
    # If no obvious start, pick a node with outgoing edges
    for node, degree in out_degree.items():
        if degree > 0:
            return node
    
    return None


def traverse_graph(adjacency_list, start_node):
    """
    Traverse the graph to find a path (pseudo-Eulerian path).
    
    Args:
        adjacency_list: Dictionary mapping nodes to their neighbors
        start_node: Starting node for traversal
        
    Returns:
        List of nodes in the traversal path
    """
    if not start_node or start_node not in adjacency_list:
        logger.warning("No valid start node for traversal")
        return []
    
    # Make a copy of the adjacency list to modify during traversal
    adj_list = {k: v.copy() for k, v in adjacency_list.items() if v}
    
    # Initialize the path with the start node
    path = [start_node]
    current = start_node
    
    # Continue traversal until we can't go further
    while current in adj_list and adj_list[current]:
        # Get the next node and remove the edge
        next_node = adj_list[current].pop()
        if not adj_list[current]:  # If no more outgoing edges
            del adj_list[current]
        
        # Add the next node to the path
        path.append(next_node)
        current = next_node
    
    logger.info("Graph traversal complete: path length = {}".format(len(path)))
    return path


def generate_contig_from_path(path):
    """
    Generate a contig sequence from a path of (k-1)-mers.
    
    Args:
        path: List of (k-1)-mers representing a path in the graph
        
    Returns:
        Assembled contig sequence
    """
    if not path:
        return ""
    
    # Start with the first node
    contig = path[0]
    
    # For each subsequent node, add only the last character
    for node in path[1:]:
        contig += node[-1]
    
    return contig


def process_connected_components(graph):
    """
    Process each weakly connected component to generate contigs.
    
    Args:
        graph: NetworkX directed graph
        
    Returns:
        List of contigs (one per connected component)
    """
    contigs = []
    
    # Find weakly connected components
    components = list(nx.weakly_connected_components(graph))
    logger.info("Found {} weakly connected components".format(len(components)))
    
    for i, component in enumerate(components):
        # Extract the subgraph for this component
        subgraph = graph.subgraph(component)
        edges = list(subgraph.edges())
        
        # Skip empty components
        if not edges:
            continue
        
        # Create adjacency list and find start node
        adj_list = create_adjacency_list(edges)
        start = find_start_node(edges)
        
        # If we have a valid start node, traverse and generate contig
        if start:
            path = traverse_graph(adj_list, start)
            contig = generate_contig_from_path(path)
            if contig:
                contigs.append(contig)
                logger.debug("Component {}: Generated contig of length {}".format(i+1, len(contig)))
    
    return contigs


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


def export_graph_to_gfa(nodes, edges, output_file):
    """
    Export the De Bruijn graph to GFA format for visualization.
    
    Args:
        nodes: Set of nodes in the graph
        edges: List of edges in the graph
        output_file: Path to output GFA file
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write nodes (segments)
        for i, node in enumerate(nodes):
            node_id = str(i+1)
            f.write("S\t{}\t{}\n".format(node_id, node))
        
        # Create a map from node names to IDs
        node_to_id = {node: str(i+1) for i, node in enumerate(nodes)}
        
        # Write edges (links)
        for source, target in edges:
            source_id = node_to_id[source]
            target_id = node_to_id[target]
            # In GFA, the orientation is + or -
            f.write("L\t{}\t+\t{}\t+\t{}M\n".format(source_id, target_id, len(source)-1))
    
    logger.info("Exported graph with {} nodes and {} edges to {}".format(
        len(nodes), len(edges), output_file))


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
    aligner = Align.PairwiseAligner()
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
            alignments.append(alignment.aligned)
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


def run_dbg_assembly(fastq_file, k, min_frequency, output_dir, reference_file=None):
    """
    Run the complete DBG assembly pipeline.
    
    Args:
        fastq_file: Path to input FASTQ file
        k: k-mer size
        min_frequency: Minimum frequency threshold for k-mers
        output_dir: Directory for output files
        reference_file: Optional path to reference genome
        
    Returns:
        Dictionary of results and metrics
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Read sequences and extract k-mers
    kmers = read_fastq_to_kmers(fastq_file, k)
    
    # 2. Filter k-mers by frequency
    if min_frequency > 1:
        kmers = filter_kmers_by_frequency(kmers, min_frequency)
    
    # 3. Build De Bruijn graph
    nodes, edges, start_node = build_debruijn_graph(kmers)
    
    # 4. Create NetworkX graph for component analysis
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    
    # 5. Process connected components to generate contigs
    contigs = process_connected_components(graph)
    
    # 6. Calculate assembly metrics
    metrics = calculate_assembly_metrics(contigs)
    
    # 7. Write outputs
    contigs_file = os.path.join(output_dir, "contigs_k{}.fasta".format(k))
    write_contigs_to_fasta(contigs, contigs_file)
    
    graph_file = os.path.join(output_dir, "assembly_graph_k{}.gfa".format(k))
    export_graph_to_gfa(nodes, edges, graph_file)
    
    # 8. Align to reference if provided
    alignment_results = {}
    if reference_file:
        alignment_results = align_contigs_to_reference(contigs, reference_file)
    
    # 9. Compile and return results
    results = {
        "kmers": {
            "total": len(kmers),
            "unique": len(set(kmers))
        },
        "graph": {
            "nodes": len(nodes),
            "edges": len(edges)
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
    results_file = os.path.join(output_dir, "assembly_results_k{}.json".format(k))
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=json_serializable)
    
    logger.info("Assembly complete. Results saved to {}".format(results_file))
    return results


def main():
    """Main entry point for the program."""
    parser = argparse.ArgumentParser(description="De Bruijn Graph Genome Assembler")
    
    # Required arguments
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output directory")
    
    # Optional arguments
    parser.add_argument("--k", type=int, default=31, help="k-mer size (default: 31)")
    parser.add_argument("--min-freq", type=int, default=1, 
                        help="Minimum k-mer frequency (default: 1)")
    parser.add_argument("--reference", help="Reference genome FASTA file (optional)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Run the assembly pipeline
    try:
        run_dbg_assembly(
            fastq_file=args.fastq,
            k=args.k,
            min_frequency=args.min_freq,
            output_dir=args.output,
            reference_file=args.reference
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