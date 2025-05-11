#!/usr/bin/env python3
"""
GFA Exporter

This script converts De Bruijn graph data to GFA format for visualization in Bandage.
"""

import os
import argparse
import json
import logging
from typing import Set, List, Tuple, Dict

import networkx as nx
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def read_fastq_to_kmers(fastq_file: str, k: int) -> List[str]:
    """
    Read sequences from a FASTQ file and convert to k-mers.
    
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
        logger.error(f"Error reading {fastq_file}: {e}")
        
    logger.info(f"Extracted {len(kmers)} k-mers from {fastq_file}")
    return kmers


def build_debruijn_graph(kmers: List[str]) -> Tuple[Set[str], List[Tuple[str, str]]]:
    """
    Construct a De Bruijn graph from k-mers.
    
    Args:
        kmers: List of k-mers
        
    Returns:
        Tuple containing:
        - Set of nodes (k-1 mers)
        - List of edges (connections between k-1 mers)
    """
    # Initialize data structures
    nodes = set()
    edges = []
    
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
    
    logger.info(f"Built De Bruijn graph with {len(nodes)} nodes and {len(edges)} edges")
    return nodes, edges


def export_to_gfa(nodes: Set[str], edges: List[Tuple[str, str]], output_file: str) -> None:
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
            f.write(f"S\t{node_id}\t{node}\n")
        
        # Create a map from node names to IDs
        node_to_id = {node: str(i+1) for i, node in enumerate(nodes)}
        
        # Write edges (links)
        for source, target in edges:
            source_id = node_to_id[source]
            target_id = node_to_id[target]
            # In GFA, the orientation is + or -
            f.write(f"L\t{source_id}\t+\t{target_id}\t+\t{len(source)-1}M\n")
    
    logger.info(f"Exported graph with {len(nodes)} nodes and {len(edges)} edges to {output_file}")


def main():
    """Main entry point for the program."""
    parser = argparse.ArgumentParser(description="Export De Bruijn graph to GFA format")
    
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output GFA file")
    parser.add_argument("--k", type=int, default=40, help="k-mer size (default: 40)")
    parser.add_argument("--min-freq", type=int, default=1, 
                        help="Minimum k-mer frequency (default: 1)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    try:
        # Read k-mers from FASTQ
        kmers = read_fastq_to_kmers(args.fastq, args.k)
        
        # Filter k-mers by frequency if requested
        if args.min_freq > 1:
            from collections import Counter
            kmer_counts = Counter(kmers)
            kmers = [kmer for kmer in kmers if kmer_counts[kmer] >= args.min_freq]
            logger.info(f"Filtered to {len(set(kmers))} unique k-mers with frequency >= {args.min_freq}")
        
        # Build graph
        nodes, edges = build_debruijn_graph(kmers)
        
        # Export to GFA
        export_to_gfa(nodes, edges, args.output)
        
        logger.info(f"GFA export completed successfully: {args.output}")
        
    except Exception as e:
        logger.error(f"GFA export failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    main()