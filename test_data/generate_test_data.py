#!/usr/bin/env python3
"""
Script to generate random DNA sequences for testing primer design.

This script creates FASTA files with random DNA sequences of various lengths
and GC contents for testing the primer designer tool.
"""

import random
import argparse
from pathlib import Path


def generate_random_dna_sequence(length: int, gc_content: float = 0.5) -> str:
    """
    Generate a random DNA sequence with specified length and GC content.
    
    Args:
        length: Length of sequence to generate
        gc_content: Target GC content (0.0 to 1.0)
        
    Returns:
        Random DNA sequence string
    """
    # Calculate number of GC and AT bases
    gc_count = int(length * gc_content)
    at_count = length - gc_count
    
    # Create base pool
    bases = ['G'] * (gc_count // 2) + ['C'] * (gc_count // 2)
    bases += ['A'] * (at_count // 2) + ['T'] * (at_count // 2)
    
    # Handle odd numbers
    if len(bases) < length:
        remaining = length - len(bases)
        if gc_count % 2 == 1:
            bases.append('G')
            remaining -= 1
        if at_count % 2 == 1 and remaining > 0:
            bases.append('A')
    
    # Shuffle to randomize
    random.shuffle(bases)
    
    return ''.join(bases)


def generate_test_fasta(
    output_file: Path,
    num_sequences: int = 10,
    min_length: int = 200,
    max_length: int = 2000,
    gc_range: tuple = (0.3, 0.7)
) -> None:
    """
    Generate a FASTA file with random test sequences.
    
    Args:
        output_file: Output FASTA file path
        num_sequences: Number of sequences to generate
        min_length: Minimum sequence length
        max_length: Maximum sequence length
        gc_range: Range of GC content (min, max)
    """
    with open(output_file, 'w') as f:
        for i in range(1, num_sequences + 1):
            # Random length and GC content
            length = random.randint(min_length, max_length)
            gc_content = random.uniform(gc_range[0], gc_range[1])
            
            # Generate sequence
            sequence = generate_random_dna_sequence(length, gc_content)
            
            # Write FASTA record
            f.write(f">test_seq_{i:03d} length={length} gc={gc_content:.2f}\n")
            
            # Write sequence with line breaks every 80 characters
            for j in range(0, len(sequence), 80):
                f.write(sequence[j:j+80] + '\n')


def generate_specific_test_cases(output_dir: Path) -> None:
    """
    Generate specific test cases for different scenarios.
    
    Args:
        output_dir: Directory to save test files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Standard test case - mixed sequences
    print("Generating standard test case...")
    generate_test_fasta(
        output_dir / "standard_test.fasta",
        num_sequences=20,
        min_length=300,
        max_length=1500,
        gc_range=(0.4, 0.6)
    )
    
    # 2. GC-rich sequences
    print("Generating GC-rich test case...")
    generate_test_fasta(
        output_dir / "gc_rich_test.fasta",
        num_sequences=10,
        min_length=400,
        max_length=1200,
        gc_range=(0.6, 0.8)
    )
    
    # 3. AT-rich sequences
    print("Generating AT-rich test case...")
    generate_test_fasta(
        output_dir / "at_rich_test.fasta",
        num_sequences=10,
        min_length=400,
        max_length=1200,
        gc_range=(0.2, 0.4)
    )
    
    # 4. Large batch test
    print("Generating large batch test case...")
    generate_test_fasta(
        output_dir / "large_batch_test.fasta",
        num_sequences=100,
        min_length=200,
        max_length=3000,
        gc_range=(0.3, 0.7)
    )
    
    # 5. Small sequences (edge case)
    print("Generating small sequences test case...")
    generate_test_fasta(
        output_dir / "small_sequences_test.fasta",
        num_sequences=15,
        min_length=100,
        max_length=300,
        gc_range=(0.4, 0.6)
    )
    
    # 6. Very long sequences
    print("Generating long sequences test case...")
    generate_test_fasta(
        output_dir / "long_sequences_test.fasta",
        num_sequences=5,
        min_length=5000,
        max_length=10000,
        gc_range=(0.45, 0.55)
    )
    
    # 7. Single sequence test
    print("Generating single sequence test case...")
    generate_test_fasta(
        output_dir / "single_sequence_test.fasta",
        num_sequences=1,
        min_length=800,
        max_length=800,
        gc_range=(0.5, 0.5)
    )


def generate_real_like_sequences(output_file: Path, num_sequences: int = 10) -> None:
    """
    Generate more realistic DNA sequences with some structure.
    
    Args:
        output_file: Output FASTA file path
        num_sequences: Number of sequences to generate
    """
    # Common motifs and patterns found in real genes
    start_codons = ["ATG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    common_motifs = [
        "TATA", "CAAT", "GCCGCC", "CGCG", "AAAAAA", "TTTTTT",
        "GAATTC", "AAGCTT", "GGATCC", "CTGCAG"
    ]
    
    with open(output_file, 'w') as f:
        for i in range(1, num_sequences + 1):
            length = random.randint(500, 2500)
            gc_content = random.uniform(0.35, 0.65)
            
            # Start with start codon
            sequence = random.choice(start_codons)
            remaining_length = length - 3
            
            # Add some structured content
            while len(sequence) < length - 50:
                if random.random() < 0.1:  # 10% chance to add a motif
                    motif = random.choice(common_motifs)
                    sequence += motif
                    remaining_length -= len(motif)
                else:
                    # Add random bases
                    chunk_size = min(random.randint(10, 50), remaining_length)
                    chunk = generate_random_dna_sequence(chunk_size, gc_content)
                    sequence += chunk
                    remaining_length -= chunk_size
            
            # Fill remaining length
            if remaining_length > 3:
                sequence += generate_random_dna_sequence(remaining_length - 3, gc_content)
            
            # End with stop codon
            sequence += random.choice(stop_codons)
            
            # Ensure exact length
            sequence = sequence[:length]
            
            # Calculate actual GC content
            actual_gc = (sequence.count('G') + sequence.count('C')) / len(sequence)
            
            # Write FASTA record
            f.write(f">gene_like_seq_{i:03d} length={len(sequence)} gc={actual_gc:.3f} synthetic_gene\n")
            
            # Write sequence with line breaks
            for j in range(0, len(sequence), 80):
                f.write(sequence[j:j+80] + '\n')


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(description="Generate test FASTA files for primer design")
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("test_data"),
        help="Output directory for test files"
    )
    parser.add_argument(
        "--all", "-a",
        action="store_true",
        help="Generate all predefined test cases"
    )
    parser.add_argument(
        "--custom",
        action="store_true",
        help="Generate custom test file"
    )
    parser.add_argument(
        "--num-sequences", "-n",
        type=int,
        default=10,
        help="Number of sequences for custom generation"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=200,
        help="Minimum sequence length"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=2000,
        help="Maximum sequence length"
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        help="Output file for custom generation"
    )
    parser.add_argument(
        "--realistic",
        action="store_true",
        help="Generate realistic gene-like sequences"
    )
    
    args = parser.parse_args()
    
    # Set random seed for reproducible results
    random.seed(42)
    
    if args.all:
        print("Generating all predefined test cases...")
        generate_specific_test_cases(args.output_dir)
        
        # Also generate realistic sequences
        print("Generating realistic gene-like sequences...")
        generate_real_like_sequences(
            args.output_dir / "realistic_genes_test.fasta",
            num_sequences=15
        )
        
        print(f"All test files generated in: {args.output_dir}")
        
    elif args.custom:
        output_file = args.output_file or args.output_dir / "custom_test.fasta"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        if args.realistic:
            print(f"Generating {args.num_sequences} realistic sequences...")
            generate_real_like_sequences(output_file, args.num_sequences)
        else:
            print(f"Generating {args.num_sequences} random sequences...")
            generate_test_fasta(
                output_file,
                num_sequences=args.num_sequences,
                min_length=args.min_length,
                max_length=args.max_length
            )
        
        print(f"Custom test file generated: {output_file}")
        
    else:
        print("Use --all to generate all test cases or --custom for custom generation")
        print("Run with --help for more options")


if __name__ == "__main__":
    main()
