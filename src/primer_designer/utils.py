"""
Utility functions for primer design calculations and logging setup.

This module provides helper functions using Biopython for sequence operations
and logging configuration.
"""

import logging
import sys
from typing import Optional
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.MeltingTemp import Tm_NN


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage using Biopython.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as percentage (0-100)
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    seq_obj = Seq(sequence.upper().strip())
    return round(gc_fraction(seq_obj) * 100, 2)


def calculate_tm_nn(sequence: str, salt_conc: float = 50.0, primer_conc: float = 50.0) -> float:
    """
    Calculate melting temperature using Biopython's nearest-neighbor method.
    
    Args:
        sequence: Primer sequence
        salt_conc: Salt concentration in mM (default: 50.0)
        primer_conc: Primer concentration in nM (default: 50.0)
        
    Returns:
        Melting temperature in °C
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    seq_obj = Seq(sequence.upper().strip())
    
    # Convert concentrations to M for Biopython
    salt_conc_m = salt_conc / 1000.0  # mM to M
    primer_conc_m = primer_conc / 1000000000.0  # nM to M
    
    tm = Tm_NN(seq_obj, Na=salt_conc_m, dnac1=primer_conc_m, dnac2=primer_conc_m)
    return round(tm, 2)


def reverse_complement(sequence: str) -> str:
    """
    Generate reverse complement using Biopython.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Reverse complement sequence
    """
    if not sequence:
        return ""
    
    seq_obj = Seq(sequence.upper().strip())
    return str(seq_obj.reverse_complement())


def validate_dna_sequence(sequence: str, min_length: int = 50) -> bool:
    """
    Validate that a sequence is a valid DNA sequence.
    
    Args:
        sequence: DNA sequence to validate
        min_length: Minimum required length
        
    Returns:
        True if sequence is valid, False otherwise
    """
    if not sequence or len(sequence) < min_length:
        return False
    
    try:
        seq_obj = Seq(sequence.upper().strip())
        # Try to calculate GC content - will raise exception if invalid bases
        gc_fraction(seq_obj)
        return True
    except Exception:
        return False


def get_molecular_weight(sequence: str) -> float:
    """
    Calculate molecular weight using Biopython.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Molecular weight in g/mol
    """
    if not sequence:
        return 0.0
    
    seq_obj = Seq(sequence.upper().strip())
    return round(molecular_weight(seq_obj, seq_type='DNA'), 2)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_string: Optional[str] = None
) -> logging.Logger:
    """
    Set up logging configuration for the primer designer application.
    
    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional path to log file
        format_string: Optional custom format string
        
    Returns:
        Configured logger instance
    """
    # Convert string level to logging constant
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    
    # Default format string
    if format_string is None:
        format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Create formatter
    formatter = logging.Formatter(format_string)
    
    # Get or create logger
    logger = logging.getLogger('primer_designer')
    logger.setLevel(numeric_level)
    
    # Remove existing handlers to avoid duplicates
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(numeric_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (optional)
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(numeric_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def format_sequence_for_display(sequence: str, line_length: int = 80) -> str:
    """
    Format a DNA sequence for display with line breaks.
    
    Args:
        sequence: DNA sequence to format
        line_length: Number of characters per line
        
    Returns:
        Formatted sequence string with line breaks
    """
    if not sequence:
        return ""
    
    lines = []
    for i in range(0, len(sequence), line_length):
        lines.append(sequence[i:i + line_length])
    
    return '\n'.join(lines)


def get_sequence_stats(sequence: str) -> dict:
    """
    Calculate basic statistics for a DNA sequence using Biopython.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Dictionary with sequence statistics
    """
    if not sequence:
        return {}
    
    seq_obj = Seq(sequence.upper().strip())
    length = len(seq_obj)
    
    # Count bases
    sequence_str = str(seq_obj)
    base_counts = {
        'A': sequence_str.count('A'),
        'T': sequence_str.count('T'),
        'G': sequence_str.count('G'),
        'C': sequence_str.count('C'),
        'N': sequence_str.count('N')
    }
    
    # Calculate percentages
    base_percentages = {
        base: round((count / length) * 100, 2) if length > 0 else 0
        for base, count in base_counts.items()
    }
    
    gc_content = gc_fraction(seq_obj) * 100
    
    return {
        'length': length,
        'base_counts': base_counts,
        'base_percentages': base_percentages,
        'gc_content': round(gc_content, 2),
        'at_content': round(100 - gc_content, 2),
        'molecular_weight': get_molecular_weight(str(seq_obj))
    }
