"""
Primer Designer Package

A batch primer design tool for PCR applications that supports:
- Batch processing of multiple DNA sequences from FASTA files
- Automatic forward and reverse primer design
- Configurable primer design parameters
- Command-line interface and Python API
- Output in CSV or JSON formats
"""

from .config import PrimerDesignerConfig, STANDARD_CONFIG, GC_RICH_CONFIG, GC_POOR_CONFIG
from .models import Primer, SequenceInfo
from .primer import PrimerDesigner
from .batch_designer import BatchDesigner
from .fasta_reader import FastaReader

__version__ = "0.1.0"
__author__ = "Primer Studio Team"

__all__ = [
    "PrimerDesignerConfig",
    "STANDARD_CONFIG", 
    "GC_RICH_CONFIG",
    "GC_POOR_CONFIG",
    "Primer",
    "SequenceInfo", 
    "PrimerDesigner",
    "BatchDesigner",
    "FastaReader",
]
