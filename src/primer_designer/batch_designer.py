"""
Batch primer design module for processing multiple sequences.

This module provides the BatchDesigner class for efficient batch processing
of multiple DNA sequences with optional parallel processing.
"""

import json
import csv
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Union, Iterator, Dict, Any

from .config import PrimerDesignerConfig, STANDARD_CONFIG
from .models import Primer, SequenceInfo
from .primer import PrimerDesigner
from .fasta_reader import FastaReader


class BatchDesigner:
    """
    Batch primer designer for processing multiple DNA sequences.
    
    This class provides functionality to process multiple sequences
    efficiently with optional parallel processing and various output formats.
    """
    
    def __init__(
        self,
        config: Optional[PrimerDesignerConfig] = None,
        max_workers: Optional[int] = None,
        use_parallel: bool = True
    ):
        """
        Initialize BatchDesigner.
        
        Args:
            config: PrimerDesignerConfig object. If None, uses STANDARD_CONFIG
            max_workers: Maximum number of worker threads for parallel processing
            use_parallel: Whether to use parallel processing
        """
        self.config = config or STANDARD_CONFIG
        self.max_workers = max_workers
        self.use_parallel = use_parallel
        self.logger = logging.getLogger('primer_designer.batch')
        
        # Initialize primer designer
        self.primer_designer = PrimerDesigner(self.config)
        
        # Initialize FASTA reader
        self.fasta_reader = FastaReader()
        
        self.logger.info(f"BatchDesigner initialized (parallel: {use_parallel})")
    
    def design_primers_from_fasta(
        self,
        fasta_file: Union[str, Path],
        output_file: Optional[Union[str, Path]] = None,
        output_format: str = "csv"
    ) -> List[Primer]:
        """
        Design primers for all sequences in a FASTA file.
        
        Args:
            fasta_file: Path to input FASTA file
            output_file: Optional path to output file
            output_format: Output format ('csv' or 'json')
            
        Returns:
            List of Primer objects
            
        Raises:
            FileNotFoundError: If FASTA file does not exist
            ValueError: If invalid output format specified
        """
        self.logger.info(f"Starting batch primer design from: {fasta_file}")
        
        # Validate output format
        if output_format not in ['csv', 'json']:
            raise ValueError("Output format must be 'csv' or 'json'")
        
        # Read sequences from FASTA file
        sequences = self.fasta_reader.read_fasta_file(fasta_file)
        self.logger.info(f"Loaded {len(sequences)} sequences from FASTA file")
        
        # Design primers
        primers = self.design_primers_batch(sequences)
        
        # Save results if output file specified
        if output_file:
            self.save_results(primers, output_file, output_format)
        
        return primers
    
    def design_primers_from_fasta_generator(
        self,
        fasta_file: Union[str, Path],
        output_file: Union[str, Path],
        output_format: str = "csv"
    ) -> Iterator[Primer]:
        """
        Design primers using generator for memory-efficient processing.
        
        Args:
            fasta_file: Path to input FASTA file
            output_file: Path to output file
            output_format: Output format ('csv' or 'json')
            
        Yields:
            Primer objects as they are processed
        """
        self.logger.info(f"Starting memory-efficient batch processing: {fasta_file}")
        
        # Validate output format
        if output_format not in ['csv', 'json']:
            raise ValueError("Output format must be 'csv' or 'json'")
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        primers_processed = 0
        primers_successful = 0
        
        # Initialize output file
        if output_format == "csv":
            csv_file = open(output_path, 'w', newline='', encoding='utf-8')
            csv_writer = None
        else:  # json
            json_results = []
        
        try:
            # Process sequences one by one
            for sequence_info in self.fasta_reader.read_fasta_generator(fasta_file):
                primers_processed += 1
                
                # Design primer
                primer = self.primer_designer.design_primer(sequence_info)
                
                if primer:
                    primers_successful += 1
                    
                    # Write to output file immediately
                    if output_format == "csv":
                        if csv_writer is None:
                            # Initialize CSV writer with headers
                            fieldnames = list(primer.to_dict().keys())
                            csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                            csv_writer.writeheader()
                        
                        csv_writer.writerow(primer.to_dict())
                        csv_file.flush()  # Ensure data is written
                    else:  # json
                        json_results.append(primer.to_dict())
                    
                    yield primer
                
                if primers_processed % 100 == 0:
                    self.logger.info(f"Processed {primers_processed} sequences, {primers_successful} successful")
        
        finally:
            # Clean up
            if output_format == "csv" and csv_file:
                csv_file.close()
            elif output_format == "json":
                with open(output_path, 'w', encoding='utf-8') as f:
                    json.dump(json_results, f, indent=2)
        
        self.logger.info(f"Batch processing completed: {primers_successful}/{primers_processed} successful")
    
    def design_primers_batch(self, sequences: List[SequenceInfo]) -> List[Primer]:
        """
        Design primers for a batch of sequences.
        
        Args:
            sequences: List of SequenceInfo objects
            
        Returns:
            List of Primer objects (may be shorter than input if some designs failed)
        """
        self.logger.info(f"Designing primers for {len(sequences)} sequences")
        
        primers = []
        
        if self.use_parallel and len(sequences) > 1:
            primers = self._design_primers_parallel(sequences)
        else:
            primers = self._design_primers_sequential(sequences)
        
        success_rate = len(primers) / len(sequences) * 100 if sequences else 0
        self.logger.info(f"Batch design completed: {len(primers)}/{len(sequences)} successful ({success_rate:.1f}%)")
        
        return primers
    
    def _design_primers_sequential(self, sequences: List[SequenceInfo]) -> List[Primer]:
        """
        Design primers sequentially (single-threaded).
        
        Args:
            sequences: List of SequenceInfo objects
            
        Returns:
            List of Primer objects
        """
        primers = []
        
        for i, sequence_info in enumerate(sequences, 1):
            self.logger.debug(f"Processing sequence {i}/{len(sequences)}: {sequence_info.sequence_id}")
            
            try:
                primer = self.primer_designer.design_primer(sequence_info)
                if primer:
                    primers.append(primer)
            except Exception as e:
                self.logger.warning(f"Failed to design primer for {sequence_info.sequence_id}: {e}")
            
            # Log progress for large batches
            if i % 50 == 0:
                self.logger.info(f"Progress: {i}/{len(sequences)} sequences processed")
        
        return primers
    
    def _design_primers_parallel(self, sequences: List[SequenceInfo]) -> List[Primer]:
        """
        Design primers in parallel using ThreadPoolExecutor.
        
        Args:
            sequences: List of SequenceInfo objects
            
        Returns:
            List of Primer objects
        """
        primers = []
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_sequence = {
                executor.submit(self._design_single_primer, seq): seq
                for seq in sequences
            }
            
            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_sequence):
                completed += 1
                sequence_info = future_to_sequence[future]
                
                try:
                    primer = future.result()
                    if primer:
                        primers.append(primer)
                except Exception as e:
                    self.logger.warning(f"Failed to design primer for {sequence_info.sequence_id}: {e}")
                
                # Log progress
                if completed % 50 == 0:
                    self.logger.info(f"Progress: {completed}/{len(sequences)} sequences processed")
        
        return primers
    
    def _design_single_primer(self, sequence_info: SequenceInfo) -> Optional[Primer]:
        """
        Design primer for a single sequence (thread-safe wrapper).
        
        Args:
            sequence_info: SequenceInfo object
            
        Returns:
            Primer object or None if design failed
        """
        try:
            # Create a new PrimerDesigner instance for thread safety
            designer = PrimerDesigner(self.config)
            return designer.design_primer(sequence_info)
        except Exception as e:
            self.logger.debug(f"Primer design failed for {sequence_info.sequence_id}: {e}")
            return None
    
    def save_results(
        self,
        primers: List[Primer],
        output_file: Union[str, Path],
        output_format: str = "csv"
    ) -> None:
        """
        Save primer design results to file.
        
        Args:
            primers: List of Primer objects
            output_file: Output file path
            output_format: Output format ('csv' or 'json')
            
        Raises:
            ValueError: If invalid output format specified
        """
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.logger.info(f"Saving {len(primers)} results to: {output_path}")
        
        if output_format == "csv":
            self._save_csv(primers, output_path)
        elif output_format == "json":
            self._save_json(primers, output_path)
        else:
            raise ValueError("Output format must be 'csv' or 'json'")
        
        self.logger.info(f"Results saved successfully to: {output_path}")
    
    def _save_csv(self, primers: List[Primer], output_path: Path) -> None:
        """Save results in CSV format."""
        if not primers:
            # Create empty CSV with headers
            with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([
                    'sequence_id', 'forward_primer', 'reverse_primer',
                    'tm_forward', 'tm_reverse', 'gc_forward', 'gc_reverse',
                    'length_forward', 'length_reverse', 'product_size',
                    'forward_start', 'reverse_start', 'penalty_forward', 'penalty_reverse'
                ])
            return
        
        # Get fieldnames from first primer
        fieldnames = list(primers[0].to_dict().keys())
        
        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for primer in primers:
                writer.writerow(primer.to_dict())
    
    def _save_json(self, primers: List[Primer], output_path: Path) -> None:
        """Save results in JSON format."""
        results = {
            'metadata': {
                'total_primers': len(primers),
                'config_summary': self.primer_designer.get_config_summary()
            },
            'primers': [primer.to_dict() for primer in primers]
        }
        
        with open(output_path, 'w', encoding='utf-8') as jsonfile:
            json.dump(results, jsonfile, indent=2)
    
    def get_batch_statistics(self, primers: List[Primer]) -> Dict[str, Any]:
        """
        Calculate statistics for a batch of primers.
        
        Args:
            primers: List of Primer objects
            
        Returns:
            Dictionary with batch statistics
        """
        if not primers:
            return {'total_primers': 0}
        
        # Calculate statistics
        tm_forward_values = [p.tm_forward for p in primers]
        tm_reverse_values = [p.tm_reverse for p in primers]
        gc_forward_values = [p.gc_forward for p in primers]
        gc_reverse_values = [p.gc_reverse for p in primers]
        product_sizes = [p.product_size for p in primers]
        
        stats = {
            'total_primers': len(primers),
            'tm_statistics': {
                'forward': {
                    'min': min(tm_forward_values),
                    'max': max(tm_forward_values),
                    'mean': sum(tm_forward_values) / len(tm_forward_values)
                },
                'reverse': {
                    'min': min(tm_reverse_values),
                    'max': max(tm_reverse_values),
                    'mean': sum(tm_reverse_values) / len(tm_reverse_values)
                }
            },
            'gc_statistics': {
                'forward': {
                    'min': min(gc_forward_values),
                    'max': max(gc_forward_values),
                    'mean': sum(gc_forward_values) / len(gc_forward_values)
                },
                'reverse': {
                    'min': min(gc_reverse_values),
                    'max': max(gc_reverse_values),
                    'mean': sum(gc_reverse_values) / len(gc_reverse_values)
                }
            },
            'product_size_statistics': {
                'min': min(product_sizes),
                'max': max(product_sizes),
                'mean': sum(product_sizes) / len(product_sizes)
            }
        }
        
        return stats
