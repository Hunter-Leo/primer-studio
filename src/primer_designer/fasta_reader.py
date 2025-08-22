"""
FASTA file reading module using Biopython.

This module provides functionality to read and parse FASTA files
containing DNA sequences for primer design.
"""

import logging
from pathlib import Path
from typing import List, Iterator, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .models import SequenceInfo
from .utils import validate_dna_sequence


class FastaReader:
    """
    FASTA file reader for DNA sequences.
    
    This class provides methods to read FASTA files and convert them
    to SequenceInfo objects for primer design processing.
    """
    
    def __init__(self, min_sequence_length: int = 50):
        """
        Initialize FastaReader.
        
        Args:
            min_sequence_length: Minimum required sequence length for processing
        """
        self.min_sequence_length = min_sequence_length
        self.logger = logging.getLogger('primer_designer.fasta_reader')
    
    def read_fasta_file(self, file_path: Union[str, Path]) -> List[SequenceInfo]:
        """
        Read FASTA file and return list of SequenceInfo objects.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            List of SequenceInfo objects
            
        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file format is invalid or no valid sequences found
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")
        
        self.logger.info(f"Reading FASTA file: {file_path}")
        
        sequences = []
        invalid_sequences = 0
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                try:
                    sequence_info = self._process_record(record)
                    if sequence_info:
                        sequences.append(sequence_info)
                        self.logger.debug(f"Processed sequence: {sequence_info.sequence_id}")
                    else:
                        invalid_sequences += 1
                except Exception as e:
                    self.logger.warning(f"Skipping invalid sequence {record.id}: {e}")
                    invalid_sequences += 1
                    continue
        
        except Exception as e:
            raise ValueError(f"Error parsing FASTA file: {e}")
        
        if not sequences:
            raise ValueError("No valid DNA sequences found in FASTA file")
        
        self.logger.info(f"Successfully read {len(sequences)} sequences from {file_path}")
        if invalid_sequences > 0:
            self.logger.warning(f"Skipped {invalid_sequences} invalid sequences")
        
        return sequences
    
    def read_fasta_generator(self, file_path: Union[str, Path]) -> Iterator[SequenceInfo]:
        """
        Read FASTA file as generator for memory-efficient processing.
        
        Args:
            file_path: Path to FASTA file
            
        Yields:
            SequenceInfo objects one at a time
            
        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file format is invalid
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")
        
        self.logger.info(f"Reading FASTA file (generator mode): {file_path}")
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                try:
                    sequence_info = self._process_record(record)
                    if sequence_info:
                        self.logger.debug(f"Yielding sequence: {sequence_info.sequence_id}")
                        yield sequence_info
                except Exception as e:
                    self.logger.warning(f"Skipping invalid sequence {record.id}: {e}")
                    continue
        
        except Exception as e:
            raise ValueError(f"Error parsing FASTA file: {e}")
    
    def _process_record(self, record: SeqRecord) -> Union[SequenceInfo, None]:
        """
        Process a single SeqRecord and convert to SequenceInfo.
        
        Args:
            record: BioPython SeqRecord object
            
        Returns:
            SequenceInfo object or None if invalid
        """
        sequence_str = str(record.seq).upper()
        
        # Validate sequence length
        if len(sequence_str) < self.min_sequence_length:
            self.logger.debug(
                f"Sequence {record.id} too short: {len(sequence_str)} bp "
                f"(minimum: {self.min_sequence_length} bp)"
            )
            return None
        
        # Validate DNA sequence
        if not validate_dna_sequence(sequence_str, self.min_sequence_length):
            self.logger.debug(f"Invalid DNA sequence: {record.id}")
            return None
        
        # Create SequenceInfo object
        try:
            sequence_info = SequenceInfo(
                sequence_id=record.id,
                sequence=sequence_str,
                description=record.description if record.description != record.id else None,
                length=len(sequence_str)
            )
            return sequence_info
        
        except Exception as e:
            self.logger.debug(f"Error creating SequenceInfo for {record.id}: {e}")
            return None
    
    def validate_fasta_file(self, file_path: Union[str, Path]) -> dict:
        """
        Validate FASTA file and return statistics.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Dictionary with validation statistics
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            return {
                'valid': False,
                'error': f"File not found: {file_path}",
                'total_sequences': 0,
                'valid_sequences': 0,
                'invalid_sequences': 0
            }
        
        stats = {
            'valid': True,
            'error': None,
            'total_sequences': 0,
            'valid_sequences': 0,
            'invalid_sequences': 0,
            'sequence_lengths': [],
            'sequence_ids': []
        }
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                stats['total_sequences'] += 1
                stats['sequence_ids'].append(record.id)
                stats['sequence_lengths'].append(len(record.seq))
                
                sequence_info = self._process_record(record)
                if sequence_info:
                    stats['valid_sequences'] += 1
                else:
                    stats['invalid_sequences'] += 1
        
        except Exception as e:
            stats['valid'] = False
            stats['error'] = f"Error parsing FASTA file: {e}"
        
        return stats
    
    @staticmethod
    def write_sequences_to_fasta(sequences: List[SequenceInfo], output_path: Union[str, Path]) -> None:
        """
        Write SequenceInfo objects to FASTA file.
        
        Args:
            sequences: List of SequenceInfo objects
            output_path: Output FASTA file path
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        records = []
        for seq_info in sequences:
            record = SeqRecord(
                seq=seq_info.sequence,
                id=seq_info.sequence_id,
                description=seq_info.description or ""
            )
            records.append(record)
        
        SeqIO.write(records, output_path, "fasta")
    
    def get_sequence_count(self, file_path: Union[str, Path]) -> int:
        """
        Get the number of sequences in a FASTA file without loading all into memory.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Number of sequences in the file
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            return 0
        
        count = 0
        try:
            for _ in SeqIO.parse(file_path, "fasta"):
                count += 1
        except Exception:
            return 0
        
        return count
