"""
Core primer design module using primer3-py.

This module provides the PrimerDesigner class that encapsulates
primer3-py functionality for designing PCR primers.
"""

import logging
from typing import Optional, Dict, Any

import primer3

from .config import PrimerDesignerConfig, STANDARD_CONFIG
from .models import Primer, SequenceInfo
from .utils import calculate_gc_content


class PrimerDesigner:
    """
    Primer designer class that wraps primer3-py functionality.
    
    This class provides an object-oriented interface to primer3-py
    for designing PCR primers with configurable parameters.
    """
    
    def __init__(self, config: Optional[PrimerDesignerConfig] = None):
        """
        Initialize PrimerDesigner with configuration.
        
        Args:
            config: PrimerDesignerConfig object. If None, uses STANDARD_CONFIG
        """
        self.config = config or STANDARD_CONFIG
        self.logger = logging.getLogger('primer_designer.primer')
        
        # Validate configuration
        self._validate_config()
        
        self.logger.info("PrimerDesigner initialized with configuration")
        self.logger.debug(f"Configuration: {self.config.model_dump()}")
    
    def _validate_config(self) -> None:
        """Validate the primer design configuration."""
        try:
            # Test conversion to primer3 dictionary
            primer3_dict = self.config.to_primer3_dict()
            self.logger.debug("Configuration validation successful")
        except Exception as e:
            raise ValueError(f"Invalid primer design configuration: {e}")
    
    def _split_primer3_args(self, primer3_input: Dict[str, Any]) -> tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Split primer3 input into seq_args and global_args.
        
        Args:
            primer3_input: Combined primer3 input dictionary
            
        Returns:
            Tuple of (seq_args, global_args)
        """
        # Sequence-specific arguments
        seq_args = {
            'SEQUENCE_ID': primer3_input.get('SEQUENCE_ID', ''),
            'SEQUENCE_TEMPLATE': primer3_input.get('SEQUENCE_TEMPLATE', ''),
        }
        
        # Global arguments (everything else)
        global_args = {k: v for k, v in primer3_input.items() 
                      if k not in ['SEQUENCE_ID', 'SEQUENCE_TEMPLATE']}
        
        return seq_args, global_args
    
    def design_primer(self, sequence_info: SequenceInfo) -> Optional[Primer]:
        """
        Design primers for a single DNA sequence.
        
        Args:
            sequence_info: SequenceInfo object containing sequence data
            
        Returns:
            Primer object with design results, or None if design failed
            
        Raises:
            ValueError: If sequence is invalid or too short
        """
        self.logger.info(f"Designing primers for sequence: {sequence_info.sequence_id}")
        
        # Validate sequence
        if len(sequence_info.sequence) < 100:
            raise ValueError(f"Sequence too short for primer design: {len(sequence_info.sequence)} bp")
        
        # Prepare primer3 input
        primer3_input = self._prepare_primer3_input(sequence_info)
        
        try:
            # Run primer3 design using the new API
            seq_args, global_args = self._split_primer3_args(primer3_input)
            result = primer3.design_primers(seq_args, global_args)
            
            # Check if primers were found
            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) == 0:
                self.logger.warning(f"No primers found for sequence: {sequence_info.sequence_id}")
                self._log_primer3_errors(result)
                return None
            
            # Extract best primer pair
            primer = self._extract_primer_from_result(sequence_info, result)
            
            self.logger.info(f"Successfully designed primers for: {sequence_info.sequence_id}")
            return primer
        
        except Exception as e:
            self.logger.error(f"Primer design failed for {sequence_info.sequence_id}: {e}")
            return None
    
    def _prepare_primer3_input(self, sequence_info: SequenceInfo) -> Dict[str, Any]:
        """
        Prepare input dictionary for primer3.
        
        Args:
            sequence_info: SequenceInfo object
            
        Returns:
            Dictionary with primer3 input parameters
        """
        # Start with configuration parameters
        primer3_input = self.config.to_primer3_dict()
        
        # Add sequence-specific parameters
        primer3_input.update({
            'SEQUENCE_ID': sequence_info.sequence_id,
            'SEQUENCE_TEMPLATE': sequence_info.sequence,
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_NUM_RETURN': 5,  # Return top 5 primer pairs
        })
        
        return primer3_input
    
    def _extract_primer_from_result(self, sequence_info: SequenceInfo, result: Dict[str, Any]) -> Primer:
        """
        Extract primer information from primer3 result.
        
        Args:
            sequence_info: Original sequence information
            result: primer3 design result dictionary
            
        Returns:
            Primer object with extracted information
        """
        # Extract best primer pair (index 0)
        left_seq = result['PRIMER_LEFT_0_SEQUENCE']
        right_seq = result['PRIMER_RIGHT_0_SEQUENCE']
        
        # Extract positions (primer3 uses 0-based indexing)
        left_pos = result['PRIMER_LEFT_0'][0]  # Start position
        right_pos = result['PRIMER_RIGHT_0'][0]  # Start position
        
        # Extract melting temperatures
        left_tm = result['PRIMER_LEFT_0_TM']
        right_tm = result['PRIMER_RIGHT_0_TM']
        
        # Calculate GC content
        left_gc = calculate_gc_content(left_seq)
        right_gc = calculate_gc_content(right_seq)
        
        # Extract penalty scores
        left_penalty = result.get('PRIMER_LEFT_0_PENALTY', 0.0)
        right_penalty = result.get('PRIMER_RIGHT_0_PENALTY', 0.0)
        
        # Calculate product size
        product_size = result['PRIMER_PAIR_0_PRODUCT_SIZE']
        
        return Primer(
            sequence_id=sequence_info.sequence_id,
            forward_primer=left_seq,
            reverse_primer=right_seq,
            tm_forward=round(left_tm, 2),
            tm_reverse=round(right_tm, 2),
            gc_forward=left_gc,
            gc_reverse=right_gc,
            length_forward=len(left_seq),
            length_reverse=len(right_seq),
            product_size=product_size,
            forward_start=left_pos,
            reverse_start=right_pos,
            penalty_forward=left_penalty,
            penalty_reverse=right_penalty
        )
    
    def _log_primer3_errors(self, result: Dict[str, Any]) -> None:
        """
        Log primer3 error messages for debugging.
        
        Args:
            result: primer3 result dictionary
        """
        # Log general errors
        if 'PRIMER_ERROR' in result:
            self.logger.warning(f"Primer3 error: {result['PRIMER_ERROR']}")
        
        # Log specific warnings
        if 'PRIMER_WARNING' in result:
            self.logger.warning(f"Primer3 warning: {result['PRIMER_WARNING']}")
        
        # Log left primer issues
        if 'PRIMER_LEFT_EXPLAIN' in result:
            self.logger.debug(f"Left primer issues: {result['PRIMER_LEFT_EXPLAIN']}")
        
        # Log right primer issues
        if 'PRIMER_RIGHT_EXPLAIN' in result:
            self.logger.debug(f"Right primer issues: {result['PRIMER_RIGHT_EXPLAIN']}")
        
        # Log pair issues
        if 'PRIMER_PAIR_EXPLAIN' in result:
            self.logger.debug(f"Primer pair issues: {result['PRIMER_PAIR_EXPLAIN']}")
    
    def design_multiple_primers(self, sequence_info: SequenceInfo, num_primers: int = 3) -> list[Primer]:
        """
        Design multiple primer pairs for a single sequence.
        
        Args:
            sequence_info: SequenceInfo object containing sequence data
            num_primers: Number of primer pairs to return (max 5)
            
        Returns:
            List of Primer objects, sorted by quality (best first)
        """
        self.logger.info(f"Designing {num_primers} primer pairs for: {sequence_info.sequence_id}")
        
        # Validate inputs
        if num_primers < 1 or num_primers > 5:
            raise ValueError("num_primers must be between 1 and 5")
        
        if len(sequence_info.sequence) < 100:
            raise ValueError(f"Sequence too short for primer design: {len(sequence_info.sequence)} bp")
        
        # Prepare primer3 input
        primer3_input = self._prepare_primer3_input(sequence_info)
        primer3_input['PRIMER_NUM_RETURN'] = num_primers
        
        try:
            # Run primer3 design using the new API
            seq_args, global_args = self._split_primer3_args(primer3_input)
            result = primer3.design_primers(seq_args, global_args)
            
            # Check if primers were found
            num_returned = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
            if num_returned == 0:
                self.logger.warning(f"No primers found for sequence: {sequence_info.sequence_id}")
                self._log_primer3_errors(result)
                return []
            
            # Extract all primer pairs
            primers = []
            for i in range(min(num_returned, num_primers)):
                try:
                    primer = self._extract_primer_from_result_index(sequence_info, result, i)
                    primers.append(primer)
                except Exception as e:
                    self.logger.warning(f"Failed to extract primer pair {i}: {e}")
                    continue
            
            self.logger.info(f"Successfully designed {len(primers)} primer pairs for: {sequence_info.sequence_id}")
            return primers
        
        except Exception as e:
            self.logger.error(f"Multiple primer design failed for {sequence_info.sequence_id}: {e}")
            return []
    
    def _extract_primer_from_result_index(self, sequence_info: SequenceInfo, result: Dict[str, Any], index: int) -> Primer:
        """
        Extract primer information from primer3 result at specific index.
        
        Args:
            sequence_info: Original sequence information
            result: primer3 design result dictionary
            index: Index of primer pair to extract
            
        Returns:
            Primer object with extracted information
        """
        # Extract primer pair at specified index
        left_seq = result[f'PRIMER_LEFT_{index}_SEQUENCE']
        right_seq = result[f'PRIMER_RIGHT_{index}_SEQUENCE']
        
        # Extract positions
        left_pos = result[f'PRIMER_LEFT_{index}'][0]
        right_pos = result[f'PRIMER_RIGHT_{index}'][0]
        
        # Extract melting temperatures
        left_tm = result[f'PRIMER_LEFT_{index}_TM']
        right_tm = result[f'PRIMER_RIGHT_{index}_TM']
        
        # Calculate GC content
        left_gc = calculate_gc_content(left_seq)
        right_gc = calculate_gc_content(right_seq)
        
        # Extract penalty scores
        left_penalty = result.get(f'PRIMER_LEFT_{index}_PENALTY', 0.0)
        right_penalty = result.get(f'PRIMER_RIGHT_{index}_PENALTY', 0.0)
        
        # Calculate product size
        product_size = result[f'PRIMER_PAIR_{index}_PRODUCT_SIZE']
        
        return Primer(
            sequence_id=sequence_info.sequence_id,
            forward_primer=left_seq,
            reverse_primer=right_seq,
            tm_forward=round(left_tm, 2),
            tm_reverse=round(right_tm, 2),
            gc_forward=left_gc,
            gc_reverse=right_gc,
            length_forward=len(left_seq),
            length_reverse=len(right_seq),
            product_size=product_size,
            forward_start=left_pos,
            reverse_start=right_pos,
            penalty_forward=left_penalty,
            penalty_reverse=right_penalty
        )
    
    def update_config(self, new_config: PrimerDesignerConfig) -> None:
        """
        Update the primer design configuration.
        
        Args:
            new_config: New PrimerDesignerConfig object
        """
        self.config = new_config
        self._validate_config()
        self.logger.info("Primer design configuration updated")
    
    def get_config_summary(self) -> dict:
        """
        Get a summary of current configuration.
        
        Returns:
            Dictionary with configuration summary
        """
        return {
            'primer_size_range': f"{self.config.primer_min_size}-{self.config.primer_max_size} bp",
            'optimal_size': f"{self.config.primer_opt_size} bp",
            'tm_range': f"{self.config.primer_min_tm}-{self.config.primer_max_tm}°C",
            'optimal_tm': f"{self.config.primer_opt_tm}°C",
            'gc_range': f"{self.config.primer_min_gc}-{self.config.primer_max_gc}%",
            'product_size_range': f"{self.config.product_size_range[0]}-{self.config.product_size_range[1]} bp",
            'salt_concentration': f"{self.config.primer_salt_monovalent} mM",
            'primer_concentration': f"{self.config.primer_dna_conc} nM"
        }
