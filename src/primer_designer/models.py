"""
Data models for primer design results and sequence information.

This module defines Pydantic models for structured data representation
of primer design results and input sequences.
"""

from typing import Optional
from pydantic import BaseModel, Field, field_validator


class SequenceInfo(BaseModel):
    """
    Model representing input DNA sequence information.
    
    Attributes:
        sequence_id: Unique identifier for the sequence
        sequence: DNA sequence string (A, T, G, C)
        description: Optional description from FASTA header
        length: Length of the sequence in base pairs
    """
    
    sequence_id: str = Field(..., min_length=1, description="Unique sequence identifier")
    sequence: str = Field(..., min_length=50, description="DNA sequence")
    description: Optional[str] = Field(default=None, description="Sequence description")
    length: int = Field(..., gt=0, description="Sequence length in bp")
    
    @field_validator('sequence')
    @classmethod
    def validate_dna_sequence(cls, v):
        """Validate that sequence contains only valid DNA bases."""
        v = v.upper().strip()
        valid_bases = set('ATGCRYSWKMBDHVN')  # Include ambiguous bases
        invalid_bases = set(v) - valid_bases
        if invalid_bases:
            raise ValueError(f'Invalid DNA bases found: {invalid_bases}')
        return v
    
    @field_validator('length')
    @classmethod
    def validate_length_matches_sequence(cls, v, info):
        """Validate that length matches actual sequence length."""
        if info.data.get('sequence') and len(info.data['sequence']) != v:
            raise ValueError('Length does not match sequence length')
        return v
    
    @classmethod
    def from_fasta_record(cls, record) -> 'SequenceInfo':
        """
        Create SequenceInfo from BioPython SeqRecord.
        
        Args:
            record: BioPython SeqRecord object
            
        Returns:
            SequenceInfo instance
        """
        return cls(
            sequence_id=record.id,
            sequence=str(record.seq),
            description=record.description,
            length=len(record.seq)
        )


class Primer(BaseModel):
    """
    Model representing primer design results for a single sequence.
    
    Attributes:
        sequence_id: ID of the template sequence
        forward_primer: Forward primer sequence
        reverse_primer: Reverse primer sequence
        tm_forward: Melting temperature of forward primer (°C)
        tm_reverse: Melting temperature of reverse primer (°C)
        gc_forward: GC content of forward primer (%)
        gc_reverse: GC content of reverse primer (%)
        length_forward: Length of forward primer (bp)
        length_reverse: Length of reverse primer (bp)
        product_size: Expected PCR product size (bp)
        forward_start: Start position of forward primer (0-based)
        reverse_start: Start position of reverse primer (0-based)
        penalty_forward: Primer3 penalty score for forward primer
        penalty_reverse: Primer3 penalty score for reverse primer
    """
    
    sequence_id: str = Field(..., description="Template sequence identifier")
    
    # Primer sequences
    forward_primer: str = Field(..., min_length=15, description="Forward primer sequence")
    reverse_primer: str = Field(..., min_length=15, description="Reverse primer sequence")
    
    # Melting temperatures
    tm_forward: float = Field(..., ge=40.0, le=80.0, description="Forward primer Tm (°C)")
    tm_reverse: float = Field(..., ge=40.0, le=80.0, description="Reverse primer Tm (°C)")
    
    # GC content
    gc_forward: float = Field(..., ge=0.0, le=100.0, description="Forward primer GC content (%)")
    gc_reverse: float = Field(..., ge=0.0, le=100.0, description="Reverse primer GC content (%)")
    
    # Lengths
    length_forward: int = Field(..., ge=15, le=50, description="Forward primer length (bp)")
    length_reverse: int = Field(..., ge=15, le=50, description="Reverse primer length (bp)")
    
    # Product information
    product_size: int = Field(..., ge=50, description="PCR product size (bp)")
    
    # Position information
    forward_start: int = Field(..., ge=0, description="Forward primer start position (0-based)")
    reverse_start: int = Field(..., ge=0, description="Reverse primer start position (0-based)")
    
    # Quality scores
    penalty_forward: float = Field(..., ge=0.0, description="Forward primer penalty score")
    penalty_reverse: float = Field(..., ge=0.0, description="Reverse primer penalty score")
    
    @field_validator('forward_primer', 'reverse_primer')
    @classmethod
    def validate_primer_sequence(cls, v):
        """Validate primer sequences contain only valid DNA bases."""
        v = v.upper().strip()
        valid_bases = set('ATGC')
        invalid_bases = set(v) - valid_bases
        if invalid_bases:
            raise ValueError(f'Invalid DNA bases in primer: {invalid_bases}')
        return v
    
    @field_validator('length_forward')
    @classmethod
    def validate_forward_length(cls, v, info):
        """Validate forward primer length matches sequence length."""
        if info.data.get('forward_primer') and len(info.data['forward_primer']) != v:
            raise ValueError('Forward primer length does not match sequence length')
        return v
    
    @field_validator('length_reverse')
    @classmethod
    def validate_reverse_length(cls, v, info):
        """Validate reverse primer length matches sequence length."""
        if info.data.get('reverse_primer') and len(info.data['reverse_primer']) != v:
            raise ValueError('Reverse primer length does not match sequence length')
        return v
    
    def to_dict(self) -> dict:
        """
        Convert primer information to dictionary format.
        
        Returns:
            Dictionary representation of primer data
        """
        return {
            'sequence_id': self.sequence_id,
            'forward_primer': self.forward_primer,
            'reverse_primer': self.reverse_primer,
            'tm_forward': round(self.tm_forward, 2),
            'tm_reverse': round(self.tm_reverse, 2),
            'gc_forward': round(self.gc_forward, 2),
            'gc_reverse': round(self.gc_reverse, 2),
            'length_forward': self.length_forward,
            'length_reverse': self.length_reverse,
            'product_size': self.product_size,
            'forward_start': self.forward_start,
            'reverse_start': self.reverse_start,
            'penalty_forward': round(self.penalty_forward, 3),
            'penalty_reverse': round(self.penalty_reverse, 3)
        }
    
    @property
    def tm_difference(self) -> float:
        """Calculate the difference between forward and reverse primer Tm values."""
        return abs(self.tm_forward - self.tm_reverse)
    
    @property
    def average_tm(self) -> float:
        """Calculate the average Tm of forward and reverse primers."""
        return (self.tm_forward + self.tm_reverse) / 2.0
    
    @property
    def gc_difference(self) -> float:
        """Calculate the difference between forward and reverse primer GC content."""
        return abs(self.gc_forward - self.gc_reverse)
    
    @property
    def total_penalty(self) -> float:
        """Calculate the total penalty score for the primer pair."""
        return self.penalty_forward + self.penalty_reverse
