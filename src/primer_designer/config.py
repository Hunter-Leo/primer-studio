"""
Configuration module for primer design parameters.

This module defines the PrimerDesignerConfig class using Pydantic v2 for
parameter validation and provides several predefined configurations for
common use cases.
"""

import json
from pathlib import Path
from typing import Optional, Union
from pydantic import BaseModel, Field, field_validator


class PrimerDesignerConfig(BaseModel):
    """
    Configuration class for primer design parameters.
    
    This class encapsulates all the parameters needed for primer3-py
    primer design, with validation and default values suitable for
    most PCR applications.
    
    Attributes:
        primer_opt_size: Optimal primer length
        primer_min_size: Minimum primer length  
        primer_max_size: Maximum primer length
        primer_opt_tm: Optimal primer melting temperature (°C)
        primer_min_tm: Minimum primer melting temperature (°C)
        primer_max_tm: Maximum primer melting temperature (°C)
        primer_opt_gc_percent: Optimal GC content percentage
        primer_min_gc: Minimum GC content percentage
        primer_max_gc: Maximum GC content percentage
        primer_max_poly_x: Maximum length of mononucleotide repeat
        primer_max_self_any: Maximum self-complementarity score
        primer_max_self_end: Maximum 3' self-complementarity score
        primer_salt_monovalent: Monovalent salt concentration (mM)
        primer_salt_divalent: Divalent salt concentration (mM)
        primer_dntp_conc: dNTP concentration (mM)
        primer_dna_conc: Primer concentration (nM)
        product_size_range: PCR product size range as tuple (min, max)
    """
    
    # Primer size parameters
    primer_opt_size: int = Field(default=20, ge=15, le=35, description="Optimal primer length")
    primer_min_size: int = Field(default=18, ge=15, le=35, description="Minimum primer length")
    primer_max_size: int = Field(default=25, ge=15, le=35, description="Maximum primer length")
    
    # Melting temperature parameters (°C)
    primer_opt_tm: float = Field(default=60.0, ge=50.0, le=70.0, description="Optimal primer Tm")
    primer_min_tm: float = Field(default=57.0, ge=50.0, le=70.0, description="Minimum primer Tm")
    primer_max_tm: float = Field(default=63.0, ge=50.0, le=70.0, description="Maximum primer Tm")
    
    # GC content parameters (%)
    primer_opt_gc_percent: Optional[float] = Field(default=50.0, ge=20.0, le=80.0, description="Optimal GC content")
    primer_min_gc: float = Field(default=40.0, ge=20.0, le=80.0, description="Minimum GC content")
    primer_max_gc: float = Field(default=60.0, ge=20.0, le=80.0, description="Maximum GC content")
    
    # Secondary structure parameters
    primer_max_poly_x: int = Field(default=4, ge=3, le=6, description="Maximum mononucleotide repeat length")
    primer_max_self_any: float = Field(default=8.0, ge=0.0, le=12.0, description="Maximum self-complementarity")
    primer_max_self_end: float = Field(default=3.0, ge=0.0, le=8.0, description="Maximum 3' self-complementarity")
    
    # Salt and concentration parameters
    primer_salt_monovalent: float = Field(default=50.0, ge=0.0, le=1000.0, description="Monovalent salt concentration (mM)")
    primer_salt_divalent: float = Field(default=1.5, ge=0.0, le=10.0, description="Divalent salt concentration (mM)")
    primer_dntp_conc: float = Field(default=0.6, ge=0.0, le=10.0, description="dNTP concentration (mM)")
    primer_dna_conc: float = Field(default=50.0, ge=0.0, le=1000.0, description="Primer concentration (nM)")
    
    # Product size range
    product_size_range: tuple[int, int] = Field(default=(100, 1000), description="PCR product size range (min, max)")
    
    @field_validator('primer_max_size')
    @classmethod
    def validate_size_range(cls, v, info):
        """Validate that size parameters are in correct order."""
        if info.data.get('primer_min_size') and v < info.data['primer_min_size']:
            raise ValueError('primer_max_size must be >= primer_min_size')
        return v
    
    @field_validator('primer_max_tm')
    @classmethod
    def validate_tm_range(cls, v, info):
        """Validate that Tm parameters are in correct order."""
        if info.data.get('primer_min_tm') and v < info.data['primer_min_tm']:
            raise ValueError('primer_max_tm must be >= primer_min_tm')
        return v
    
    @field_validator('primer_max_gc')
    @classmethod
    def validate_gc_range(cls, v, info):
        """Validate that GC parameters are in correct order."""
        if info.data.get('primer_min_gc') and v < info.data['primer_min_gc']:
            raise ValueError('primer_max_gc must be >= primer_min_gc')
        return v
    
    @field_validator('product_size_range')
    @classmethod
    def validate_product_size_range(cls, v):
        """Validate that product size range is valid."""
        if len(v) != 2:
            raise ValueError('product_size_range must be a tuple of (min, max)')
        if v[0] >= v[1]:
            raise ValueError('product_size_range max must be > min')
        if v[0] < 50:
            raise ValueError('minimum product size should be >= 50 bp')
        return v
    
    def to_primer3_dict(self) -> dict:
        """
        Convert configuration to primer3-py compatible dictionary.
        
        Returns:
            Dictionary with primer3-py parameter names and values
        """
        return {
            'PRIMER_OPT_SIZE': self.primer_opt_size,
            'PRIMER_MIN_SIZE': self.primer_min_size,
            'PRIMER_MAX_SIZE': self.primer_max_size,
            'PRIMER_OPT_TM': self.primer_opt_tm,
            'PRIMER_MIN_TM': self.primer_min_tm,
            'PRIMER_MAX_TM': self.primer_max_tm,
            'PRIMER_OPT_GC_PERCENT': self.primer_opt_gc_percent,
            'PRIMER_MIN_GC': self.primer_min_gc,
            'PRIMER_MAX_GC': self.primer_max_gc,
            'PRIMER_MAX_POLY_X': self.primer_max_poly_x,
            'PRIMER_MAX_SELF_ANY': self.primer_max_self_any,
            'PRIMER_MAX_SELF_END': self.primer_max_self_end,
            'PRIMER_SALT_MONOVALENT': self.primer_salt_monovalent,
            'PRIMER_SALT_DIVALENT': self.primer_salt_divalent,
            'PRIMER_DNTP_CONC': self.primer_dntp_conc,
            'PRIMER_DNA_CONC': self.primer_dna_conc,
            'PRIMER_PRODUCT_SIZE_RANGE': [list(self.product_size_range)],
        }
    
    def to_json_file(self, file_path: Union[str, Path]) -> None:
        """
        Export configuration to JSON file.
        
        Args:
            file_path: Path where to save the JSON configuration file
            
        Raises:
            IOError: If file cannot be written
        """
        file_path = Path(file_path)
        
        # Create directory if it doesn't exist
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to dictionary and add metadata
        config_dict = {
            "_metadata": {
                "description": "Primer Designer Configuration File",
                "version": "1.0",
                "created_by": "primer-designer CLI"
            },
            "config": self.model_dump()
        }
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(config_dict, f, indent=2, ensure_ascii=False)
        except IOError as e:
            raise IOError(f"Failed to write configuration to {file_path}: {e}")
    
    @classmethod
    def from_json_file(cls, file_path: Union[str, Path]) -> 'PrimerDesignerConfig':
        """
        Load configuration from JSON file.
        
        Args:
            file_path: Path to the JSON configuration file
            
        Returns:
            PrimerDesignerConfig instance loaded from file
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If JSON is invalid or configuration is malformed
            IOError: If file cannot be read
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {file_path}")
        
        if not file_path.is_file():
            raise ValueError(f"Path is not a file: {file_path}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in configuration file {file_path}: {e}")
        except IOError as e:
            raise IOError(f"Failed to read configuration file {file_path}: {e}")
        
        # Handle both new format (with metadata) and simple format
        if isinstance(data, dict) and "config" in data:
            config_data = data["config"]
        else:
            config_data = data
        
        try:
            return cls(**config_data)
        except Exception as e:
            raise ValueError(f"Invalid configuration data in {file_path}: {e}")
    
    def to_json_string(self, indent: int = 2) -> str:
        """
        Convert configuration to JSON string.
        
        Args:
            indent: Number of spaces for indentation
            
        Returns:
            JSON string representation of the configuration
        """
        config_dict = {
            "_metadata": {
                "description": "Primer Designer Configuration",
                "version": "1.0"
            },
            "config": self.model_dump()
        }
        return json.dumps(config_dict, indent=indent, ensure_ascii=False)
    
    @classmethod
    def create_template(cls) -> 'PrimerDesignerConfig':
        """
        Create a template configuration with default values and helpful comments.
        
        Returns:
            PrimerDesignerConfig instance with default values
        """
        return cls()  # Uses all default values from Field definitions
    
    @classmethod
    def get_field_info(cls) -> dict:
        """
        Get detailed information about all configuration fields.
        
        Returns:
            Dictionary containing field information including types, constraints, and descriptions
        """
        # Manual mapping of field constraints (extracted from Field definitions above)
        field_constraints = {
            "primer_opt_size": {"min": 15, "max": 35, "unit": "bp"},
            "primer_min_size": {"min": 15, "max": 35, "unit": "bp"},
            "primer_max_size": {"min": 15, "max": 35, "unit": "bp"},
            "primer_opt_tm": {"min": 50.0, "max": 70.0, "unit": "°C"},
            "primer_min_tm": {"min": 50.0, "max": 70.0, "unit": "°C"},
            "primer_max_tm": {"min": 50.0, "max": 70.0, "unit": "°C"},
            "primer_opt_gc_percent": {"min": 20.0, "max": 80.0, "unit": "%"},
            "primer_min_gc": {"min": 20.0, "max": 80.0, "unit": "%"},
            "primer_max_gc": {"min": 20.0, "max": 80.0, "unit": "%"},
            "primer_max_poly_x": {"min": 3, "max": 6, "unit": "bases"},
            "primer_max_self_any": {"min": 0.0, "max": 12.0, "unit": "score"},
            "primer_max_self_end": {"min": 0.0, "max": 8.0, "unit": "score"},
            "primer_salt_monovalent": {"min": 0.0, "max": 1000.0, "unit": "mM"},
            "primer_salt_divalent": {"min": 0.0, "max": 10.0, "unit": "mM"},
            "primer_dntp_conc": {"min": 0.0, "max": 10.0, "unit": "mM"},
            "primer_dna_conc": {"min": 0.0, "max": 1000.0, "unit": "nM"},
            "product_size_range": {"format": "[min, max] tuple", "unit": "bp", "min_value": 50}
        }
        
        field_info = {}
        
        for field_name, field_info_obj in cls.model_fields.items():
            # Get basic info
            info = {
                "description": field_info_obj.description or "No description available",
                "type": str(field_info_obj.annotation).replace("typing.", "").replace("Union[", "").replace(", NoneType]", " | None"),
                "default": field_info_obj.default,
                "constraints": {}
            }
            
            # Add constraint information
            if field_name in field_constraints:
                constraints = field_constraints[field_name]
                if "min" in constraints and "max" in constraints:
                    info["constraints"]["range"] = f"{constraints['min']}-{constraints['max']}"
                elif "format" in constraints:
                    info["constraints"]["format"] = constraints["format"]
                
                if "unit" in constraints:
                    info["constraints"]["unit"] = constraints["unit"]
                
                if "min_value" in constraints:
                    info["constraints"]["min_value"] = constraints["min_value"]
            
            # Special handling for specific field types
            if field_name == "product_size_range":
                info["type"] = "tuple[int, int]"
            elif "Optional" in str(field_info_obj.annotation):
                info["type"] = info["type"].replace("Optional[", "").replace("]", " | None")
            
            field_info[field_name] = info
        
        return field_info
    
    @classmethod
    def create_documented_template(cls) -> dict:
        """
        Create a template with inline documentation for all parameters.
        
        Returns:
            Dictionary with configuration and parameter documentation
        """
        field_info = cls.get_field_info()
        template_config = cls().model_dump()
        
        documented_config = {
            "_metadata": {
                "description": "Primer Designer Configuration Template with Parameter Documentation",
                "version": "1.0",
                "created_by": "primer-designer CLI",
                "parameter_documentation": "See _parameter_info section for detailed parameter descriptions"
            },
            "_parameter_info": {},
            "config": template_config
        }
        
        # Add parameter documentation
        for field_name, info in field_info.items():
            doc = {
                "description": info["description"],
                "type": info["type"],
                "default_value": info["default"],
                "current_value": template_config[field_name]
            }
            
            if info["constraints"]:
                doc["constraints"] = info["constraints"]
            
            documented_config["_parameter_info"][field_name] = doc
        
        return documented_config


# Predefined configurations for common use cases

STANDARD_CONFIG = PrimerDesignerConfig(
    primer_opt_size=20,
    primer_min_size=18,
    primer_max_size=25,
    primer_opt_tm=60.0,
    primer_min_tm=57.0,
    primer_max_tm=63.0,
    primer_min_gc=40.0,
    primer_max_gc=60.0,
    product_size_range=(100, 1000)
)
"""Standard configuration suitable for most PCR applications."""

GC_RICH_CONFIG = PrimerDesignerConfig(
    primer_opt_size=22,
    primer_min_size=20,
    primer_max_size=27,
    primer_opt_tm=62.0,
    primer_min_tm=59.0,
    primer_max_tm=65.0,
    primer_min_gc=50.0,
    primer_max_gc=70.0,
    primer_max_poly_x=3,
    product_size_range=(150, 800)
)
"""Configuration optimized for GC-rich templates."""

GC_POOR_CONFIG = PrimerDesignerConfig(
    primer_opt_size=18,
    primer_min_size=16,
    primer_max_size=22,
    primer_opt_tm=58.0,
    primer_min_tm=55.0,
    primer_max_tm=61.0,
    primer_min_gc=30.0,
    primer_max_gc=50.0,
    primer_max_poly_x=5,
    product_size_range=(100, 1200)
)
"""Configuration optimized for AT-rich templates."""
