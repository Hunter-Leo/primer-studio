# Primer Studio

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![GitHub issues](https://img.shields.io/github/issues/Hunter-Leo/primer-studio)](https://github.com/Hunter-Leo/primer-studio/issues)
[![GitHub stars](https://img.shields.io/github/stars/Hunter-Leo/primer-studio)](https://github.com/Hunter-Leo/primer-studio/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/Hunter-Leo/primer-studio)](https://github.com/Hunter-Leo/primer-studio/network)

A batch primer design tool for PCR applications that supports automatic primer design for multiple DNA sequences with configurable parameters and high-performance processing.

## Demo

![Primer Studio Demo](primer-demo.gif)

## Features

- **Batch Processing**: Design primers for multiple sequences from FASTA files
- **Configurable Parameters**: Pre-defined configurations for different template types (standard, GC-rich, AT-rich)
- **🆕 Custom Configuration**: Create and use custom JSON configuration files with full parameter control
- **🆕 Parameter Documentation**: Built-in parameter information with types, ranges, and constraints
- **Multiple Output Formats**: Export results in CSV or JSON format
- **Memory Efficient**: Support for memory-efficient processing of large datasets
- **Parallel Processing**: Multi-threaded processing for improved performance
- **Command Line Interface**: Easy-to-use CLI with rich output formatting and emoji support
- **Python API**: Can be used as a Python package for integration into other tools
- **Comprehensive Validation**: Input validation and detailed error reporting

## Installation

### Prerequisites

Make sure you have [uv](https://docs.astral.sh/uv/getting-started/installation/) installed on your system.

### Install from Source (Development)

```bash
git clone https://github.com/Hunter-Leo/primer-studio.git
cd primer-studio
uv sync
uv build
uv pip install -e .
```

### Install as Tool (Recommended for Users)

```bash
# Install directly from repository
uv tool install git+https://github.com/Hunter-Leo/primer-studio.git

# Upgrade to latest version
uv tool upgrade primer-studio
```

### Alternative Installation Methods

#### Using pip (Traditional)

```bash
git clone https://github.com/Hunter-Leo/primer-studio.git
cd primer-studio
pip install -e .
```

#### Using uvx (Run without Installation)

```bash
# Run directly without installing
uvx --from primer-studio primer-designer --help
```

#### Direct Usage
You can directly use the `primer-designer` command after installation with any of the above methods.

```bash
primer-designer --help
```

## Quick Start

### Command Line Usage

1. **Design primers for a FASTA file:**
```bash
uv run python -m primer_designer.cli design input.fasta --output results.csv
```

2. **Use different configuration presets:**
```bash
# For GC-rich templates
uv run python -m primer_designer.cli design input.fasta --config gc-rich --output results.csv

# For AT-rich templates  
uv run python -m primer_designer.cli design input.fasta --config gc-poor --output results.csv
```

3. **🆕 Use custom configuration:**
```bash
# Export a configuration template
uv run python -m primer_designer.cli export-template -o my_config.json

# Edit my_config.json with your custom parameters
# Then use it for primer design
uv run python -m primer_designer.cli design input.fasta --custom-config my_config.json
```

4. **Memory-efficient processing for large files:**
```bash
uv run python -m primer_designer.cli design large_file.fasta --memory-efficient --output results.csv
```

5. **Get batch statistics:**
```bash
uv run python -m primer_designer.cli design input.fasta --output results.csv --stats
```

6. **Validate FASTA file:**
```bash
uv run python -m primer_designer.cli validate input.fasta
```

7. **List available configurations:**
```bash
uv run python -m primer_designer.cli list-configs
```

8. **🆕 View parameter information:**
```bash
uv run python -m primer_designer.cli export-template --show-params
```

### 🆕 Custom Configuration Workflow

1. **Export configuration template:**
```bash
# Export default template
uv run python -m primer_designer.cli export-template -o my_config.json

# Export template based on existing preset
uv run python -m primer_designer.cli export-template -o custom.json --preset gc-rich

# Export with parameter documentation
uv run python -m primer_designer.cli export-template -o documented.json --documented
```

2. **View parameter specifications:**
```bash
# Show all parameters with types, ranges, and descriptions
uv run python -m primer_designer.cli export-template --show-params
```

3. **Edit configuration file:**
Edit the JSON file with your desired parameter values. All parameters include validation constraints.

4. **Use custom configuration:**
```bash
uv run python -m primer_designer.cli design input.fasta --custom-config my_config.json
```

### Python API Usage

```python
from primer_designer import BatchDesigner, FastaReader, STANDARD_CONFIG, GC_RICH_CONFIG

# Initialize batch designer
batch_designer = BatchDesigner(config=STANDARD_CONFIG)

# Design primers from FASTA file
primers = batch_designer.design_primers_from_fasta("input.fasta", "output.csv")

# Use different configuration
batch_designer_gc = BatchDesigner(config=GC_RICH_CONFIG)
primers = batch_designer_gc.design_primers_from_fasta("gc_rich_sequences.fasta")

# 🆕 Use custom configuration from JSON file
from primer_designer import PrimerDesignerConfig

custom_config = PrimerDesignerConfig.from_json_file("my_config.json")
batch_designer_custom = BatchDesigner(config=custom_config)
primers = batch_designer_custom.design_primers_from_fasta("input.fasta")

# Process individual sequences
from primer_designer import PrimerDesigner, SequenceInfo

designer = PrimerDesigner(STANDARD_CONFIG)
sequence = SequenceInfo(
    sequence_id="test_seq",
    sequence="ATGCGATCGATCG...",  # Your DNA sequence
    length=len("ATGCGATCGATCG...")
)
primer = designer.design_primer(sequence)
```

## Configuration Presets

### Standard Configuration
- **Use case**: General purpose PCR primers
- **Primer size**: 18-25 bp (optimal: 20 bp)
- **Tm range**: 57-63°C (optimal: 60°C)
- **GC content**: 40-60%
- **Product size**: 100-1000 bp

### GC-Rich Configuration
- **Use case**: GC-rich templates (>60% GC)
- **Primer size**: 20-27 bp (optimal: 22 bp)
- **Tm range**: 59-65°C (optimal: 62°C)
- **GC content**: 50-70%
- **Product size**: 150-800 bp

### GC-Poor Configuration
- **Use case**: AT-rich templates (<40% GC)
- **Primer size**: 16-22 bp (optimal: 18 bp)
- **Tm range**: 55-61°C (optimal: 58°C)
- **GC content**: 30-50%
- **Product size**: 100-1200 bp

## 🆕 Custom Configuration

### Overview
Create fully customized primer design parameters using JSON configuration files. This feature allows precise control over all primer design parameters beyond the predefined presets.

### Parameter Categories

#### 🧬 Primer Size Parameters
- `primer_opt_size` (15-35 bp): Optimal primer length
- `primer_min_size` (15-35 bp): Minimum primer length  
- `primer_max_size` (15-35 bp): Maximum primer length

#### 🌡️ Melting Temperature Parameters
- `primer_opt_tm` (50.0-70.0°C): Optimal primer Tm
- `primer_min_tm` (50.0-70.0°C): Minimum primer Tm
- `primer_max_tm` (50.0-70.0°C): Maximum primer Tm

#### 🔬 GC Content Parameters
- `primer_opt_gc_percent` (20.0-80.0%): Optimal GC content
- `primer_min_gc` (20.0-80.0%): Minimum GC content
- `primer_max_gc` (20.0-80.0%): Maximum GC content

#### 🔗 Secondary Structure Parameters
- `primer_max_poly_x` (3-6 bases): Maximum mononucleotide repeat
- `primer_max_self_any` (0.0-12.0): Maximum self-complementarity score
- `primer_max_self_end` (0.0-8.0): Maximum 3' self-complementarity score

#### 🧪 Salt and Concentration Parameters
- `primer_salt_monovalent` (0.0-1000.0 mM): Monovalent salt concentration
- `primer_salt_divalent` (0.0-10.0 mM): Divalent salt concentration
- `primer_dntp_conc` (0.0-10.0 mM): dNTP concentration
- `primer_dna_conc` (0.0-1000.0 nM): Primer concentration

#### 📏 Product Size Parameters
- `product_size_range` ([min, max] bp): PCR product size range

### Configuration File Format
```json
{
  "_metadata": {
    "description": "Custom Primer Designer Configuration",
    "version": "1.0"
  },
  "config": {
    "primer_opt_size": 22,
    "primer_min_size": 20,
    "primer_max_size": 28,
    "primer_opt_tm": 65.0,
    "primer_min_tm": 62.0,
    "primer_max_tm": 68.0,
    "primer_opt_gc_percent": 55.0,
    "primer_min_gc": 45.0,
    "primer_max_gc": 65.0,
    "primer_max_poly_x": 3,
    "primer_max_self_any": 6.0,
    "primer_max_self_end": 2.0,
    "primer_salt_monovalent": 75.0,
    "primer_salt_divalent": 2.0,
    "primer_dntp_conc": 0.8,
    "primer_dna_conc": 100.0,
    "product_size_range": [200, 800]
  }
}
```

### Usage Examples

#### View All Parameters
```bash
# Display detailed parameter information
uv run python -m primer_designer.cli export-template --show-params
```

#### Create Custom Configuration
```bash
# Export template with documentation
uv run python -m primer_designer.cli export-template -o my_config.json --documented

# Edit my_config.json with your custom values
# Use custom configuration
uv run python -m primer_designer.cli design input.fasta --custom-config my_config.json
```

#### Base on Existing Preset
```bash
# Start with GC-rich preset and customize
uv run python -m primer_designer.cli export-template -o custom.json --preset gc-rich
# Edit custom.json as needed
uv run python -m primer_designer.cli design input.fasta --custom-config custom.json
```

## Output Format

### CSV Output
```csv
sequence_id,forward_primer,reverse_primer,tm_forward,tm_reverse,gc_forward,gc_reverse,length_forward,length_reverse,product_size,forward_start,reverse_start,penalty_forward,penalty_reverse
test_seq_001,TCTTGCTGGCCGTCAGAATT,TTCGCCCCTACTAGACGTGA,59.96,60.04,50.0,55.0,20,20,640,20,659,0.037,0.035
```

### JSON Output
```json
{
  "metadata": {
    "total_primers": 1,
    "config_summary": {
      "primer_size_range": "18-25 bp",
      "tm_range": "57.0-63.0°C",
      "gc_range": "40.0-60.0%"
    }
  },
  "primers": [
    {
      "sequence_id": "test_seq_001",
      "forward_primer": "TCTTGCTGGCCGTCAGAATT",
      "reverse_primer": "TTCGCCCCTACTAGACGTGA",
      "tm_forward": 59.96,
      "tm_reverse": 60.04,
      "gc_forward": 50.0,
      "gc_reverse": 55.0,
      "length_forward": 20,
      "length_reverse": 20,
      "product_size": 640,
      "forward_start": 20,
      "reverse_start": 659,
      "penalty_forward": 0.037,
      "penalty_reverse": 0.035
    }
  ]
}
```

## Command Line Options

### Design Command
```bash
uv run python -m primer_designer.cli design [OPTIONS] INPUT_FILE

Options:
  --output, -o PATH           Output file path
  --format, -f TEXT          Output format (csv/json) [default: csv]
  --config, -c TEXT          Configuration preset [default: standard]
  --custom-config PATH       🆕 Path to custom JSON configuration file (overrides --config)
  --parallel/--no-parallel   Use parallel processing [default: parallel]
  --workers, -w INTEGER      Maximum worker threads
  --memory-efficient         Use memory-efficient processing
  --log-level TEXT           Logging level [default: INFO]
  --log-file PATH            Log file path
  --stats                    Show batch statistics
  --validate-only            Only validate input file
```

### 🆕 Export Template Command
```bash
uv run python -m primer_designer.cli export-template [OPTIONS]

Options:
  --output, -o PATH          Output JSON file path for configuration template
  --preset, -p TEXT          Base template on existing preset (standard/gc-rich/gc-poor)
  --show                     Display template in console instead of saving to file
  --show-params              Display detailed parameter information with types and ranges
  --documented               Include parameter documentation in exported template
```

### Other Commands
```bash
# Validate FASTA file
uv run python -m primer_designer.cli validate INPUT_FILE

# Show configuration details
uv run python -m primer_designer.cli config-info [PRESET_NAME]

# List all configurations
uv run python -m primer_designer.cli list-configs

# 🆕 Export configuration template
uv run python -m primer_designer.cli export-template [OPTIONS]
```

## Performance

- **Standard mode**: Processes sequences in parallel using multiple threads
- **Memory-efficient mode**: Processes sequences one by one, suitable for very large files
- **Typical performance**: 
  - Small files (1-100 sequences): < 1 minute
  - Medium files (100-1000 sequences): 1-10 minutes
  - Large files (1000+ sequences): 10+ minutes (use memory-efficient mode)

## Dependencies

- **Python**: >=3.12
- **BioPython**: >=1.85 (for sequence handling)
- **primer3-py**: >=2.2.0 (for primer design)
- **Pydantic**: >=2.0.0 (for data validation)
- **Typer**: >=0.9.0 (for CLI)
- **Rich**: >=13.0.0 (for rich terminal output)

## Testing

Generate test data:
```bash
python generate_test_data.py --all
```

Run tests:
```bash
uv run python -m pytest tests/ -v
```

Test with sample data:
```bash
# Test single sequence
uv run python -m primer_designer.cli design test_data/single_sequence_test.fasta --stats

# Test batch processing
uv run python -m primer_designer.cli design test_data/standard_test.fasta --stats

# Test different configurations
uv run python -m primer_designer.cli design test_data/gc_rich_test.fasta --config gc-rich

# 🆕 Test custom configuration
uv run python -m primer_designer.cli export-template -o test_config.json --preset gc-rich
uv run python -m primer_designer.cli design test_data/standard_test.fasta --custom-config test_config.json

# 🆕 View parameter information
uv run python -m primer_designer.cli export-template --show-params
```

## Architecture

The project follows a modular, object-oriented design:

```
src/primer_designer/
├── __init__.py          # Package initialization
├── config.py            # Configuration management
├── models.py            # Pydantic data models
├── primer.py            # Core primer design logic
├── fasta_reader.py      # FASTA file handling
├── batch_designer.py    # Batch processing
├── utils.py             # Utility functions
└── cli.py               # Command line interface
```

### Key Classes

- **PrimerDesignerConfig**: Configuration management with validation and JSON import/export
- **PrimerDesigner**: Core primer design using primer3-py
- **BatchDesigner**: Batch processing with parallel support
- **FastaReader**: FASTA file reading and validation
- **Primer/SequenceInfo**: Data models with Pydantic validation

### 🆕 New Configuration Features

- **JSON Import/Export**: Load and save custom configurations
- **Parameter Documentation**: Auto-generated parameter info from Pydantic models
- **Template Generation**: Create configuration templates based on presets
- **Validation**: Comprehensive parameter validation with helpful error messages

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Run the test suite
6. Submit a pull request

## License

This project is licensed under the GNU General Public License v2.0 (GPL-2.0) - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Primer3**: Core primer design algorithm
- **BioPython**: Sequence analysis tools
- **Pydantic**: Data validation framework
- **Typer**: Modern CLI framework

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=Hunter-Leo/primer-studio&type=Date)](https://www.star-history.com/#Hunter-Leo/primer-studio&Date)
