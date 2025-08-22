"""
Command-line interface for the primer designer tool.

This module provides a Typer-based CLI for batch primer design
with configurable parameters and output options.
"""

import json
import sys
from pathlib import Path
from typing import Optional, Annotated

import typer
import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich import print as rprint

from .config import PrimerDesignerConfig, STANDARD_CONFIG, GC_RICH_CONFIG, GC_POOR_CONFIG
from .batch_designer import BatchDesigner
from .fasta_reader import FastaReader
from .utils import setup_logging


# Initialize Typer app
app = typer.Typer(
    name="primer-designer",
    help="""🧬 High-performance batch primer design tool for PCR applications.

Primer Studio provides automated primer design for multiple DNA sequences with
configurable parameters, parallel processing, and comprehensive validation.

✨ FEATURES:

\b
• 📁 Batch processing of FASTA files with multiple sequences
• ⚙️ Three optimized configuration presets (standard, GC-rich, GC-poor)  
• 🚀 Parallel processing for improved performance
• 💾 Memory-efficient mode for large datasets
• ✅ Comprehensive input validation and error reporting
• 📊 Multiple output formats (CSV, JSON) with detailed statistics
• 🔬 Built on primer3 algorithm for reliable primer design

🚀 QUICK START:

\b
  primer-designer design sequences.fasta                    # Basic usage
  primer-designer design input.fasta --config gc-rich      # Use GC-rich preset
  primer-designer validate input.fasta                     # Validate FASTA file
  primer-designer list-configs                             # Show all presets

💡 For detailed help on any command, use: primer-designer COMMAND --help""",
    add_completion=False
)

# Initialize Rich console
console = Console()

# Predefined configurations
CONFIGS = {
    "standard": STANDARD_CONFIG,
    "gc-rich": GC_RICH_CONFIG,
    "gc-poor": GC_POOR_CONFIG
}


@app.command()
def design(
    input_file: Annotated[Path, typer.Argument(
        help="📁 Input FASTA file containing DNA sequences (minimum 50 bp each) for primer design"
    )],
    output_file: Annotated[Optional[Path], typer.Option(
        "--output", "-o", 
        help="💾 Output file path. If not specified, defaults to '<input_file>.primers.<format>'"
    )] = None,
    output_format: Annotated[str, typer.Option(
        "--format", "-f", 
        help="📊 Output format for primer results",
        click_type=click.Choice(["csv", "json"], case_sensitive=False)
    )] = "csv",
    config_preset: Annotated[str, typer.Option(
        "--config", "-c", 
        help="⚙️ Primer design configuration preset (ignored if --custom-config is used)",
        click_type=click.Choice(["standard", "gc-rich", "gc-poor"], case_sensitive=False)
    )] = "standard",
    custom_config: Annotated[Optional[Path], typer.Option(
        "--custom-config", 
        help="📄 Path to custom JSON configuration file (overrides --config preset)"
    )] = None,
    parallel: Annotated[bool, typer.Option(
        "--parallel/--no-parallel", 
        help="🚀 Enable parallel processing using multiple CPU cores for faster execution"
    )] = True,
    max_workers: Annotated[Optional[int], typer.Option(
        "--workers", "-w", 
        help="👥 Maximum number of worker threads (1-32). Auto-detects CPU cores if not specified",
        min=1, max=32
    )] = None,
    memory_efficient: Annotated[bool, typer.Option(
        "--memory-efficient", 
        help="💾 Process sequences one by one to reduce memory usage for very large FASTA files"
    )] = False,
    log_level: Annotated[str, typer.Option(
        "--log-level", 
        help="📝 Set logging verbosity level for debugging and monitoring",
        click_type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False)
    )] = "INFO",
    log_file: Annotated[Optional[Path], typer.Option(
        "--log-file", 
        help="📄 Write log messages to specified file instead of console"
    )] = None,
    show_stats: Annotated[bool, typer.Option(
        "--stats", 
        help="📊 Display detailed statistics including Tm ranges, GC content, and product sizes"
    )] = False,
    validate_only: Annotated[bool, typer.Option(
        "--validate-only", 
        help="✅ Only validate FASTA file format and sequence content without primer design"
    )] = False
):
    """
    🧬 Design PCR primers for DNA sequences in a FASTA file.
    
    This command processes all sequences in the input FASTA file and designs
    optimal PCR primer pairs for each sequence using primer3 algorithm.
    
    ⚙️ CONFIGURATION OPTIONS:
    
    \b
    • 🎯 Presets: Use predefined configurations (--config)
    • 📄 Custom: Use your own JSON configuration file (--custom-config)
    
    🎯 PREDEFINED PRESETS:
    
    \b
    • 🎯 standard: General purpose (Tm: 57-63°C, Size: 18-25bp, GC: 40-60%)
    • 🔥 gc-rich: For GC-rich templates >60% (Tm: 59-65°C, Size: 20-27bp, GC: 50-70%)  
    • ❄️  gc-poor: For AT-rich templates <40% (Tm: 55-61°C, Size: 16-22bp, GC: 30-50%)
    
    📄 CUSTOM CONFIGURATION:
    
    \b
    • Export template: primer-designer export-template -o my_config.json
    • Edit parameters: Modify JSON file with your custom values
    • Use in design: primer-designer design input.fasta --custom-config my_config.json
    
    📊 OUTPUT FORMATS:
    
    \b
    • 📋 csv: Comma-separated values with primer sequences, Tm, GC%, positions
    • 📄 json: Structured JSON with metadata and detailed primer information
    
    💡 EXAMPLES:
    
    \b
    Basic usage with preset:
        primer-designer design sequences.fasta
    
    Use custom configuration:
        primer-designer design input.fasta --custom-config my_config.json
    
    Export template and customize:
        primer-designer export-template -o custom.json
        # Edit custom.json with your parameters
        primer-designer design input.fasta --custom-config custom.json
    
    Specify output file and format:
        primer-designer design input.fasta -o results.json -f json
    
    Memory-efficient processing:
        primer-designer design large.fasta --memory-efficient
    """
    # Setup logging
    logger = setup_logging(level=log_level, log_file=log_file)
    
    # Validate input file
    if not input_file.exists():
        console.print(f"[red]❌ Error: Input file not found: {input_file}[/red]")
        console.print("[yellow]💡 Tip: Check the file path and ensure the file exists[/yellow]")
        raise typer.Exit(1)
    
    # Validate file is readable
    if not input_file.is_file():
        console.print(f"[red]❌ Error: Path is not a file: {input_file}[/red]")
        raise typer.Exit(1)
    
    # Validate output format (redundant with Choice but provides better error message)
    if output_format not in ["csv", "json"]:
        console.print(f"[red]❌ Error: Invalid output format '{output_format}'. Must be 'csv' or 'json'[/red]")
        raise typer.Exit(1)
    
    # Validate configuration preset (redundant with Choice but provides better error message)
    if not custom_config and config_preset not in CONFIGS:
        console.print(f"[red]❌ Error: Invalid configuration preset '{config_preset}'[/red]")
        console.print(f"[yellow]📋 Available presets: {', '.join(CONFIGS.keys())}[/yellow]")
        console.print("[blue]💡 Use 'primer-designer list-configs' to see detailed information[/blue]")
        raise typer.Exit(1)
    
    # Validate custom config file if provided
    if custom_config:
        if not custom_config.exists():
            console.print(f"[red]❌ Error: Custom configuration file not found: {custom_config}[/red]")
            console.print("[yellow]💡 Use 'primer-designer export-template' to create a template[/yellow]")
            raise typer.Exit(1)
        
        if not custom_config.is_file():
            console.print(f"[red]❌ Error: Path is not a file: {custom_config}[/red]")
            raise typer.Exit(1)
    
    # Set default output file if not provided
    if output_file is None:
        output_file = input_file.with_suffix(f".primers.{output_format}")
    
    # Warn if format doesn't match file extension
    file_ext = output_file.suffix.lower()
    if output_format == "csv" and file_ext != ".csv":
        console.print(f"[yellow]⚠️ Warning: Output format is CSV but file extension is '{file_ext}'. Consider using .csv extension.[/yellow]")
    elif output_format == "json" and file_ext != ".json":
        console.print(f"[yellow]⚠️ Warning: Output format is JSON but file extension is '{file_ext}'. Consider using .json extension.[/yellow]")
    
    try:
        # Validate FASTA file
        console.print(f"[blue]🔍 Validating FASTA file: {input_file}[/blue]")
        fasta_reader = FastaReader()
        validation_stats = fasta_reader.validate_fasta_file(input_file)
        
        if not validation_stats['valid']:
            console.print(f"[red]❌ Error: {validation_stats['error']}[/red]")
            raise typer.Exit(1)
        
        # Display validation results
        _display_validation_stats(validation_stats)
        
        if validate_only:
            console.print("[green]✅ FASTA file validation completed[/green]")
            return
        
        # Get configuration
        if custom_config:
            try:
                console.print(f"[blue]📄 Loading custom configuration from: {custom_config}[/blue]")
                config = PrimerDesignerConfig.from_json_file(custom_config)
                config_name = f"custom ({custom_config.name})"
            except Exception as e:
                console.print(f"[red]❌ Error loading custom configuration: {e}[/red]")
                console.print("[yellow]💡 Use 'primer-designer export-template' to create a valid template[/yellow]")
                raise typer.Exit(1)
        else:
            config = CONFIGS[config_preset]
            config_name = config_preset
        
        # Display configuration summary
        _display_config_summary(config_name, config)
        
        # Initialize batch designer
        batch_designer = BatchDesigner(
            config=config,
            max_workers=max_workers,
            use_parallel=parallel
        )
        
        # Process sequences
        if memory_efficient:
            primers = _process_memory_efficient(batch_designer, input_file, output_file, output_format)
        else:
            primers = _process_standard(batch_designer, input_file, output_file, output_format)
        
        # Display results
        _display_results(primers, output_file, show_stats, batch_designer)
        
    except Exception as e:
        logger.error(f"Primer design failed: {e}")
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def validate(
    input_file: Annotated[Path, typer.Argument(
        help="📁 Input FASTA file to validate for primer design compatibility (sequences must be ≥50 bp)"
    )],
    log_level: Annotated[str, typer.Option(
        "--log-level", 
        help="📝 Set logging verbosity level for validation process",
        click_type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False)
    )] = "INFO"
):
    """
    ✅ Validate a FASTA file for primer design compatibility.
    
    This command performs comprehensive validation of the input FASTA file:
    
    🔍 VALIDATION CHECKS:
    
    \b
    • 📁 File format and structure
    • 🏷️ Sequence IDs are unique and valid
    • 🧬 DNA sequences contain only valid bases (A, T, G, C, plus ambiguous codes)
    • 📏 Minimum sequence length (≥50 bp required for primer design)
    • 🎯 Sequence quality and composition
    
    📊 VALIDATION RESULTS:
    
    \b
    • 📈 Total number of sequences found
    • ✅❌ Count of valid vs invalid sequences  
    • 📏 Sequence length statistics (min, max, average)
    • 🚨 Detailed error messages for any issues found
    
    💡 EXAMPLES:
    
    \b
    Basic validation:
        primer-designer validate sequences.fasta
    
    Verbose validation with debug output:
        primer-designer validate input.fasta --log-level DEBUG
    """
    # Setup logging
    setup_logging(level=log_level)
    
    if not input_file.exists():
        console.print(f"[red]❌ Error: Input file not found: {input_file}[/red]")
        console.print("[yellow]💡 Tip: Check the file path and ensure the FASTA file exists[/yellow]")
        raise typer.Exit(1)
    
    if not input_file.is_file():
        console.print(f"[red]❌ Error: Path is not a file: {input_file}[/red]")
        raise typer.Exit(1)
    
    try:
        console.print(f"[blue]🔍 Validating FASTA file: {input_file}[/blue]")
        
        fasta_reader = FastaReader()
        validation_stats = fasta_reader.validate_fasta_file(input_file)
        
        if validation_stats['valid']:
            _display_validation_stats(validation_stats)
            console.print("[green]✅ FASTA file is valid for primer design[/green]")
        else:
            console.print(f"[red]❌ FASTA file validation failed: {validation_stats['error']}[/red]")
            raise typer.Exit(1)
            
    except Exception as e:
        console.print(f"[red]❌ Error during validation: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def config_info(
    config_preset: Annotated[str, typer.Argument(
        help="⚙️ Configuration preset name to display detailed information for",
        click_type=click.Choice(["standard", "gc-rich", "gc-poor"], case_sensitive=False)
    )] = "standard"
):
    """
    ⚙️ Display detailed information about a configuration preset.
    
    Shows comprehensive parameter details for the specified configuration preset
    including all primer design constraints and optimization targets.
    
    📋 CONFIGURATION DETAILS SHOWN:
    
    \b
    • 📏 Primer size parameters (optimal, minimum, maximum lengths)
    • 🌡️ Melting temperature ranges and targets (°C)
    • 🧬 GC content constraints and optimal values (%)
    • 🔗 Secondary structure limits (self-complementarity, poly-X repeats)
    • 🧪 Salt and concentration conditions (monovalent, divalent, dNTP, primer)
    • 📊 PCR product size range constraints
    
    🎯 AVAILABLE PRESETS:
    
    \b
    • 🎯 standard: Balanced parameters for general PCR applications
    • 🔥 gc-rich: Optimized for templates with high GC content (>60%)
    • ❄️ gc-poor: Optimized for templates with low GC content (<40%)
    
    💡 EXAMPLES:
    
    \b
    Show standard configuration:
        primer-designer config-info standard
    
    Show GC-rich configuration:
        primer-designer config-info gc-rich
    """
    if config_preset not in CONFIGS:
        console.print(f"[red]❌ Error: Invalid config preset. Available: {', '.join(CONFIGS.keys())}[/red]")
        raise typer.Exit(1)
    
    config = CONFIGS[config_preset]
    _display_detailed_config(config_preset, config)


@app.command()
def export_template(
    output_file: Annotated[Optional[Path], typer.Option(
        "--output", "-o", 
        help="📄 Output JSON file path for the configuration template"
    )] = None,
    preset: Annotated[str, typer.Option(
        "--preset", "-p",
        help="⚙️ Base the template on an existing preset instead of defaults",
        click_type=click.Choice(["standard", "gc-rich", "gc-poor"], case_sensitive=False)
    )] = None,
    show_json: Annotated[bool, typer.Option(
        "--show", 
        help="📺 Display the JSON template in console instead of saving to file"
    )] = False,
    show_params: Annotated[bool, typer.Option(
        "--show-params", 
        help="📋 Display detailed parameter information with types, ranges, and descriptions"
    )] = False,
    documented: Annotated[bool, typer.Option(
        "--documented", 
        help="📚 Include parameter documentation in the exported template"
    )] = False
):
    """
    📄 Export a JSON configuration template for custom primer design parameters.
    
    Creates a JSON template file that can be modified and used with the --custom-config
    option in the design command. The template includes all available parameters with
    their default values and validation constraints.
    
    🎯 TEMPLATE OPTIONS:
    
    \b
    • 📋 Default template: Uses standard default values for all parameters
    • ⚙️ Preset-based: Uses an existing preset as the starting point
    • 📺 Console output: Display template without saving to file
    • 📚 Documented: Include parameter documentation in template
    • 📋 Parameter info: Show detailed parameter specifications
    
    📝 USAGE WORKFLOW:
    
    \b
    1. Export template: primer-designer export-template -o my_config.json
    2. Edit parameters: Modify the JSON file with your custom values
    3. Use in design: primer-designer design input.fasta --custom-config my_config.json
    
    💡 EXAMPLES:
    
    \b
    Export default template:
        primer-designer export-template -o custom_config.json
    
    Show parameter information:
        primer-designer export-template --show-params
    
    Export documented template:
        primer-designer export-template -o config.json --documented
    
    Base template on GC-rich preset:
        primer-designer export-template -o gc_custom.json --preset gc-rich
    
    Show template in console:
        primer-designer export-template --show
    """
    try:
        # Show parameter information if requested
        if show_params:
            _display_parameter_info()
            return
        
        # Get base configuration
        if preset:
            if preset not in CONFIGS:
                console.print(f"[red]❌ Error: Invalid preset '{preset}'[/red]")
                raise typer.Exit(1)
            config = CONFIGS[preset]
            console.print(f"[blue]📋 Using '{preset}' preset as template base[/blue]")
        else:
            config = PrimerDesignerConfig.create_template()
            console.print("[blue]📋 Using default values as template base[/blue]")
        
        # Handle output
        if show_json:
            # Display in console
            console.print("\n[bold cyan]📄 Configuration Template JSON:[/bold cyan]\n")
            if documented:
                documented_template = PrimerDesignerConfig.create_documented_template()
                if preset:
                    # Update with preset values
                    documented_template["config"] = config.model_dump()
                    for field_name, value in config.model_dump().items():
                        if field_name in documented_template["_parameter_info"]:
                            documented_template["_parameter_info"][field_name]["current_value"] = value
                json_str = json.dumps(documented_template, indent=2, ensure_ascii=False)
            else:
                json_str = config.to_json_string(indent=2)
            console.print(json_str)
        else:
            # Save to file
            if output_file is None:
                if documented:
                    output_file = Path("primer_config_documented_template.json")
                else:
                    output_file = Path("primer_config_template.json")
            
            if documented:
                documented_template = PrimerDesignerConfig.create_documented_template()
                if preset:
                    # Update with preset values
                    documented_template["config"] = config.model_dump()
                    for field_name, value in config.model_dump().items():
                        if field_name in documented_template["_parameter_info"]:
                            documented_template["_parameter_info"][field_name]["current_value"] = value
                
                # Save documented template
                output_file.parent.mkdir(parents=True, exist_ok=True)
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(documented_template, f, indent=2, ensure_ascii=False)
            else:
                config.to_json_file(output_file)
            
            console.print(f"[green]✅ Configuration template exported to: {output_file}[/green]")
            if documented:
                console.print("[yellow]📚 Template includes parameter documentation and constraints[/yellow]")
            console.print(f"[yellow]💡 Edit the file and use it with: --custom-config {output_file}[/yellow]")
            
    except Exception as e:
        console.print(f"[red]❌ Error exporting template: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def list_configs():
    """
    📋 List all available configuration presets with summary information.
    
    Displays a comprehensive table showing all predefined configuration presets
    with their key characteristics and recommended use cases.
    
    📊 INFORMATION DISPLAYED:
    
    \b
    • 🏷️ Preset name and description
    • 📏 Primer size range (minimum-maximum bp)
    • 🌡️ Melting temperature range (°C)  
    • 🧬 GC content range (%)
    • 🎯 Recommended template types
    
    🎯 USE CASES:
    
    \b
    • 🎯 standard: Most PCR applications, balanced GC content templates
    • 🔥 gc-rich: High GC content templates (>60%), difficult amplification
    • ❄️ gc-poor: Low GC content templates (<40%), AT-rich sequences
    
    💡 EXAMPLE:
    
    \b
        primer-designer list-configs
    """
    table = Table(title="Available Configuration Presets")
    table.add_column("Preset", style="cyan", no_wrap=True)
    table.add_column("Description", style="magenta")
    table.add_column("Primer Size", style="green")
    table.add_column("Tm Range", style="yellow")
    table.add_column("GC Range", style="blue")
    
    descriptions = {
        "standard": "General purpose PCR primers",
        "gc-rich": "Optimized for GC-rich templates",
        "gc-poor": "Optimized for AT-rich templates"
    }
    
    for preset_name, config in CONFIGS.items():
        table.add_row(
            preset_name,
            descriptions.get(preset_name, ""),
            f"{config.primer_min_size}-{config.primer_max_size} bp",
            f"{config.primer_min_tm}-{config.primer_max_tm}°C",
            f"{config.primer_min_gc}-{config.primer_max_gc}%"
        )
    
    console.print(table)


def _display_parameter_info():
    """Display detailed parameter information from Pydantic model."""
    console.print("\n[bold cyan]📋 Primer Designer Configuration Parameters[/bold cyan]\n")
    
    field_info = PrimerDesignerConfig.get_field_info()
    
    # Group parameters by category
    categories = {
        "🧬 Primer Size Parameters": [
            "primer_opt_size", "primer_min_size", "primer_max_size"
        ],
        "🌡️ Melting Temperature Parameters": [
            "primer_opt_tm", "primer_min_tm", "primer_max_tm"
        ],
        "🔬 GC Content Parameters": [
            "primer_opt_gc_percent", "primer_min_gc", "primer_max_gc"
        ],
        "🔗 Secondary Structure Parameters": [
            "primer_max_poly_x", "primer_max_self_any", "primer_max_self_end"
        ],
        "🧪 Salt and Concentration Parameters": [
            "primer_salt_monovalent", "primer_salt_divalent", 
            "primer_dntp_conc", "primer_dna_conc"
        ],
        "📏 Product Size Parameters": [
            "product_size_range"
        ]
    }
    
    for category, param_names in categories.items():
        console.print(f"[bold yellow]{category}[/bold yellow]")
        
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("Parameter", style="cyan", width=25)
        table.add_column("Type", style="green", width=15)
        table.add_column("Default", style="yellow", width=12)
        table.add_column("Range/Constraints", style="blue", width=20)
        table.add_column("Description", style="white", width=40)
        
        for param_name in param_names:
            if param_name in field_info:
                info = field_info[param_name]
                
                # Format type
                param_type = info["type"].replace("typing.", "").replace("Union[", "").replace(", NoneType]", " | None")
                
                # Format default value
                default_val = str(info["default"])
                if param_name == "product_size_range":
                    default_val = f"[{info['default'][0]}, {info['default'][1]}]"
                
                # Format constraints
                constraints_str = ""
                if info["constraints"]:
                    constraints = info["constraints"]
                    if "range" in constraints:
                        constraints_str = constraints["range"]
                    elif "format" in constraints:
                        constraints_str = constraints["format"]
                    
                    if "unit" in constraints:
                        if constraints_str:
                            constraints_str += f" {constraints['unit']}"
                        else:
                            constraints_str = constraints["unit"]
                
                # Format description
                description = info["description"]
                if len(description) > 35:
                    description = description[:32] + "..."
                
                table.add_row(
                    param_name,
                    param_type,
                    default_val,
                    constraints_str,
                    description
                )
        
        console.print(table)
        console.print()  # Add spacing between categories
    
    # Add usage examples
    console.print("[bold cyan]💡 Usage Examples:[/bold cyan]")
    console.print()
    console.print("[yellow]Export template with documentation:[/yellow]")
    console.print("  primer-designer export-template -o config.json --documented")
    console.print()
    console.print("[yellow]Create custom configuration:[/yellow]")
    console.print("  1. primer-designer export-template -o my_config.json")
    console.print("  2. Edit my_config.json with your desired values")
    console.print("  3. primer-designer design input.fasta --custom-config my_config.json")
    console.print()


def _process_standard(batch_designer, input_file, output_file, output_format):
    """Process sequences using standard method."""
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("🧬 Designing primers...", total=None)
        
        primers = batch_designer.design_primers_from_fasta(
            input_file, output_file, output_format
        )
        
        progress.update(task, completed=True)
    
    return primers


def _process_memory_efficient(batch_designer, input_file, output_file, output_format):
    """Process sequences using memory-efficient method."""
    primers = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("💾 Processing sequences (memory-efficient)...", total=None)
        
        for primer in batch_designer.design_primers_from_fasta_generator(
            input_file, output_file, output_format
        ):
            primers.append(primer)
        
        progress.update(task, completed=True)
    
    return primers


def _display_validation_stats(stats):
    """Display FASTA file validation statistics."""
    table = Table(title="FASTA File Validation Results")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="magenta")
    
    table.add_row("Total sequences", str(stats['total_sequences']))
    table.add_row("Valid sequences", str(stats['valid_sequences']))
    table.add_row("Invalid sequences", str(stats['invalid_sequences']))
    
    if stats['sequence_lengths']:
        table.add_row("Min sequence length", f"{min(stats['sequence_lengths'])} bp")
        table.add_row("Max sequence length", f"{max(stats['sequence_lengths'])} bp")
        table.add_row("Avg sequence length", f"{sum(stats['sequence_lengths']) // len(stats['sequence_lengths'])} bp")
    
    console.print(table)


def _display_config_summary(preset_name, config):
    """Display configuration summary."""
    table = Table(title=f"Configuration: {preset_name}")
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="magenta")
    
    summary = {
        "Primer size range": f"{config.primer_min_size}-{config.primer_max_size} bp (opt: {config.primer_opt_size})",
        "Tm range": f"{config.primer_min_tm}-{config.primer_max_tm}°C (opt: {config.primer_opt_tm})",
        "GC content range": f"{config.primer_min_gc}-{config.primer_max_gc}%",
        "Product size range": f"{config.product_size_range[0]}-{config.product_size_range[1]} bp",
        "Salt concentration": f"{config.primer_salt_monovalent} mM",
        "Primer concentration": f"{config.primer_dna_conc} nM"
    }
    
    for param, value in summary.items():
        table.add_row(param, value)
    
    console.print(table)


def _display_detailed_config(preset_name, config):
    """Display detailed configuration information."""
    console.print(f"\n[bold cyan]Configuration Preset: {preset_name}[/bold cyan]\n")
    
    # Size parameters
    size_table = Table(title="Primer Size Parameters")
    size_table.add_column("Parameter", style="cyan")
    size_table.add_column("Value", style="magenta")
    
    size_table.add_row("Optimal size", f"{config.primer_opt_size} bp")
    size_table.add_row("Minimum size", f"{config.primer_min_size} bp")
    size_table.add_row("Maximum size", f"{config.primer_max_size} bp")
    
    console.print(size_table)
    
    # Temperature parameters
    tm_table = Table(title="Melting Temperature Parameters")
    tm_table.add_column("Parameter", style="cyan")
    tm_table.add_column("Value", style="magenta")
    
    tm_table.add_row("Optimal Tm", f"{config.primer_opt_tm}°C")
    tm_table.add_row("Minimum Tm", f"{config.primer_min_tm}°C")
    tm_table.add_row("Maximum Tm", f"{config.primer_max_tm}°C")
    
    console.print(tm_table)
    
    # GC content parameters
    gc_table = Table(title="GC Content Parameters")
    gc_table.add_column("Parameter", style="cyan")
    gc_table.add_column("Value", style="magenta")
    
    if config.primer_opt_gc_percent:
        gc_table.add_row("Optimal GC%", f"{config.primer_opt_gc_percent}%")
    gc_table.add_row("Minimum GC%", f"{config.primer_min_gc}%")
    gc_table.add_row("Maximum GC%", f"{config.primer_max_gc}%")
    
    console.print(gc_table)


def _display_results(primers, output_file, show_stats, batch_designer):
    """Display processing results."""
    console.print(f"\n[green]🎉 Primer design completed![/green]")
    console.print(f"[blue]💾 Results saved to: {output_file}[/blue]")
    console.print(f"[yellow]✅ Successfully designed primers for {len(primers)} sequences[/yellow]")
    
    if show_stats and primers:
        stats = batch_designer.get_batch_statistics(primers)
        _display_batch_statistics(stats)


def _display_batch_statistics(stats):
    """Display batch processing statistics."""
    console.print("\n[bold cyan]📊 Batch Statistics[/bold cyan]")
    
    # Tm statistics
    tm_table = Table(title="🌡️ Melting Temperature Statistics")
    tm_table.add_column("Primer", style="cyan")
    tm_table.add_column("Min Tm", style="green")
    tm_table.add_column("Max Tm", style="red")
    tm_table.add_column("Mean Tm", style="yellow")
    
    tm_table.add_row(
        "Forward",
        f"{stats['tm_statistics']['forward']['min']:.1f}°C",
        f"{stats['tm_statistics']['forward']['max']:.1f}°C",
        f"{stats['tm_statistics']['forward']['mean']:.1f}°C"
    )
    tm_table.add_row(
        "Reverse",
        f"{stats['tm_statistics']['reverse']['min']:.1f}°C",
        f"{stats['tm_statistics']['reverse']['max']:.1f}°C",
        f"{stats['tm_statistics']['reverse']['mean']:.1f}°C"
    )
    
    console.print(tm_table)
    
    # Product size statistics
    product_table = Table(title="📏 Product Size Statistics")
    product_table.add_column("Metric", style="cyan")
    product_table.add_column("Value", style="magenta")
    
    product_table.add_row("Minimum", f"{stats['product_size_statistics']['min']} bp")
    product_table.add_row("Maximum", f"{stats['product_size_statistics']['max']} bp")
    product_table.add_row("Average", f"{stats['product_size_statistics']['mean']:.0f} bp")
    
    console.print(product_table)


if __name__ == "__main__":
    app()
