# EarlGrey Pangenome Pipeline (ParTEA)

> **Repository:** https://github.com/TobyBaril/EarlGreyParTEA  
> **Depends on:** [EarlGrey](https://github.com/TobyBaril/EarlGrey) ≥7.0.3

A Snakemake-based multi-genome transposable element analysis pipeline for building pangenome TE libraries and performing comparative TE annotation across multiple genomes.

## Installation

### Via conda/mamba (recommended)

```bash
# Install earlgrey-partea (automatically installs earlgrey as dependency)
mamba install -c conda-forge -c bioconda earlgrey-partea

# Verify installation
earlGreyParTEA --help
```

**Version Compatibility**: ParTEA automatically detects installed EarlGrey versions (7.x, 8.x, or higher) and requires no configuration changes when EarlGrey is updated.

### Development Installation

```bash
git clone https://github.com/TobyBaril/EarlGreyParTEA.git
cd EarlGreyParTEA
chmod +x earlGreyParTEA*
export PATH="$PWD:$PATH"
```

## Quick Start

### 1. Generate a config file

```bash
earlGreyParTEA --generate-config my_config.yaml
```

### 2. Edit the config file with your genome paths

```yaml
genome:
  species1: /path/to/genome1.fasta
  species2: /path/to/genome2.fasta
  species3: /path/to/genome3.fasta

species:
  - species1
  - species2
  - species3

output_dir: /path/to/output
```

### 3. Run the pipeline

```bash
earlGreyParTEA -c my_config.yaml -t 16
```

## Pipeline Modes

### Full Pipeline (`earlGreyParTEA`)

Runs the complete analysis: library construction → clustering → annotation

```bash
earlGreyParTEA -c config.yaml -t 16
```

**Output:**
- Pangenome TE library (clustered across all genomes)
- TE annotations for each genome (BED, GFF)
- Divergence analysis
- Summary charts and statistics

### Library Construction Only (`earlGreyParTEA_LibConstruct`)

Builds a pangenome TE library without performing annotation.

```bash
earlGreyParTEA_LibConstruct -c config.yaml -t 16
```

**Output:**
- `{output_dir}/combinedLibraries/combined_all_species.clstrd.fa`

**Use case:** Build a TE library from multiple genomes to use for annotating other genomes.

### Annotation Only (`earlGreyParTEA_AnnotationOnly`)

Annotates genomes using a pre-existing TE library (skips library construction).

```bash
earlGreyParTEA_AnnotationOnly -c config.yaml -t 16
```

**Requirements:**
- Must specify `annotation_library` in config.yaml
- Library should be in fasta format

**Use case:** Annotate multiple genomes with a curated TE library from a previous run or external source.

## Command-Line Options

| Option | Short | Description |
|--------|-------|-------------|
| `--config FILE` | `-c` | Config file (required) |
| `--threads INT` | `-t` | Number of threads (required) |
| `--memory INT` | `-m` | Max memory in MB (optional) |
| `--dry-run` | `-n` | Show what would run without executing |
| `--generate-config FILE` | - | Generate example config template |
| `--unlock` | - | Unlock directory after crash |
| `--rerun-incomplete` | - | Rerun incomplete jobs |
| `--help` | `-h` | Show help message |

## Configuration Parameters

### Required Parameters

```yaml
genome:                    # Dictionary of genome paths
  species1: /path/to/genome1.fasta
  
species: [species1]       # List of species to analyze

output_dir: /path/to/out  # Output directory
```

**Note:** The EarlGrey `script_dir` parameter is **automatically detected** and does not need to be specified in your config file. ParTEA will find the correct EarlGrey installation regardless of version (7.x, 8.x, etc.). Only set `script_dir` manually if you have a custom installation location.

### Library Construction Parameters

```yaml
iterations: 10            # BLAST-extend-align cycles
flank: 1000              # Flanking basepairs to extract
max_consensus_seqs: 20   # Max sequences for consensus
min_consensus_seqs: 3    # Min sequences for consensus
```

### Initial Masking (Optional)

Choose ONE or leave both empty:

```yaml
repeatmasker_species: "fungi"        # Use RepeatMasker database
# OR
custom_library: "/path/to/lib.fa"   # Use custom library
```

### Clustering Options

```yaml
skip_clustering: false     # Set true to skip clustering
clustering_identity: 0.8   # cd-hit identity threshold (0.0-1.0)
clustering_coverage: 0.8   # cd-hit coverage threshold (0.0-1.0)
```

### Output Options

```yaml
softmask: false           # Generate softmasked genomes
margin: false             # Remove short TEs (<100bp)
run_heliano: true         # Run HELIANO for Helitron detection
```

### Visualization Options

```yaml
generate_dag: true        # Generate workflow DAG graphs
dag_format: "svg"         # Format: svg, png, or pdf
```

## Output Structure

```
output_dir/
├── combinedLibraries/
│   ├── combined_all_species.clstrd.fa      # Pangenome TE library
│   └── combined_all_species.nonclstrd.fa   # Unclustered library
│
├── species1_EarlGrey/
│   ├── species1_Database/              # RepeatModeler database
│   ├── species1_RepeatModeler/         # RepeatModeler working files
│   ├── species1_strainer/              # TEstrainer output
│   ├── species1_RepeatMasker_Against_Custom_Library/
│   ├── species1_mergedRepeats/         # Merged annotations
│   └── species1_summaryFiles/          # Final outputs
│       ├── species1.filteredRepeats.bed
│       ├── species1.filteredRepeats.gff
│       ├── species1.highLevelCount.txt
│       ├── species1.summaryPie.pdf
│       ├── species1_divergence_summary_table.tsv
│       └── species1.softmasked.fasta (if enabled)
│
├── species2_EarlGrey/
│   └── ...
│
├── workflow_visualization/
│   ├── dag_full_mode.svg               # Workflow DAG visualization
│   └── dag_full_mode_rulegraph.svg     # Simplified rule graph
│
└── validated_config.yaml               # Config used for run
```

## Example Workflows

### Example 1: Full Analysis of Multiple Genomes

```bash
# Generate config
earlGreyParTEA --generate-config analysis.yaml

# Edit config with genome paths
# Then run
earlGreyParTEA -c analysis.yaml -t 32 -m 128000
```

### Example 2: Build Pangenome Library

```bash
# Generate config for library construction
earlGreyParTEA_LibConstruct --generate-config build_lib.yaml

# Edit config, then build library
earlGreyParTEA_LibConstruct -c build_lib.yaml -t 16

# Output: build_lib_output/combinedLibraries/combined_all_species.clstrd.fa
```

### Example 3: Annotate with Pre-existing Library

```bash
# Generate config for annotation
earlGreyParTEA_AnnotationOnly --generate-config annotate.yaml

# Edit config and set annotation_library parameter
# annotation_library: "/path/to/combined_all_species.clstrd.fa"

# Run annotation
earlGreyParTEA_AnnotationOnly -c annotate.yaml -t 16
```

### Example 4: Dry Run to Check Pipeline

```bash
earlGreyParTEA -c config.yaml -t 16 --dry-run
```

## Dynamic Resource Allocation

The pipeline automatically distributes threads across genomes:

| Cores | Genomes | Threads/Genome | Parallel Jobs |
|-------|---------|----------------|---------------|
| 8     | 2       | 4              | 2 genomes     |
| 16    | 4       | 4              | 4 genomes     |
| 32    | 2       | 16             | 2 genomes     |
| 64    | 8       | 8              | 8 genomes     |

The pipeline maximizes efficiency by:
- Running multiple genomes in parallel when cores available
- Using fewer threads per genome when many genomes analyzed
- Capping threads at optimal levels for each tool

## Requirements

- Snakemake ≥7.0
- Python ≥3.9
- EarlGrey dependencies (installed automatically with conda/mamba)
- Graphviz (optional, for DAG visualization)

## Troubleshooting

### Error: "Config file required"

Make sure you specify the config file:
```bash
earlGreyParTEA -c config.yaml -t 16
```

### Error: "Both RepeatMasker species and custom library specified"

Choose only ONE initial masking method in your config:
```yaml
# Either
repeatmasker_species: "fungi"
custom_library: ""

# Or
repeatmasker_species: ""
custom_library: "/path/to/library.fa"
```

### Error: "Pipeline mode 'annotate' requires 'annotation_library'"

For annotation-only mode, you must specify a TE library:
```yaml
pipeline_mode: "annotate"
annotation_library: "/path/to/TE_library.fasta"
```

### Pipeline stops early or has incomplete output

Try rerunning incomplete jobs:
```bash
earlGreyParTEA -c config.yaml -t 16 --rerun-incomplete
```

### Snakemake directory locked after crash

Unlock the directory:
```bash
earlGreyParTEA -c config.yaml -t 16 --unlock
```

### Error: "Script directory not found" or "TEstrainer module not found"

This means ParTEA couldn't auto-detect your EarlGrey installation. This usually happens with custom installations. Check your EarlGrey is installed:

```bash
# Check if EarlGrey is available
which earlGrey

# Check conda environment
conda list | grep earlgrey
```

If installed, manually specify the script directory in your config:
```yaml
script_dir: "/path/to/earlgrey/scripts"
```

For conda installations, this is typically:
```yaml
script_dir: "$CONDA_PREFIX/share/earlgrey-7.0.3-0/scripts"  # Adjust version
```

### DAG visualization not generated

Install graphviz:
```bash
mamba install graphviz
```

Or disable DAG generation in config:
```yaml
generate_dag: false
```

## Citation

If you use this pipeline, please cite:

Baril, T., Galbraith, J. and Hayward, A., 2024. Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline. *Molecular Biology and Evolution*, 41(4), p.msae068.

## Support

- GitHub Issues: https://github.com/TobyBaril/EarlGreyParTEA/issues
- Email: tobias.baril[at]unine.ch
- Documentation: https://github.com/TobyBaril/EarlGreyParTEA

## License

See LICENSE file in the EarlGrey repository.
