#!/usr/bin/env python3
"""
generate_config.py — helper for generating EarlGreyParTEA config files
from a genome directory (--genome-dir) or a CSV file (--from-csv).

Called by the earlGreyParTEA* bash entry points when the user passes
--genome-dir or --from-csv alongside --generate-config.
"""

import argparse
import csv
import os
import re
import sys

FASTA_EXTENSIONS = ('.fa.gz', '.fasta.gz', '.fna.gz', '.fa', '.fasta', '.fna')


def sanitize_name(name):
    """Replace characters unsafe in YAML keys or species prefixes with '_'."""
    sanitized = re.sub(r'[^A-Za-z0-9_]', '_', name)
    # A leading digit is valid in YAML only when quoted; avoid it.
    if sanitized and sanitized[0].isdigit():
        sanitized = 's_' + sanitized
    return sanitized


def strip_fasta_ext(filename):
    for ext in FASTA_EXTENSIONS:
        if filename.endswith(ext):
            return filename[:-len(ext)]
    return filename


def genomes_from_dir(genome_dir):
    """Return list of (species_name, abs_path) from a genome directory."""
    genome_dir = os.path.abspath(genome_dir)
    if not os.path.isdir(genome_dir):
        sys.exit(f"[ERROR] --genome-dir '{genome_dir}' is not a directory or does not exist")

    raw_entries = [
        fname for fname in sorted(os.listdir(genome_dir))
        if any(fname.endswith(ext) for ext in FASTA_EXTENSIONS)
    ]
    if not raw_entries:
        sys.exit(
            f"[ERROR] No FASTA files found in '{genome_dir}'.\n"
            f"        Expected extensions: {', '.join(FASTA_EXTENSIONS)}"
        )

    entries = []
    seen = {}
    for fname in raw_entries:
        raw_name = strip_fasta_ext(fname)
        species = sanitize_name(raw_name)
        if species != raw_name:
            print(f"[WARN] Species name '{raw_name}' sanitized to '{species}'", file=sys.stderr)
        # Deduplicate after sanitization
        if species in seen:
            seen[species] += 1
            new_name = f"{species}_{seen[species]}"
            print(
                f"[WARN] Duplicate species name '{species}' after sanitization; "
                f"renamed to '{new_name}'",
                file=sys.stderr,
            )
            species = new_name
        else:
            seen[species] = 1
        entries.append((species, os.path.join(genome_dir, fname)))
    return entries


def genomes_from_csv(csv_path):
    """Return list of (species_name, abs_path) from a CSV file."""
    csv_path = os.path.abspath(csv_path)
    if not os.path.isfile(csv_path):
        sys.exit(f"[ERROR] --from-csv '{csv_path}' does not exist")

    csv_dir = os.path.dirname(csv_path)

    with open(csv_path, newline='') as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            sys.exit("[ERROR] CSV file is empty or has no header row")

        # Case-insensitive column lookup with accepted aliases
        fl = {f.strip().lower(): f for f in reader.fieldnames}
        species_col = next(
            (fl[k] for k in ('species',) if k in fl), None
        )
        path_col = next(
            (fl[k] for k in ('genome_path', 'genome', 'path') if k in fl), None
        )
        if species_col is None:
            sys.exit(
                f"[ERROR] CSV must contain a 'species' column.\n"
                f"        Columns found: {list(reader.fieldnames)}"
            )
        if path_col is None:
            sys.exit(
                f"[ERROR] CSV must contain a 'genome_path' column "
                f"(also accepted: 'genome', 'path').\n"
                f"        Columns found: {list(reader.fieldnames)}"
            )

        entries = []
        for i, row in enumerate(reader, start=2):
            raw_name = row[species_col].strip()
            raw_path = row[path_col].strip()
            if not raw_name or not raw_path:
                print(f"[WARN] Skipping row {i}: empty species or genome_path", file=sys.stderr)
                continue
            species = sanitize_name(raw_name)
            if species != raw_name:
                print(f"[WARN] Species name '{raw_name}' sanitized to '{species}'", file=sys.stderr)
            # Resolve relative paths relative to the CSV file's location
            if not os.path.isabs(raw_path):
                raw_path = os.path.join(csv_dir, raw_path)
            entries.append((species, os.path.abspath(raw_path)))

    if not entries:
        sys.exit("[ERROR] No valid entries found in CSV file")
    return entries


def build_genome_block(entries):
    lines = ['genome:']
    for name, path in entries:
        lines.append(f'  {name}: {path}')
    lines.append('')
    lines.append('species:')
    for name, _ in entries:
        lines.append(f'  - {name}')
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# YAML config templates — mirror the inline heredocs in the bash scripts but
# with genome/species/output_dir filled in dynamically.
# ---------------------------------------------------------------------------

FULL_TEMPLATE = """\
# EarlGrey Pangenome Pipeline Configuration
# Full Mode: Complete library construction and annotation

# Input genomes
{genome_block}

# Output directory for all results
output_dir: {output_dir}

# Pipeline mode (do not change for earlGreyParTEA)
pipeline_mode: "full"

# Initial masking with known repeats (choose ONE or leave both empty)
repeatmasker_species: ""  # e.g., "fungi", "arthropoda", "viridiplantae"
custom_library: ""        # path to custom TE library in fasta format

# Library construction parameters
iterations: 10           # Number of BLAST-extend-align cycles
flank: 1000             # Flanking basepairs to extract
max_consensus_seqs: 20  # Max sequences for consensus building
min_consensus_seqs: 3   # Min sequences required for consensus

# Clustering options (for combining TE libraries from multiple genomes)
skip_clustering: false          # Set to true to skip clustering (just concatenate)
clustering_identity: 0.8        # cd-hit sequence identity threshold (0.0-1.0)
clustering_coverage: 0.8        # cd-hit alignment coverage (0.0-1.0)

# Output options
softmask: false         # Generate softmasked genome for each input
margin: false           # Remove short TE sequences (<100bp)
run_heliano: true       # Run HELIANO for Helitron detection

# Workflow visualization (requires graphviz/dot installed)
generate_dag: true      # Generate workflow DAG visualizations
dag_format: "svg"       # Options: "svg", "png", "pdf"

# Saturation analysis options
saturation_permutations: 100  # Number of random genome-addition permutations to average
                              # over for the TE family saturation plot. Higher values
                              # give smoother confidence intervals but increase runtime.

# ---------------------------------------------------------------------------
# Optional analysis modules (v0.1.6+)
# ---------------------------------------------------------------------------

# Shared/unique TE content analysis (full and annotate modes only)
run_shared_unique: false

# BUSCO-based phylogenomics
run_busco_phylo: false
busco_lineage: ""        # REQUIRED if run_busco_phylo: true  (e.g. "fungi_odb10")
busco_prefix: "busco"   # Prefix for BUSCO run directory names
busco_min_occupancy: 0.5  # Min fraction of species a gene must appear in (0.0-1.0)

# Advanced options (usually not needed)
# script_dir: "/path/to/earlgrey/scripts"  # Auto-detected if installed via conda/mamba

# SLURM cluster settings (only used with --slurm flag)
slurm_partition: ""   # partition/queue to submit to (required when using --slurm; can be set here instead of --slurm-partition)
slurm_account: ""     # account string (leave empty if not required)
slurm_extra: ""       # any extra sbatch flags, e.g. "--constraint=avx2"
"""

LIBCONSTRUCT_TEMPLATE = """\
# EarlGrey Pangenome Pipeline Configuration
# LibConstruct Mode: Library construction only (stops after clustering)

# Input genomes
{genome_block}

# Output directory for all results
output_dir: {output_dir}

# Pipeline mode (do not change for earlGreyParTEA_LibConstruct)
pipeline_mode: "libconstruct"

# Initial masking with known repeats (choose ONE or leave both empty)
repeatmasker_species: ""  # e.g., "fungi", "arthropoda", "viridiplantae"
custom_library: ""        # path to custom TE library in fasta format

# Library construction parameters
iterations: 10           # Number of BLAST-extend-align cycles
flank: 1000             # Flanking basepairs to extract
max_consensus_seqs: 20  # Max sequences for consensus building
min_consensus_seqs: 3   # Min sequences required for consensus

# Clustering options (for combining TE libraries from multiple genomes)
skip_clustering: false          # Set to true to skip clustering (just concatenate)
clustering_identity: 0.8        # cd-hit sequence identity threshold (0.0-1.0)
clustering_coverage: 0.8        # cd-hit alignment coverage (0.0-1.0)

# Output options (not applicable in libconstruct mode)
softmask: false
margin: false
run_heliano: false

# Workflow visualization (requires graphviz/dot installed)
generate_dag: true      # Generate workflow DAG visualizations
dag_format: "svg"       # Options: "svg", "png", "pdf"

# Saturation analysis options
saturation_permutations: 100  # Number of random genome-addition permutations to average
                              # over for the TE family saturation plot. Higher values
                              # give smoother confidence intervals but increase runtime.

# ---------------------------------------------------------------------------
# Optional analysis modules (v0.1.6+)
# ---------------------------------------------------------------------------

# Shared/unique TE content: not applicable in libconstruct mode (no annotation outputs)
run_shared_unique: false

# BUSCO-based phylogenomics
run_busco_phylo: false
busco_lineage: ""        # REQUIRED if run_busco_phylo: true  (e.g. "fungi_odb10")
busco_prefix: "busco"   # Prefix for BUSCO run directory names
busco_min_occupancy: 0.5

# Advanced options (usually not needed)
# script_dir: "/path/to/earlgrey/scripts"  # Auto-detected if installed via conda/mamba

# SLURM cluster settings (only used with --slurm flag)
slurm_partition: ""   # partition/queue to submit to (required when using --slurm; can be set here instead of --slurm-partition)
slurm_account: ""     # account string (leave empty if not required)
slurm_extra: ""       # any extra sbatch flags, e.g. "--constraint=avx2"
"""

ANNOTATE_TEMPLATE = """\
# EarlGrey Pangenome Pipeline Configuration
# AnnotationOnly Mode: Annotation with pre-existing TE library

# Input genomes
{genome_block}

# Output directory for all results
output_dir: {output_dir}

# Pipeline mode (do not change for earlGreyParTEA_AnnotationOnly)
pipeline_mode: "annotate"

# REQUIRED: Path to pre-existing TE library
annotation_library: "/path/to/your/TE_library.fasta"

# Library construction parameters (not used in annotation mode)
repeatmasker_species: ""
custom_library: ""
iterations: 10
flank: 1000
max_consensus_seqs: 20
min_consensus_seqs: 3
skip_clustering: false
clustering_identity: 0.8
clustering_coverage: 0.8

# Output options
softmask: true          # Generate softmasked genome for each input
margin: false           # Remove short TE sequences (<100bp)
run_heliano: true       # Run HELIANO for Helitron detection

# Workflow visualization (requires graphviz/dot installed)
generate_dag: true      # Generate workflow DAG visualizations
dag_format: "svg"       # Options: "svg", "png", "pdf"

# ---------------------------------------------------------------------------
# Optional analysis modules (v0.1.6+)
# ---------------------------------------------------------------------------

# Shared/unique TE content analysis: requires GFF outputs from annotation
run_shared_unique: false

# BUSCO-based phylogenomics
run_busco_phylo: false
busco_lineage: ""        # REQUIRED if run_busco_phylo: true  (e.g. "fungi_odb10")
busco_prefix: "busco"   # Prefix for BUSCO run directory names
busco_min_occupancy: 0.5  # Min fraction of species a gene must appear in (0.0-1.0)

# Advanced options (usually not needed)
# script_dir: "/path/to/earlgrey/scripts"  # Auto-detected if installed via conda/mamba

# SLURM cluster settings (only used with --slurm flag)
slurm_partition: ""   # partition/queue to submit to (required when using --slurm; can be set here instead of --slurm-partition)
slurm_account: ""     # account string (leave empty if not required)
slurm_extra: ""       # any extra sbatch flags, e.g. "--constraint=avx2"
"""

TEMPLATES = {
    'full': FULL_TEMPLATE,
    'libconstruct': LIBCONSTRUCT_TEMPLATE,
    'annotate': ANNOTATE_TEMPLATE,
}


def main():
    parser = argparse.ArgumentParser(
        description='Generate an EarlGreyParTEA config file from a genome directory or CSV.'
    )
    parser.add_argument('--output', required=True, help='Output config YAML filename')
    parser.add_argument(
        '--mode', required=True, choices=['full', 'libconstruct', 'annotate'],
        help='Pipeline mode this config is for'
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--genome-dir', metavar='DIR',
                        help='Directory containing FASTA genome files')
    source.add_argument('--from-csv', metavar='FILE',
                        help='CSV file with species and genome_path columns')
    parser.add_argument('--output-dir', metavar='DIR', default=None,
                        help='Value to write as output_dir in the config')
    args = parser.parse_args()

    if args.genome_dir:
        entries = genomes_from_dir(args.genome_dir)
    else:
        entries = genomes_from_csv(args.from_csv)

    genome_block = build_genome_block(entries)
    out_dir = os.path.abspath(args.output_dir) if args.output_dir else '/path/to/output/directory'

    content = TEMPLATES[args.mode].format(genome_block=genome_block, output_dir=out_dir)

    with open(args.output, 'w') as fh:
        fh.write(content)

    print(f"[INFO] Config file generated: {args.output}")
    print(f"[INFO] {len(entries)} genome(s) added to config")
    if args.output_dir is None:
        print("[INFO] Remember to set 'output_dir' in the config before running")
    if args.mode == 'annotate':
        print("[INFO] IMPORTANT: Set 'annotation_library' to your TE library path before running")


if __name__ == '__main__':
    main()
