# EarlGrey ParTEA - Technical Documentation

This directory contains technical documentation for developers, maintainers, and those interested in the implementation details of EarlGrey ParTEA.

## Documentation Files

### For Package Maintainers

**[PACKAGING.md](PACKAGING.md)** - Conda/Mamba packaging guide
- Package structure and dependencies
- Build scripts (meta.yaml, build.sh)
- Version-agnostic installation details
- Testing procedures for different EarlGrey versions
- `snakemake-executor-plugin-slurm` run dependency (required for SLURM mode)

### For Developers

**[VERSION_COMPATIBILITY.md](VERSION_COMPATIBILITY.md)** - Version robustness documentation
- How version-agnostic detection works
- Glob pattern usage for finding any EarlGrey version
- Testing results and upgrade scenarios
- Implementation details for wrapper scripts

**[AUTO_DETECTION.md](AUTO_DETECTION.md)** - Auto-detection feature documentation
- How script_dir auto-detection works
- Search path prioritization
- User experience and configuration details
- Troubleshooting auto-detection failures

**[DAG_VISUALIZATION.md](DAG_VISUALIZATION.md)** - Workflow visualization guide
- Understanding generated DAG files
- Viewing and interpreting workflow graphs
- Graphviz installation and configuration
- Customizing DAG output formats

## User-Facing Documentation

For general usage, installation, and getting started, see the main [README.md](../README.md) in the repository root.

New in v0.1.5:
- **SLURM cluster submission** — `--slurm` flag and related options; see the [SLURM Cluster Submission](../README.md#%EF%B8%8F-slurm-cluster-submission) section in the main README.
- **Auto-populate configs** — `--genome-dir` and `--from-csv` flags; see [Command-Line Options](../README.md#command-line-options) and [Example Workflows](../README.md#example-workflows) in the main README.

New in v0.1.7:
- **Per-rule log files** — all major rules write their tool output to a dedicated log file placed alongside the relevant output directory; see [Changes in v0.1.7](../README.md#previous-release-v017) in the main README.
- **RepeatModeler small-genome robustness** — automatic `-genomeSampleSizeMax` capping; see [Troubleshooting](../README.md#troubleshooting).
- **Extended `repeatmasker_warmup`** — species-specific BLAST cache pre-building and repair before parallel jobs start.
- **Removal of `-norna` flag** — RNA-derived repeats are now included in all RepeatMasker annotation runs.

New in v0.1.8:
- **CD-HIT `-s` cutoff** (`clustering_length_diff`) — prevents very short sequences merging with full-length cluster representatives.
- **CD-HIT `-aL` filter** (`clustering_coverage_long`) — optional coverage requirement on the longer sequence in each cluster pair.
- **Post-clustering chimera detection** (`split_chimeras`) — detects and splits chimeric cluster representatives; outputs `chimera_split.fa` and `chimera_detection_summary.tsv`.
- **Dynamic `CONS-*` cache directory discovery** — `repeatmasker_warmup` now works on Dfam-only installations.

For contribution guidelines, see [CONTRIBUTING.md](../CONTRIBUTING.md).

## Quick Links

- **Main Repository:** https://github.com/TobyBaril/EarlGreyParTEA
- **Issue Tracker:** https://github.com/TobyBaril/EarlGreyParTEA/issues
- **EarlGrey (Dependency):** https://github.com/TobyBaril/EarlGrey
