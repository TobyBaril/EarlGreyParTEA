import sys
import os

sys.path.append(os.path.dirname(workflow.snakefile))
from scripts.on_start_functions import (
    running_tea, 
    validate_parameters,
    make_directories,
)
from scripts.generate_dag import generate_dag


# Apply at parse time
config = validate_parameters(
    config, 
    outfile = os.path.join(config["output_dir"], "validated_config.yaml")
    )

GENOME = config["genome"]
SPECIES_LIST = config["species"]
OUTDIR = config["output_dir"] #os.path.join(config["output_dir"], f"{SPECIES}_EarlGrey")
REPSPEC = config["repeatmasker_species"]
CUSTOM_LIB = config["custom_library"]
ITER = config["iterations"]
FLANK = config["flank"]
MAX_SEQ = config["max_consensus_seqs"]
MIN_SEQ = config["min_consensus_seqs"]
SCRIPT_DIR = config["script_dir"]
HELI = config.get("run_heliano", False)
# Convert boolean to string for shell compatibility
_softmask_val = config.get("softmask", False)
SOFTMASK = "yes" if (_softmask_val is True or _softmask_val == "yes") else "no"
_margin_val = config.get("margin", False)
MARGIN = "yes" if (_margin_val is True or _margin_val == "yes") else "no"

# Get pipeline mode
PIPELINE_MODE = config.get("pipeline_mode", "full")
ANNOTATION_LIB = config.get("annotation_library", "")

# Optional analysis modules
RUN_SHARED_UNIQUE = config.get("run_shared_unique", False)
RUN_BUSCO_PHYLO   = config.get("run_busco_phylo",   False)
BUSCO_LINEAGE     = config.get("busco_lineage",     "")
BUSCO_PREFIX      = config.get("busco_prefix",      "busco")
BUSCO_MIN_OCC     = config.get("busco_min_occupancy", 0.5)

# DAG visualization settings
GENERATE_DAG = config.get("generate_dag", True)  # Generate DAG by default
DAG_FORMAT = config.get("dag_format", "svg")  # svg, png, or pdf

# Include rules based on pipeline mode
# lib_construct.smk always included (for prep_genome rule needed by annotation)
include: "rules/lib_construct.smk"

if PIPELINE_MODE in ["full", "libconstruct"]:
    # Clustering needed for library construction modes
    include: "rules/clustering.smk"
    include: "rules/saturation.smk"

if PIPELINE_MODE in ["full", "annotate"]:
    # Annotation rules needed
    include: "rules/annotate.smk"

# Optional: shared/unique TE content analysis
if RUN_SHARED_UNIQUE and PIPELINE_MODE in ["full", "annotate"]:
    include: "rules/shared_unique_content.smk"

# Optional: BUSCO phylogenomics
if RUN_BUSCO_PHYLO:
    include: "rules/busco_phylo.smk"

onstart:
    # Display tea art
    running_tea("Earl Grey Pipeline Starting")

    # Directories are created by individual rules as needed
    # make_directories is not needed for pangenome pipeline
    
    # Generate DAG visualization
    if GENERATE_DAG:
        dag_dir = os.path.join(OUTDIR, "workflow_visualization")
        generate_dag(
            snakefile=workflow.snakefile,
            configfile=workflow.configfiles[0] if workflow.configfiles else "config/config.yaml",
            output_dir=dag_dir,
            cores=workflow.cores,
            format=DAG_FORMAT,
            filename_prefix=f"dag_{PIPELINE_MODE}_mode"
        )


# Conditional rule all based on pipeline mode
if PIPELINE_MODE == "libconstruct":
    rule all:
        input:
            # Library construction only - request the clustered library
            f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa",
            f"{OUTDIR}/combinedLibraries/saturation_plot.pdf",
            f"{OUTDIR}/combinedLibraries/saturation_data.tsv",
            # shared/unique not applicable for libconstruct (no annotation)

elif PIPELINE_MODE == "annotate":
    rule all:
        input:
            # Annotation only - request final annotation outputs (skips library construction)
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}_divergence_summary_table.tsv",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.softmasked.fasta",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ) if SOFTMASK is True or SOFTMASK == 'yes' else [],
            # Optional: shared/unique (presence/absence mode in annotate)
            *(
                [
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_families.pdf",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.pdf",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_families.tsv",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
                ]
                if RUN_SHARED_UNIQUE else []
            ),
            # Optional: BUSCO phylo
            *(
                [
                    f"{OUTDIR}/buscoPhylo/busco_completeness.pdf",
                    f"{OUTDIR}/buscoPhylo/busco_completeness.tsv",
                    f"{OUTDIR}/buscoPhylo/busco_completeness_phylo.pdf",
                    f"{OUTDIR}/buscoPhylo/species_tree.nwk",
                ]
                if RUN_BUSCO_PHYLO else []
            ),
            # Optional: phylo-ordered shared/unique + QC scatter (both modules)
            *(
                [
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_families_phylo.pdf",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage_phylo.pdf",
                    f"{OUTDIR}/buscoPhylo/busco_te_qc.pdf",
                    f"{OUTDIR}/buscoPhylo/busco_te_qc.tsv",
                ]
                if (RUN_SHARED_UNIQUE and RUN_BUSCO_PHYLO) else []
            )

else:  # "full" mode (default)
    rule all:
        input:
            # Full pipeline - request all final annotation outputs
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}_divergence_summary_table.tsv",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ),
            expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.softmasked.fasta",
                outdir=OUTDIR,
                species=SPECIES_LIST
            ) if SOFTMASK is True or SOFTMASK == 'yes' else [],
            f"{OUTDIR}/combinedLibraries/saturation_plot.pdf",
            f"{OUTDIR}/combinedLibraries/saturation_data.tsv",
            # Optional: cluster-based shared/unique
            *(
                [
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_families.pdf",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.pdf",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_families.tsv",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
                ]
                if RUN_SHARED_UNIQUE else []
            ),
            # Optional: BUSCO phylogenomics
            *(
                [
                    f"{OUTDIR}/buscoPhylo/busco_completeness.pdf",
                    f"{OUTDIR}/buscoPhylo/busco_completeness.tsv",
                    f"{OUTDIR}/buscoPhylo/busco_completeness_phylo.pdf",
                    f"{OUTDIR}/buscoPhylo/species_tree.nwk",
                ]
                if RUN_BUSCO_PHYLO else []
            ),
            # Optional: phylo-ordered shared/unique + QC scatter (both active)
            *(
                [
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_families_phylo.pdf",
                    f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage_phylo.pdf",
                    f"{OUTDIR}/buscoPhylo/busco_te_qc.pdf",
                    f"{OUTDIR}/buscoPhylo/busco_te_qc.tsv",
                ]
                if (RUN_SHARED_UNIQUE and RUN_BUSCO_PHYLO) else []
            )


# Rule to symlink user-provided library for annotation mode
if PIPELINE_MODE == "annotate":
    rule symlink_annotation_library:
        input:
            ANNOTATION_LIB
        output:
            f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa"
        shell:
            """
            mkdir -p {OUTDIR}/combinedLibraries
            ln -sf $(realpath {input}) {output}
            """
