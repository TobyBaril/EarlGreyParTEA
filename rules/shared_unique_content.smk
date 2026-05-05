# shared_unique_content.smk — shared/unique TE family and coverage analysis.
#
# Provides four rules, of which at most two will be active per run:
#
#   shared_unique_plot        — cluster-based (full mode), no phylo tree.
#   shared_unique_plot_phylo  — cluster-based (full mode), with phylo tree
#                               (only when run_busco_phylo is also true).
#   shared_unique_pa_plot        — presence/absence (annotate mode), no tree.
#   shared_unique_pa_plot_phylo  — presence/absence (annotate mode) + tree.
#
# The Snakefile selects the correct pair by conditioning on PIPELINE_MODE and
# RUN_BUSCO_PHYLO before calling include:.
#
# Inputs that must exist before these rules run
# -----------------------------------------------
# Cluster mode only:
#   {OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa.clstr
# All modes:
#   {outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff
#   {outdir}/{species}_EarlGrey/{species}.prep
#
# Outputs (written to {OUTDIR}/sharedUniqueContent/)
# ---------------------------------------------------
#   shared_unique_families.pdf / .tsv
#   shared_unique_coverage.pdf / .tsv
#   (+ phylo variants when a tree is available)

if PIPELINE_MODE == "full" and not RUN_BUSCO_PHYLO:
    rule shared_unique_plot:
        """Cluster-based shared/unique analysis (full pipeline mode, no phylo tree)."""
        input:
            clstr=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa.clstr",
            gffs=expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
            preps=expand(
                "{outdir}/{species}_EarlGrey/{species}.prep",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
        output:
            fam_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.pdf",
            fam_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.tsv",
            cov_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.pdf",
            cov_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
        log:
            f"{OUTDIR}/sharedUniqueContent/shared_unique_plot.log"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 8000 * attempt,
            runtime=60,
        params:
            species=SPECIES_LIST,
            detection_mode="cluster",
            has_phylo_tree=False,
            outdir=OUTDIR,
        script:
            "../scripts/shared_unique_plot.py"

elif PIPELINE_MODE == "full" and RUN_BUSCO_PHYLO:
    rule shared_unique_plot_phylo:
        """Cluster-based shared/unique analysis with phylogenetic ordering."""
        input:
            clstr=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa.clstr",
            gffs=expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
            preps=expand(
                "{outdir}/{species}_EarlGrey/{species}.prep",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
            tree=f"{OUTDIR}/buscoPhylo/species_tree.nwk",
        output:
            fam_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.pdf",
            fam_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.tsv",
            cov_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.pdf",
            cov_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
            fam_phylo_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_families_phylo.pdf",
            cov_phylo_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage_phylo.pdf",
        log:
            f"{OUTDIR}/sharedUniqueContent/shared_unique_plot.log"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 8000 * attempt,
            runtime=60,
        params:
            species=SPECIES_LIST,
            detection_mode="cluster",
            has_phylo_tree=True,
            outdir=OUTDIR,
        script:
            "../scripts/shared_unique_plot.py"

elif PIPELINE_MODE == "annotate" and not RUN_BUSCO_PHYLO:
    rule shared_unique_pa_plot:
        """Presence/absence shared/unique analysis (annotate pipeline mode, no phylo tree)."""
        input:
            gffs=expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
            preps=expand(
                "{outdir}/{species}_EarlGrey/{species}.prep",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
        output:
            fam_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.pdf",
            fam_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.tsv",
            cov_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.pdf",
            cov_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
        log:
            f"{OUTDIR}/sharedUniqueContent/shared_unique_plot.log"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 8000 * attempt,
            runtime=60,
        params:
            species=SPECIES_LIST,
            detection_mode="presence_absence",
            has_phylo_tree=False,
            outdir=OUTDIR,
        script:
            "../scripts/shared_unique_plot.py"

elif PIPELINE_MODE == "annotate" and RUN_BUSCO_PHYLO:
    rule shared_unique_pa_plot_phylo:
        """Presence/absence shared/unique analysis with phylogenetic ordering."""
        input:
            gffs=expand(
                "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
            preps=expand(
                "{outdir}/{species}_EarlGrey/{species}.prep",
                outdir=OUTDIR, species=SPECIES_LIST,
            ),
            tree=f"{OUTDIR}/buscoPhylo/species_tree.nwk",
        output:
            fam_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.pdf",
            fam_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_families.tsv",
            cov_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.pdf",
            cov_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
            fam_phylo_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_families_phylo.pdf",
            cov_phylo_pdf=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage_phylo.pdf",
        log:
            f"{OUTDIR}/sharedUniqueContent/shared_unique_plot.log"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 8000 * attempt,
            runtime=60,
        params:
            species=SPECIES_LIST,
            detection_mode="presence_absence",
            has_phylo_tree=True,
            outdir=OUTDIR,
        script:
            "../scripts/shared_unique_plot.py"
