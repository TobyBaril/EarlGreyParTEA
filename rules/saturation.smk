# Saturation plot rule for EarlGreyParTEA
#
# Visualises how many unique TE families accumulate as each additional genome
# is added to the pipeline.  Depends on the cd-hit-est .clstr file produced
# by cluster_all_species so that cluster attribution is derived from the same
# clustering run used for annotation (no re-clustering required).
#
# Outputs (written to {OUTDIR}/combinedLibraries/):
#   saturation_plot.pdf   - saturation curve (mean ± 95% CI across permutations)
#   saturation_data.tsv   - tabular statistics for downstream use

rule saturation_plot:
    input:
        clstr=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa.clstr",
        strained=expand(
            "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained",
            outdir=OUTDIR,
            species=SPECIES_LIST,
        ),
    output:
        plot=f"{OUTDIR}/combinedLibraries/saturation_plot.pdf",
        table=f"{OUTDIR}/combinedLibraries/saturation_data.tsv",
    log:
        f"{OUTDIR}/combinedLibraries/saturation_plot.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    params:
        species=SPECIES_LIST,
        repspec=REPSPEC,
        has_custom=bool(CUSTOM_LIB),
        permutations=lambda wildcards: config.get("saturation_permutations", 100),
        outdir=OUTDIR,
    script:
        "../scripts/saturation_plot.py"
