rule cluster_all_species:
    input:
        strained=expand("{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained", 
                       outdir=OUTDIR, species=SPECIES_LIST),
        replib=f"{OUTDIR}/{REPSPEC}.RepeatMasker.lib" if REPSPEC else [],
        custom=CUSTOM_LIB if CUSTOM_LIB else []
    output:
        combined="{outdir}/combinedLibraries/combined_all_species.clstrd.fa",
        clstr="{outdir}/combinedLibraries/combined_all_species.clstrd.fa.clstr",
        combined_fa=temp("{outdir}/combinedLibraries/combined_all_species.fa")
    log:
        "{outdir}/combinedLibraries/cluster_all_species.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 64))  # cd-hit: scales 1-64 threads (single job, runs once)
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
        runtime=480
    params:
        outdir=OUTDIR,
        species=SPECIES_LIST,
        repspec=REPSPEC,
        skip_clustering=config.get("skip_clustering", False),
        cluster_identity=config.get("clustering_identity", 0.8),
        cluster_coverage=config.get("clustering_coverage", 0.8),
        cluster_coverage_long=config.get("clustering_coverage_long", 0.0),
        cluster_length_diff=config.get("clustering_length_diff", 0.5)
    run:
        import os
        shell("mkdir -p {params.outdir}/combinedLibraries >> " + str(log) + " 2>&1")
        
        # Define combined file path
        combined_file = f"{params.outdir}/combinedLibraries/combined_all_species.fa"
        
        # Create combined file with genome-prefixed sequences
        with open(combined_file, 'w') as outf:
            # Add strained sequences with genome prefix
            for i, species in enumerate(params.species):
                with open(input.strained[i], 'r') as inf:
                    for line in inf:
                        if line.startswith('>'):
                            outf.write(f">{species}_{line[1:]}")
                        else:
                            outf.write(line)
            
            # Add RepeatMasker library if present with REPMASKER_{species} prefix
            if input.replib and os.path.exists(input.replib):
                with open(input.replib, 'r') as inf:
                    for line in inf:
                        if line.startswith('>'):
                            outf.write(f">REPMASKER_{params.repspec}_{line[1:]}")
                        else:
                            outf.write(line)
            
            # Add custom library if present with CUSTOM prefix
            if input.custom and os.path.exists(input.custom):
                with open(input.custom, 'r') as inf:
                    for line in inf:
                        if line.startswith('>'):
                            outf.write(f">CUSTOM_{line[1:]}")
                        else:
                            outf.write(line)
        
        # Cluster or skip based on config
        if params.skip_clustering:
            # Just copy/rename the combined file as the output (no clustering)
            shell(f"cp {combined_file} {{output.combined}} >> " + str(log) + " 2>&1")
            # cd-hit-est is not run, so no .clstr is produced.
            # Touch an empty sentinel so Snakemake dependency tracking still works.
            # The saturation_plot script detects an empty file and uses a fallback.
            shell("touch {output.clstr} >> " + str(log) + " 2>&1")
        else:
            # Run cd-hit-est clustering with config parameters
            shell(f"cd-hit-est -d 0 -aS {{params.cluster_coverage}} -aL {{params.cluster_coverage_long}} "
                  f"-c {{params.cluster_identity}} "
                  f"-s {{params.cluster_length_diff}} -G 0 -g 1 -b 500 -r 1 "
                  f"-i {combined_file} -o {{output.combined}} "
                  f"-M {{resources.mem_mb}} -T {{threads}} >> " + str(log) + " 2>&1")
            # combined_fa is declared as temp() — Snakemake will delete it once all
            # downstream rules that use it (e.g. split_chimeras) have completed.
        
        # Ensure file is fully written to disk before dependent jobs start
        shell("sync >> " + str(log) + " 2>&1")
        shell("sleep 2")


# ── Post-clustering chimera detection (only when clustering is active) ────────
if not config.get('skip_clustering', False) and config.get('split_chimeras', False):
    rule split_chimeras:
        input:
            clstr=rules.cluster_all_species.output.clstr,
            clustered_fa=rules.cluster_all_species.output.combined,
            combined_fa=rules.cluster_all_species.output.combined_fa
        output:
            fasta="{outdir}/combinedLibraries/combined_all_species.chimera_split.fa",
            summary="{outdir}/combinedLibraries/chimera_detection_summary.tsv"
        log:
            "{outdir}/combinedLibraries/split_chimeras.log"
        params:
            overlap_min=config.get('chimera_overlap_min', 50),
            min_members=config.get('chimera_min_members', 3),
            min_component_span=config.get('chimera_min_component_span', 0.1)
        script:
            "../scripts/split_chimeras.py"