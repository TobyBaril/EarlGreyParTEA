rule cluster_all_species:
    input:
        strained=expand("{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained", 
                       outdir=OUTDIR, species=SPECIES_LIST),
        replib=f"{OUTDIR}/{REPSPEC}.RepeatMasker.lib" if REPSPEC else [],
        custom=CUSTOM_LIB if CUSTOM_LIB else []
    output:
        combined="{outdir}/combinedLibraries/combined_all_species.clstrd.fa",
        clstr="{outdir}/combinedLibraries/combined_all_species.clstrd.fa.clstr"
    threads: lambda wildcards: max(1, min(workflow.cores, 32))  # cd-hit: scales 1-32 threads (single job, runs once)
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt  # 16GB for cd-hit, scales with retries
    params:
        outdir=OUTDIR,
        species=SPECIES_LIST,
        repspec=REPSPEC,
        skip_clustering=config.get("skip_clustering", False),
        cluster_identity=config.get("clustering_identity", 0.8),
        cluster_coverage=config.get("clustering_coverage", 0.8)
    run:
        import os
        shell("mkdir -p {params.outdir}/combinedLibraries")
        
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
            shell(f"cp {combined_file} {{output.combined}}")
            shell(f"rm -f {combined_file}")
            # cd-hit-est is not run, so no .clstr is produced.
            # Touch an empty sentinel so Snakemake dependency tracking still works.
            # The saturation_plot script detects an empty file and uses a fallback.
            shell("touch {output.clstr}")
        else:
            # Run cd-hit-est clustering with config parameters
            shell(f"cd-hit-est -d 0 -aS {{params.cluster_coverage}} -c {{params.cluster_identity}} "
                  f"-G 0 -g 1 -b 500 -r 1 "
                  f"-i {combined_file} -o {{output.combined}} "
                  f"-M {{resources.mem_mb}} -T {{threads}}")
            
            # Clean up intermediate files
            shell(f"rm -f {combined_file}")
        
        # Ensure file is fully written to disk before dependent jobs start
        shell("sync")
        shell("sleep 2")