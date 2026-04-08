import os

# Configuration variables (SOFTMASK and MARGIN are defined in main Snakefile)
GENOME = config["genome"]
SPECIES_LIST = config["species"]
OUTDIR = config["output_dir"]
REPSPEC = config.get("repeatmasker_species", None)
# Handle both heliano and run_heliano config keys, convert boolean to yes/no
_heliano_val = config.get("heliano", config.get("run_heliano", False))
HELIANO = "yes" if (_heliano_val is True or _heliano_val == "yes") else "no"
SCRIPT_DIR = config["script_dir"]  # Path to EarlGrey scripts directory

rule repeatmasker_annotation:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        library=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa",
        _cache=f"{OUTDIR}/.repeatmasker_cache_ready"
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.masked",
        out="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.out",
        tbl="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.tbl"
    threads: lambda wildcards: max(1, min(workflow.cores, 64)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=10080
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library",
        rm_threads=lambda wildcards, threads: max(1, threads // 4)  # RepeatMasker -pa value (uses 4x this)
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        RepeatMasker -lib $(realpath {input.library}) -norna -no_is -lcambig -s -a \
            -pa {params.rm_threads} -dir {params.outdir} $(realpath {input.genome})
        """

rule heliano_detection:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep"
    output:
        helitron_gff="{outdir}/{species}_EarlGrey/{species}_heliano/RC.representative.gff"
    threads: lambda wildcards: max(1, min(workflow.cores, 32)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 32))
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=480
    params:
        heliano_dir="{outdir}/{species}_EarlGrey/{species}_heliano"
    shell:
        """
        if [ "{HELIANO}" == "yes" ]; then
            mkdir -p {params.heliano_dir}
            cd {params.heliano_dir}
            timestamp=$(date +"%Y%m%d_%H%M")
            heliano -g {input.genome} --nearest -dn 6000 -flank_sim 0.5 \
                    -o {params.heliano_dir}/HEL_$timestamp -w 10000 -n {threads}
            awk '{{OFS="\t"}}{{print $1, "HELIANO", "RC/Helitron", $2+1, $3, $5, $6, ".", "ID="$9"_"$11";shortTE=F"}}' \
                {params.heliano_dir}/HEL_$timestamp/RC.representative.bed > {output.helitron_gff}
        else
            mkdir -p {params.heliano_dir}
            touch {output.helitron_gff}
        fi
        """

rule merge_repeats:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        dict="{outdir}/{species}_EarlGrey/{species}.dict",
        out="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.out",
        tbl="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.tbl",
        helitron_gff="{outdir}/{species}_EarlGrey/{species}_heliano/RC.representative.gff" if HELIANO == "yes" else []
    output:
        bed="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.bed",
        gff="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.gff",
        summary="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.summary"
    threads: lambda wildcards: max(1, min(workflow.cores, 16)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 16))
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=240
    params:
        script_dir=SCRIPT_DIR,
        outdir="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge",
        margin=MARGIN,
        helitron_param=lambda wildcards, input: f"-e {input.helitron_gff}" if HELIANO == "yes" and os.path.getsize(input.helitron_gff if HELIANO == "yes" else "/dev/null") > 0 else ""
    shell:
        """
        mkdir -p {params.outdir}
        
        # Try loose merge first
        {params.script_dir}/rcMergeRepeatsLoose -f {input.genome} -s {wildcards.species} \
            -d {params.outdir} -u {input.out} -q {input.tbl} -t {threads} \
            -b {input.dict} -m {params.margin} {params.helitron_param}
        
        # Fix GFF formatting
        if [ -f "{output.gff}" ]; then
            awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $8, toupper($9)}}' {output.gff} > {output.gff}.tmp
            mv {output.gff}.tmp {output.gff}
        fi
        
        # If loose merge failed, try strict merge
        if [ ! -f "{output.bed}" ]; then
            echo "Loose merge failed, trying strict merge..."
            {params.script_dir}/rcMergeRepeats -f {input.genome} -s {wildcards.species} \
                -d {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_mergedRepeats \
                -u {input.out} -q {input.tbl} -t {threads} \
                -b {input.dict} -m {params.margin} {params.helitron_param}
            
            # Move strict merge results to expected location if loose merge failed
            if [ -f "{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_mergedRepeats/{wildcards.species}.filteredRepeats.bed" ]; then
                mkdir -p {params.outdir}
                mv {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_mergedRepeats/{wildcards.species}.filteredRepeats.* {params.outdir}/
            fi
        fi
        """

rule generate_summary_charts:
    input:
        summary="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.summary",
        tbl="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.tbl"
    output:
        pie="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
        highLevelCount="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    params:
        script_dir=SCRIPT_DIR,
        outdir="{outdir}/{species}_EarlGrey/{species}_summaryFiles"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        {params.script_dir}/autoPie.sh -i {input.summary} -t {input.tbl} \
                                       -p {output.pie} -o {output.highLevelCount}
        """

rule calculate_divergence:
    input:
        library=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa",
        genome_orig=lambda wildcards: GENOME[wildcards.species],
        gff="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.gff"
    output:
        div_gff="{outdir}/{species}_EarlGrey/{species}_RepeatLandscape/{species}.filteredRepeats.withDivergence.gff",
        div_summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}_divergence_summary_table.tsv"
    threads: lambda wildcards: max(1, min(workflow.cores, 16)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 16))
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=480
    params:
        script_dir=SCRIPT_DIR,
        landscape_dir="{outdir}/{species}_EarlGrey/{species}_RepeatLandscape",
        summary_dir="{outdir}/{species}_EarlGrey/{species}_summaryFiles",
        divcalc_tmp="/tmp/egdiv_{species}",
    shell:
        """
        mkdir -p {params.landscape_dir}
        cd {params.landscape_dir}
        
        # divergence_calc.py calls pybedtools.set_tempdir() which sets
        # tempfile.tempdir globally. tempfile.tempdir takes priority over TMPDIR,
        # so a long output path would exceed the 108-char AF_UNIX socket limit.
        # Pass -tmp with a short per-species path in /tmp to avoid this.
        mkdir -p {params.divcalc_tmp}
        
        # Calculate divergence
        python {params.script_dir}/divergenceCalc/divergence_calc.py \
            -l {input.library} -g {input.genome_orig} -i {input.gff} \
            -o {output.div_gff} -t {threads} \
            -tmp {params.divcalc_tmp}
        
        # Generate divergence plots
        Rscript {params.script_dir}/divergenceCalc/divergence_plot.R \
            -s {wildcards.species} -g {output.div_gff} -o {params.landscape_dir}
        
        # Copy results to summary directory
        mkdir -p {params.summary_dir}
        cp {params.landscape_dir}/*.pdf {params.summary_dir}/ || true
        cp {params.landscape_dir}/*_summary_table.tsv {output.div_summary} || true
        
        # Update main GFF with divergence info (copy instead of move to keep output file)
        cp {output.div_gff} {input.gff}
        
        # Cleanup
        rm -rf {params.divcalc_tmp} || true
        """

rule sweep_up_files:
    input:
        bed="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.bed",
        gff="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.gff",
        highLevelCount="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt",
        pie="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
        div_summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}_divergence_summary_table.tsv"
    output:
        summary_bed="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed",
        summary_gff="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15
    params:
        summary_dir="{outdir}/{species}_EarlGrey/{species}_summaryFiles"
    shell:
        """
        # Copy final results to summary directory
        cp {input.bed} {output.summary_bed}
        cp {input.gff} {output.summary_gff}
        """

rule generate_softmasked_genome:
    input:
        bed="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed",
        backup="{outdir}/{species}_EarlGrey/{species}.bak.gz"
    output:
        softmasked="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.softmasked.fasta"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        """
        if [ "{SOFTMASK}" == "yes" ]; then
            gunzip -c {input.backup} > {input.backup}.tmp
            bedtools maskfasta -fi {input.backup}.tmp -bed {input.bed} \
                              -fo {output.softmasked} -soft
            rm -f {input.backup}.tmp
        fi
        """