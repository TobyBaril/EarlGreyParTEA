import os

GENOME = config["genome"]
SPECIES_LIST = config["species"]
OUTDIR = config["output_dir"] #os.path.join(config["output_dir"], f"{SPECIES}_EarlGrey")
REPSPEC = config["repeatmasker_species"]
CUSTOM_LIB = config["custom_library"]
ITER = config["iterations"]
FLANK = config["flank"]
MAX_SEQ = config["max_consensus_seqs"]
MIN_SEQ = config["min_consensus_seqs"]
HELI = config.get("run_heliano", False)
SCRIPT_DIR = config["script_dir"]

# Rule priority: if CUSTOM_LIB is set, use repeatmasker_custom; otherwise use repeatmasker
if CUSTOM_LIB:
    ruleorder: repeatmasker_custom > repeatmasker
elif REPSPEC:
    ruleorder: repeatmasker > repeatmasker_custom

rule repeatmasker_warmup:
    """
    Pre-build the RepeatMasker general and species-specific library BLAST caches
    before parallel genome jobs start. On a freshly configured conda environment,
    the caches do not yet exist. Running multiple RepeatMasker processes
    simultaneously causes a race where each tries to write the same cache
    directory and all but one fail. Running a single warmup job first ensures
    both caches are fully in place (including refineableHash.dat and
    speciesMeta.pm) so that all parallel RepeatMasker jobs can proceed without
    conflict.
    """
    input:
        f"{OUTDIR}/.repeatmasker_configuration_done"
    output:
        sentinel=touch(f"{OUTDIR}/.repeatmasker_cache_ready")
    log:
        f"{OUTDIR}/.repeatmasker_warmup.log"
    threads: 1
    resources:
        mem_mb=32000,
        runtime=120
    params:
        rep_spec=REPSPEC
    shell:
        """
        exec > {log} 2>&1
        RM_SHARE=$(which RepeatMasker | sed 's|/bin/RepeatMasker$|/share/RepeatMasker|')

        # Warm up the general library cache if not already built
        if ! find "$RM_SHARE/Libraries" -maxdepth 2 -type d -name "general" 2>/dev/null | grep -q .; then
            echo "Warming up RepeatMasker general library cache..." >&2
            tmp=$(mktemp -d)
            printf '>dummy\\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\\n' > "$tmp/dummy.fa"
            RepeatMasker -lib "$tmp/dummy.fa" -no_is -pa 1 -dir "$tmp" "$tmp/dummy.fa" > /dev/null 2>&1 || true
            rm -rf "$tmp"
        fi

        # Warm up the species-specific library cache if a species is specified
        # This prevents a race condition where parallel jobs all try to build
        # the same species cache (including refineableHash.dat / speciesMeta.pm)
        # simultaneously, causing all but one to fail.
        if [ -n "{params.rep_spec}" ]; then
            SPECIES_WORD=$(echo "{params.rep_spec}" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
            # Dynamically discover the CONS cache parent — its name varies depending on
            # whether RepBase is included (e.g. CONS-Dfam_withRBRM_3.9 vs CONS-Dfam_3.9).
            CACHE_PARENT=$(find "$RM_SHARE/Libraries" -maxdepth 1 -type d -name "CONS-*" 2>/dev/null | head -n 1)
            if [ -z "$CACHE_PARENT" ]; then
                echo "WARNING: No CONS-* cache directory found under $RM_SHARE/Libraries — skipping species cache check." >&2
                exit 0
            fi
            CACHE_DIR="$CACHE_PARENT/$SPECIES_WORD"

            # If the cache directory exists but refineableHash.dat is missing, it is
            # incomplete (e.g. a previous OOM-killed makeblastdb run). RepeatMasker
            # considers *.nhr sufficient to skip rebuilding, so we must delete the
            # incomplete directory to force a fresh build.
            if [ -d "$CACHE_DIR" ] && [ ! -f "$CACHE_DIR/refineableHash.dat" ]; then
                echo "Incomplete species cache detected (missing refineableHash.dat). Removing $CACHE_DIR to force rebuild..." >&2
                rm -rf "$CACHE_DIR"
            fi
            # Also remove any stale .working directory left by a previous aborted run.
            rm -rf "$CACHE_DIR.working" 2>/dev/null || true

            if [ ! -f "$CACHE_DIR/refineableHash.dat" ]; then
                echo "Warming up RepeatMasker species library cache for {params.rep_spec}..." >&2
                tmp=$(mktemp -d)
                printf '>dummy\\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\\n' > "$tmp/dummy.fa"
                RepeatMasker -species "{params.rep_spec}" -no_is -pa 1 -dir "$tmp" "$tmp/dummy.fa" > /dev/null 2>&1 || true
                rm -rf "$tmp"
                if [ ! -f "$CACHE_DIR/refineableHash.dat" ]; then
                    echo "ERROR: species cache build failed — refineableHash.dat still missing in $CACHE_DIR" >&2
                    exit 1
                fi
                echo "Species cache built successfully." >&2
            fi
        fi
        """

def get_masked_genome_input(wildcards):
    """Return appropriate input based on config settings"""
    if REPSPEC:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.masked"
    elif CUSTOM_LIB:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.masked"
    else:
        # No masking needed, use prep genome directly
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}.prep"

rule prep_genome:
    input:
        genome=lambda wildcards: GENOME[wildcards.species]
    output:
        gen_prep="{OUTDIR}/{species}_EarlGrey/{species}.prep",
        gen_dict="{OUTDIR}/{species}_EarlGrey/{species}.dict",
        backup="{OUTDIR}/{species}_EarlGrey/{species}.bak.gz"
    log:
        "{OUTDIR}/{species}_EarlGrey/{species}.prep_genome.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    params:
        script_dir=SCRIPT_DIR
    shell:
        """
        exec > {log} 2>&1
        # Create backup of original genome in output directory
        # Remove any existing .orig file first (in case of rerun with read-only permissions)
        rm -f {output.gen_prep}.orig
        cp {input.genome} {output.gen_prep}.orig
        gzip -c {output.gen_prep}.orig > {output.backup}

        # Process genome
        sed '/>/ s/[[:space:]].*//g; /^$/d' {output.gen_prep}.orig > {output.gen_prep}.tmp
        {params.script_dir}/headSwap.sh -i {output.gen_prep}.tmp -o {output.gen_prep}
        rm -f {output.gen_prep}.tmp {output.gen_prep}.orig

        # Move dictionary file
        mv {output.gen_prep}.tmp.dict {output.gen_dict}

        # Replace ambiguous nucleotides
        sed -i.bak '/^>/! s/[DVHBPE]/N/g' {output.gen_prep}
        rm -f {output.gen_prep}.bak
        """

def get_masked_genome_input(wildcards):
    """Return appropriate input based on config settings"""
    if REPSPEC:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.prep.masked"
    elif CUSTOM_LIB:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.prep.masked"
    else:
        # No masking needed, use prep genome directly
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}.prep"

rule repeatmasker:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        _cache=f"{OUTDIR}/.repeatmasker_cache_ready"
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.prep.masked"
    log:
        "{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.repeatmasker.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 64)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=10080
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_RepeatMasker",
        rep_spec=REPSPEC,
        rm_threads=lambda wildcards, threads: max(1, threads // 4)  # RepeatMasker -pa value (uses 4x this)
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        cd {params.outdir}
        RepeatMasker -species {params.rep_spec} -no_is -lcambig -s -a \
            -pa {params.rm_threads} -dir {params.outdir} $(realpath {input.genome})
        """

rule repeatmasker_custom:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        lib=CUSTOM_LIB if CUSTOM_LIB else [],
        _cache=f"{OUTDIR}/.repeatmasker_cache_ready"
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.prep.masked"
    log:
        "{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.repeatmasker.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 64)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=10080
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_RepeatMasker",
        rm_threads=lambda wildcards, threads: max(1, threads // 4)  # RepeatMasker -pa value (uses 4x this)
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        cd {params.outdir}
        RepeatMasker -lib $(realpath {input.lib}) -no_is -lcambig -s -a \
            -pa {params.rm_threads} -dir {params.outdir} $(realpath {input.genome})
        """

rule extract_repeatmasker_library:
    output:
        replib=f"{{outdir}}/{REPSPEC}.RepeatMasker.lib"
    log:
        f"{{outdir}}/{REPSPEC}.extract_repeatmasker_library.log"
    params:
        repspec=REPSPEC,
        outdir=OUTDIR
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        """
        exec > {log} 2>&1
        # Determine RepeatMasker library path
        if [[ $(which RepeatMasker) == *"bin"* ]]; then
            libpath="$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/Libraries/famdb/|')"
            export PATH=$PATH:"$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/|g')"
        else
            libpath="$(which RepeatMasker | sed 's|/[^/]*$||g')/Libraries/famdb/"
        fi

        # Create output directory
        mkdir -p {params.outdir}

        # Extract RepeatMasker library for specified species/clade
        famdb.py -i $libpath families -f fasta_name --include-class-in-name -a -d --curated {params.repspec} > {output.replib}
        """

rule build_db:
    input:
        masked=get_masked_genome_input
    output:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nhr",
        nin="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nin",
        nsq="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nsq"
    log:
        "{outdir}/{species}_EarlGrey/{species}_Database/{species}.build_db.log"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=120
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_Database",
        name="{species}"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        cd {params.outdir}
        BuildDatabase -name {params.name} {input.masked}
        """

rule repeatmodeler:
    input:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nhr",
        nin="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nin",
        nsq="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nsq",
        prep="{outdir}/{species}_EarlGrey/{species}.prep"
    output:
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
    log:
        "{outdir}/{species}_EarlGrey/{species}_RepeatModeler/{species}.repeatmodeler.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 64)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
        runtime=10080
    params:
        db_dir="{outdir}/{species}_EarlGrey/{species}_Database",
        rm_dir="{outdir}/{species}_EarlGrey/{species}_RepeatModeler",
        db_name="{species}"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.rm_dir}
        cd {params.rm_dir}

        # Compute sampable genome size: sum of contigs >= 40 kb only.
        # RepeatModeler discards contigs shorter than 40 kb during sampling, so using
        # the total genome size would overestimate the sequence available per round.
        GENOME_SIZE=$(awk '/^>/{{if(len>=40000)sum+=len; len=0; next}}{{len+=length($0)}} END{{if(len>=40000)sum+=len; print sum+0}}' {input.prep})
        echo "Samplable genome size (contigs >= 40 kb): $GENOME_SIZE bp"

        # Set -genomeSampleSizeMax to the highest RECON round threshold the genome
        # can support. Thresholds are cumulative across rounds (r2=3M, r3=9M, r4=27M,
        # r5=81M, r6=243M), so the sampable size must cover all rounds up to the cap:
        #   >= 363M (3+9+27+81+243): no cap  -> all 6 rounds
        #   >= 120M (3+9+27+81):     cap 81M  -> rounds 1-5
        #   >=  39M (3+9+27):        cap 27M  -> rounds 1-4
        #   >=  12M (3+9):           cap  9M  -> rounds 1-3
        #   <   12M:                 cap  3M  -> rounds 1-2
        SAMPLE_FLAG=""
        if [ "$GENOME_SIZE" -ge 363000000 ] 2>/dev/null; then
            : # genome large enough for all rounds; use RepeatModeler default
        elif [ "$GENOME_SIZE" -ge 120000000 ] 2>/dev/null; then
            echo "Capping at round 5 (-genomeSampleSizeMax 81000000)"
            SAMPLE_FLAG="-genomeSampleSizeMax 81000000"
        elif [ "$GENOME_SIZE" -ge 39000000 ] 2>/dev/null; then
            echo "Capping at round 4 (-genomeSampleSizeMax 27000000)"
            SAMPLE_FLAG="-genomeSampleSizeMax 27000000"
        elif [ "$GENOME_SIZE" -ge 12000000 ] 2>/dev/null; then
            echo "Capping at round 3 (-genomeSampleSizeMax 9000000)"
            SAMPLE_FLAG="-genomeSampleSizeMax 9000000"
        else
            echo "Capping at round 2 (-genomeSampleSizeMax 3000000)"
            SAMPLE_FLAG="-genomeSampleSizeMax 3000000"
        fi

        RepeatModeler -threads {threads} -database {params.db_dir}/{params.db_name} $SAMPLE_FLAG
        """

rule testrainer:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
    output:
        strained="{outdir}/{species}_EarlGrey/{species}_strainer/{species}-families.fa.strained",
        summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained"
    log:
        "{outdir}/{species}_EarlGrey/{species}_strainer/{species}.testrainer.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 64)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
    resources:
        mem_mb=lambda wildcards, attempt: config.get("total_memory_mb", 32000 * attempt),
        runtime=10080
    params:
        outdir=OUTDIR,
        flank=FLANK,
        iter=ITER,
        max_seq=MAX_SEQ,
        min_seq=MIN_SEQ,
        script_dir=SCRIPT_DIR,
        strainer_dir="{outdir}/{species}_EarlGrey/{species}_strainer"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.strainer_dir}
        cd {params.strainer_dir}
        {params.script_dir}/TEstrainer/TEstrainer_for_earlGrey.sh \
           -g {input.genome} -l {input.families} \
           -t {threads} -f {params.flank} \
           -r {params.iter} -n {params.max_seq} \
           -m {params.min_seq}

        # Find and copy the latest TEstrainer output from subdirectory
        latestDir=$(ls -td {params.strainer_dir}/*/ 2>/dev/null | head -n 1)
        if [ -n "$latestDir" ]; then
            latestFile="${{latestDir}}{wildcards.species}-families.fa.strained"
            if [ -f "$latestFile" ]; then
                cp "$latestFile" {output.strained}
            fi
        fi

        # Add species name to fasta headers
        sed -i.bak "s/>/>{wildcards.species}_/g" {output.strained}

        mv {output.strained}.bak {output.summary}
        """
