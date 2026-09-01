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
    Pre-build the RepeatMasker general and species-specific library caches
    before parallel genome jobs start. On a freshly configured conda environment,
    the caches do not yet exist. Running multiple RepeatMasker processes
    simultaneously causes a race where each tries to write the same cache
    directory and all but one fail. Running a single warmup job first ensures
    both caches are fully in place so that all parallel RepeatMasker jobs can
    proceed without conflict.

    RepeatMasker picks the first *writable* directory from a fixed priority
    list (Libraries/ under the install, then ~/.RepeatMaskerCache, then a
    throwaway temp dir). Depending on permissions on the shared conda env,
    it can land in either the install's Libraries/ dir or the user's home
    cache — so this rule checks both locations rather than assuming one.
    """
    output:
        sentinel=touch(f"{OUTDIR}/.repeatmasker_cache_ready")
    log:
        f"{OUTDIR}/.repeatmasker_warmup.log"
    threads: 1
    resources:
        mem_mb=6000,
        runtime=120
    params:
        rep_spec=REPSPEC
    shell:
        """
        exec > {log} 2>&1
        RM_SHARE=$(which RepeatMasker | sed 's|/bin/RepeatMasker$|/share/RepeatMasker|')
        CACHE_ROOTS="$RM_SHARE/Libraries $HOME/.RepeatMaskerCache"

        # Find a CONS-* library parent dir (consensus/BLAST-format cache root)
        # across every candidate cache location.
        find_cons_parents() {{
            for root in $CACHE_ROOTS; do
                find "$root" -maxdepth 1 -type d -name "CONS-*" 2>/dev/null
            done
        }}

        # Does a completed species cache exist anywhere? Print its path if so.
        find_completed_species_dir() {{
            local species="$1"
            for parent in $(find_cons_parents); do
                dir="$parent/$species"
                if [ -f "$dir/speciesMeta.pm" ] || [ -f "$dir/refineableHash.dat" ]; then
                    echo "$dir"
                    return 0
                fi
            done
            return 1
        }}

        # Any existing (possibly incomplete) species cache dirs, across all roots.
        find_all_species_dirs() {{
            local species="$1"
            for parent in $(find_cons_parents); do
                dir="$parent/$species"
                [ -d "$dir" ] && echo "$dir"
            done
        }}

        # ---- Warm up the general library cache if not already built, in either root ----
        GENERAL_FOUND=0
        for root in $CACHE_ROOTS; do
            if find "$root" -maxdepth 2 -type d -name "general" 2>/dev/null | grep -q .; then
                GENERAL_FOUND=1
                break
            fi
        done
        if [ "$GENERAL_FOUND" -eq 0 ]; then
            echo "Warming up RepeatMasker general library cache..." >&2
            tmp=$(mktemp -d)
            printf '>dummy\\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\\n' > "$tmp/dummy.fa"
            RepeatMasker -lib "$tmp/dummy.fa" -no_is -pa 1 -dir "$tmp" "$tmp/dummy.fa" || true
            rm -rf "$tmp"
        fi

        # ---- Warm up the species-specific library cache if a species is specified ----
        if [ -n "{params.rep_spec}" ]; then
            SPECIES_WORD=$(echo "{params.rep_spec}" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')

            if [ -z "$(find_cons_parents)" ]; then
                echo "WARNING: No CONS-* cache directory found under any of: $CACHE_ROOTS — skipping species cache check." >&2
                exit 0
            fi

            # Clean up any incomplete species cache dirs (e.g. from an OOM-killed
            # makeblastdb run) in EITHER location, so a stale half-built dir in
            # one root doesn't shadow a fresh build attempt.
            for dir in $(find_all_species_dirs "$SPECIES_WORD"); do
                if [ ! -f "$dir/speciesMeta.pm" ] && [ ! -f "$dir/refineableHash.dat" ]; then
                    echo "Incomplete species cache detected. Removing $dir to force rebuild..." >&2
                    rm -rf "$dir"
                fi
                rm -rf "$dir.working" 2>/dev/null || true
            done

            EXISTING=$(find_completed_species_dir "$SPECIES_WORD" || true)
            if [ -n "$EXISTING" ]; then
                echo "Species cache already present at $EXISTING — skipping build." >&2
            else
                echo "Warming up RepeatMasker species library cache for {params.rep_spec}..." >&2
                tmp=$(mktemp -d)
                printf '>dummy\\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\\n' > "$tmp/dummy.fa"
                RepeatMasker -species "{params.rep_spec}" -no_is -pa 1 -dir "$tmp" "$tmp/dummy.fa" || true
                rm -rf "$tmp"

                BUILT=$(find_completed_species_dir "$SPECIES_WORD" || true)
                if [ -z "$BUILT" ]; then
                    echo "ERROR: species cache build failed — no completed cache found for '$SPECIES_WORD' in any of: $CACHE_ROOTS" >&2
                    echo "Checked for CONS-*/$SPECIES_WORD dirs containing speciesMeta.pm or refineableHash.dat." >&2
                    exit 1
                fi
                echo "Species cache built successfully at $BUILT" >&2
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
    threads: lambda wildcards: max(1, min(workflow.cores, 128)) if config.get("slurm_mode", False) or config.get("lsf_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 128))
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
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
    threads: lambda wildcards: max(1, min(workflow.cores, 128)) if config.get("slurm_mode", False) or config.get("lsf_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 128))
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
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
        # Locate famdb library path — supports both:
        #   FamDB 3.0.0+ standalone (share/famdb-*/Libraries/famdb/): Dfam 4.0+
        #   FamDB < 3.0.0 embedded  (share/RepeatMasker/Libraries/famdb/):   Dfam 3.9
        CONDA_PREFIX_RM=$(which RepeatMasker | sed 's|/bin/RepeatMasker$||')
        FAMDB_SHARE=$(find "$CONDA_PREFIX_RM/share" -maxdepth 1 -type d -name "famdb-*" 2>/dev/null | sort -V | tail -n 1)
        if [ -n "$FAMDB_SHARE" ]; then
            libpath="$FAMDB_SHARE/Libraries/famdb/"
        else
            libpath="$CONDA_PREFIX_RM/share/RepeatMasker/Libraries/famdb/"
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
    threads: lambda wildcards: max(1, min(workflow.cores, 128)) if config.get("slurm_mode", False) or config.get("lsf_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 128))
    resources:
        mem_mb=lambda wildcards, attempt: 64000 * attempt,
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
    threads: lambda wildcards: max(1, min(workflow.cores, 128)) if config.get("slurm_mode", False) or config.get("lsf_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 128))
    resources:
        mem_mb=lambda wildcards, attempt: config.get("total_memory_mb", 64000 * attempt),
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

