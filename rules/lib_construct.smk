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
    Pre-build the RepeatMasker general library BLAST cache before parallel genome
    jobs start. On a freshly configured conda environment, the cache does not yet
    exist. Running multiple RepeatMasker processes simultaneously causes a race
    where each tries to write the same cache directory and all but one fail.
    Running a single warmup job first ensures the cache is in place so that all
    parallel RepeatMasker jobs can proceed without conflict.
    """
    output:
        sentinel=touch(f"{OUTDIR}/.repeatmasker_cache_ready")
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        """
        RM_SHARE=$(which RepeatMasker | sed 's|/bin/RepeatMasker$|/share/RepeatMasker|')
        if ! find "$RM_SHARE/Libraries" -maxdepth 2 -type d -name "general" 2>/dev/null | grep -q .; then
            echo "Warming up RepeatMasker general library cache..." >&2
            tmp=$(mktemp -d)
            printf '>dummy\\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\\n' > "$tmp/dummy.fa"
            RepeatMasker -lib "$tmp/dummy.fa" -norna -no_is -pa 1 -dir "$tmp" "$tmp/dummy.fa" > /dev/null 2>&1 || true
            rm -rf "$tmp"
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
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    params:
        script_dir=SCRIPT_DIR
    shell:
        """
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
        mkdir -p {params.outdir}
        cd {params.outdir}
        RepeatMasker -species {params.rep_spec} -norna -no_is -lcambig -s -a \
            -pa {params.rm_threads} -dir {params.outdir} $(realpath {input.genome})
        """

rule repeatmasker_custom:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        lib=CUSTOM_LIB if CUSTOM_LIB else [],
        _cache=f"{OUTDIR}/.repeatmasker_cache_ready"
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.prep.masked"
    threads: lambda wildcards: max(1, min(workflow.cores, 64)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=10080
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_RepeatMasker",
        rm_threads=lambda wildcards, threads: max(1, threads // 4)  # RepeatMasker -pa value (uses 4x this)
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        RepeatMasker -lib $(realpath {input.lib}) -norna -no_is -lcambig -s -a \
            -pa {params.rm_threads} -dir {params.outdir} $(realpath {input.genome})
        """

rule extract_repeatmasker_library:
    output:
        replib=f"{{outdir}}/{REPSPEC}.RepeatMasker.lib"
    params:
        repspec=REPSPEC,
        outdir=OUTDIR
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        """
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
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=120
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_Database",
        name="{species}"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        BuildDatabase -name {params.name} {input.masked}
        """

rule repeatmodeler:
    input:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nhr",
        nin="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nin",
        nsq="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nsq"
    output:
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
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
        mkdir -p {params.rm_dir}
        cd {params.rm_dir}
        RepeatModeler -threads {threads} -database {params.db_dir}/{params.db_name}
        """

rule testrainer:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
    output:
        strained="{outdir}/{species}_EarlGrey/{species}_strainer/{species}-families.fa.strained",
        summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained"
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

