# busco_phylo.smk — BUSCO-based phylogenomics pipeline for EarlGreyParTEA.
#
# Workflow
# --------
# 1. run_busco            — run BUSCO genome mode on each prepped genome.
# 2. busco_summary_table  — parse all short_summary.txt files → completeness PDF + TSV.
# 3. extract_busco_aa     — CHECKPOINT: collect single-copy AA FASTAs per gene,
#                           filter by occupancy, write one faa per gene to
#                           {OUTDIR}/buscoPhylo/busco_genes/.
# 4. align_busco_gene     — run mafft on each per-gene faa → trimmed MSA via clipkit.
# 5. create_supermatrix   — concatenate trimmed alignments into a supermatrix.
# 6. run_fasttree         — infer ML tree from the supermatrix with FastTree.
# 7. busco_te_qc          — scatter plot of BUSCO completeness vs TE content.
#
# Configuration keys (all optional with defaults in on_start_functions.py)
# -------------------------------------------------------------------------
#   busco_lineage         : str  — BUSCO lineage dataset (e.g. "fungi_odb10")
#   busco_prefix          : str  — prefix for BUSCO run directory names
#   busco_min_occupancy   : float — fraction of species a gene must be in (default 0.5)
#
# Outputs (all under {OUTDIR}/buscoPhylo/)
# -----------------------------------------
#   busco_completeness.pdf / .tsv
#   busco_genes/           (per-gene FASTAs; checkpoint output)
#   aligned/<gene>.clipkit.fa
#   supermatrix.fa
#   species_tree.nwk
#   busco_te_qc.pdf / .tsv

# ---------------------------------------------------------------------------
# 0. Download BUSCO lineage dataset once (before parallel runs)
# ---------------------------------------------------------------------------

rule fetch_busco_db:
    """Download the BUSCO lineage dataset to a shared path before any run_busco
    jobs start.  All run_busco instances depend on this sentinel so the
    download only happens once regardless of how many species are processed
    in parallel."""
    output:
        sentinel=f"{OUTDIR}/buscoPhylo/.busco_db_ready",
    log:
        f"{OUTDIR}/buscoPhylo/fetch_busco_db.log"
    params:
        lineage=BUSCO_LINEAGE,
        db_path=f"{OUTDIR}/buscoPhylo/busco_db",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=120,
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.db_path}
        busco --download {params.lineage} --download_path {params.db_path}
        touch {output.sentinel}
        """


# ---------------------------------------------------------------------------
# 1. Run BUSCO for each species
# ---------------------------------------------------------------------------

rule run_busco:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        db_ready=f"{OUTDIR}/buscoPhylo/.busco_db_ready",
    output:
        summary=f"{{outdir}}/{{species}}_EarlGrey/{{species}}_busco/short_summary.specific.{BUSCO_LINEAGE}.{{species}}_busco.txt",
        busco_dir=directory("{outdir}/{species}_EarlGrey/{species}_busco"),
    log:
        "{outdir}/{species}_EarlGrey/{species}_busco/{species}.run_busco.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 16)) if config.get("slurm_mode", False) or config.get("lsf_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 16))
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=480,
    params:
        lineage=BUSCO_LINEAGE,
        busco_prefix=BUSCO_PREFIX,
        outdir="{outdir}/{species}_EarlGrey",
        db_path=f"{OUTDIR}/buscoPhylo/busco_db",
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        cd {params.outdir}
        busco \
            -i $(realpath {input.genome}) \
            -l {params.lineage} \
            -o {wildcards.species}_busco \
            -m genome \
            -c {threads} \
            --force \
            --download_path {params.db_path}
        """


# ---------------------------------------------------------------------------
# 2. BUSCO completeness summary table and plot
# ---------------------------------------------------------------------------

rule busco_summary_table:
    input:
        summaries=expand(
            "{outdir}/{species}_EarlGrey/{species}_busco/"
            "short_summary.specific.{lineage}.{species}_busco.txt",
            outdir=OUTDIR,
            species=SPECIES_LIST,
            lineage=BUSCO_LINEAGE,
        ),
    output:
        pdf=f"{OUTDIR}/buscoPhylo/busco_completeness.pdf",
        tsv=f"{OUTDIR}/buscoPhylo/busco_completeness.tsv",
    log:
        f"{OUTDIR}/buscoPhylo/busco_summary_table.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        species=SPECIES_LIST,
    script:
        "../scripts/busco_summary_plot.py"


# ---------------------------------------------------------------------------
# 3. Checkpoint: extract per-gene single-copy amino acid FASTAs
# ---------------------------------------------------------------------------

checkpoint extract_busco_aa:
    input:
        busco_dirs=expand(
            "{outdir}/{species}_EarlGrey/{species}_busco",
            outdir=OUTDIR,
            species=SPECIES_LIST,
        ),
    output:
        gene_dir=directory(f"{OUTDIR}/buscoPhylo/busco_genes"),
        occupancy_tsv=f"{OUTDIR}/buscoPhylo/busco_gene_occupancy.tsv",
    log:
        f"{OUTDIR}/buscoPhylo/extract_busco_aa.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        species=SPECIES_LIST,
        busco_prefix=BUSCO_PREFIX,
        busco_min_occupancy=BUSCO_MIN_OCC,
        outdir=OUTDIR,
    script:
        "../scripts/extract_busco_aa.py"


# ---------------------------------------------------------------------------
# 4. Align each gene and trim with clipkit
# ---------------------------------------------------------------------------

rule align_busco_gene:
    input:
        faa=f"{OUTDIR}/buscoPhylo/busco_genes/{{gene_id}}.faa",
    output:
        aln=f"{OUTDIR}/buscoPhylo/aligned/{{gene_id}}.clipkit.fa",
    log:
        f"{OUTDIR}/buscoPhylo/aligned/{{gene_id}}.align.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
    shell:
        """
        exec > {log} 2>&1
        mkdir -p $(dirname {output.aln})
        mafft --auto --quiet {input.faa} > {output.aln}.mafft.fa
        clipkit {output.aln}.mafft.fa -m smart-gap -o {output.aln}
        rm -f {output.aln}.mafft.fa
        # Failsafe: clipkit exits 0 even when all columns are trimmed away.
        # If the output has no alignment columns (no non-header lines), remove it
        # so it is excluded from supermatrix construction.
        if [ ! -s {output.aln} ] || ! grep -q -v '^>' {output.aln} 2>/dev/null; then
            echo "[align_busco_gene] WARNING: empty alignment after trimming, removing {output.aln}" >&2
            rm -f {output.aln}
            touch {output.aln}  # keep a zero-byte sentinel so Snakemake is satisfied
        fi
        """


def _get_aligned_genes(wildcards):
    """Resolve checkpoint output to a list of aligned gene paths."""
    checkpoint_output = checkpoints.extract_busco_aa.get(**wildcards).output.gene_dir
    gene_ids, = glob_wildcards(os.path.join(checkpoint_output, "{gene_id}.faa"))
    return expand(
        f"{OUTDIR}/buscoPhylo/aligned/{{gene_id}}.clipkit.fa",
        gene_id=gene_ids,
    )


# ---------------------------------------------------------------------------
# 5. Concatenate trimmed alignments → supermatrix
# ---------------------------------------------------------------------------

rule create_supermatrix:
    input:
        alignments=_get_aligned_genes,
    output:
        supermatrix=f"{OUTDIR}/buscoPhylo/supermatrix.fa",
        partition=f"{OUTDIR}/buscoPhylo/supermatrix.partitions",
    log:
        f"{OUTDIR}/buscoPhylo/create_supermatrix.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        aligned_dir=f"{OUTDIR}/buscoPhylo/aligned",
        species=SPECIES_LIST,
    run:
        import os
        import glob

        aln_files = sorted(glob.glob(os.path.join(params.aligned_dir, "*.clipkit.fa")))
        if not aln_files:
            raise RuntimeError("No aligned gene files found for supermatrix construction.")

        # Read each alignment into a dict {species: seq}
        def _read_fasta(path):
            seqs = {}
            current = None
            with open(path) as fh:
                for line in fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        current = line[1:].strip()
                        seqs[current] = []
                    elif current:
                        seqs[current].append(line)
            return {k: "".join(v) for k, v in seqs.items()}

        species_seqs = {sp: [] for sp in params.species}
        partitions = []
        pos = 1

        for aln_f in aln_files:
            # Skip zero-byte sentinel files left by align_busco_gene failsafe
            if os.path.getsize(aln_f) == 0:
                print(f"[create_supermatrix] Skipping empty alignment: {os.path.basename(aln_f)}",
                      flush=True)
                continue
            gene_aln = _read_fasta(aln_f)
            gene_id = os.path.basename(aln_f).replace(".clipkit.fa", "")
            # Determine alignment length from first entry
            if not gene_aln:
                continue
            aln_len = len(next(iter(gene_aln.values())))
            # Gaps for missing species
            for sp in params.species:
                seq = gene_aln.get(sp, "-" * aln_len)
                # Pad or truncate to aln_len if necessary
                if len(seq) < aln_len:
                    seq = seq + "-" * (aln_len - len(seq))
                elif len(seq) > aln_len:
                    seq = seq[:aln_len]
                species_seqs[sp].append(seq)
            partitions.append(f"AA, {gene_id} = {pos}-{pos + aln_len - 1}")
            pos += aln_len

        with open(output.supermatrix, "w") as ofh:
            for sp in params.species:
                cat_seq = "".join(species_seqs[sp])
                ofh.write(f">{sp}\n{cat_seq}\n")

        with open(output.partition, "w") as pfh:
            pfh.write("\n".join(partitions) + "\n")

        print(
            f"[create_supermatrix] Wrote supermatrix with {len(aln_files)} genes, "
            f"{pos - 1} positions.",
            flush=True,
        )


# ---------------------------------------------------------------------------
# 6. Infer phylogenetic tree with FastTree
# ---------------------------------------------------------------------------

rule run_fasttree:
    input:
        supermatrix=f"{OUTDIR}/buscoPhylo/supermatrix.fa",
    output:
        tree=f"{OUTDIR}/buscoPhylo/species_tree.nwk",
    log:
        f"{OUTDIR}/buscoPhylo/fasttree.log"
    threads: lambda wildcards: max(1, min(workflow.cores, 8))
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt,
        runtime=240,
    shell:
        """
        export OMP_NUM_THREADS={threads}
        FastTree -lg -gamma {input.supermatrix} > {output.tree} 2>> {log}
        """


# ---------------------------------------------------------------------------
# 7. BUSCO completeness plot reordered by phylogeny
# ---------------------------------------------------------------------------

rule busco_completeness_phylo:
    input:
        tsv=f"{OUTDIR}/buscoPhylo/busco_completeness.tsv",
        tree=f"{OUTDIR}/buscoPhylo/species_tree.nwk",
    output:
        pdf=f"{OUTDIR}/buscoPhylo/busco_completeness_phylo.pdf",
    log:
        f"{OUTDIR}/buscoPhylo/busco_completeness_phylo.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        species=SPECIES_LIST,
    script:
        "../scripts/busco_summary_plot.py"


# ---------------------------------------------------------------------------
# 8. BUSCO completeness vs TE content QC scatter
# ---------------------------------------------------------------------------

rule busco_te_qc:
    input:
        busco_summaries=expand(
            "{outdir}/{species}_EarlGrey/{species}_busco/"
            "short_summary.specific.{lineage}.{species}_busco.txt",
            outdir=OUTDIR,
            species=SPECIES_LIST,
            lineage=BUSCO_LINEAGE,
        ),
        cov_tsv=f"{OUTDIR}/sharedUniqueContent/shared_unique_coverage.tsv",
    output:
        pdf=f"{OUTDIR}/buscoPhylo/busco_te_qc.pdf",
        tsv=f"{OUTDIR}/buscoPhylo/busco_te_qc.tsv",
    log:
        f"{OUTDIR}/buscoPhylo/busco_te_qc.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        species=SPECIES_LIST,
    script:
        "../scripts/busco_te_qc_plot.py"
