"""extract_busco_aa.py — checkpoint script: gather per-gene BUSCO amino acid FASTAs.

This script is invoked as a Snakemake *checkpoint* (script: directive).  It
scans each species' BUSCO run output directory for single-copy complete BUSCO
genes, collects the amino-acid FASTA from each species that has that gene, and
writes one multi-FASTA per gene to {outdir}/buscoPhylo/busco_genes/.

Only genes present in ≥ busco_min_occupancy fraction of species are written.
Genes that are below this threshold are skipped with a log message.

After this checkpoint completes, downstream rules use glob_wildcards to
discover the set of genes that were actually written.

snakemake.params.species           : list[str]
snakemake.params.busco_prefix      : str  — prefix used in BUSCO run dirs
                                     (e.g. 'busco' → dirs named busco_{species})
snakemake.params.busco_min_occupancy : float (0–1)
snakemake.params.outdir            : str
snakemake.input.busco_dirs         : list of paths to BUSCO run directories
snakemake.output.gene_dir          : directory that will contain per-gene FASTAs
snakemake.output.occupancy_tsv     : TSV summarising gene occupancy
"""

import os
import glob
import shutil

AMINO_ACID_SUBDIR = "run_{lineage}/busco_sequences/single_copy_busco_sequences"


def _find_aa_fasta(busco_dir, gene_id):
    """Return path to the amino-acid FASTA for *gene_id* in *busco_dir*, or None."""
    # BUSCO v5: single_copy_busco_sequences/{gene_id}.faa
    for subdir in glob.glob(os.path.join(busco_dir, "run_*", "busco_sequences",
                                          "single_copy_busco_sequences")):
        faa = os.path.join(subdir, gene_id + ".faa")
        if os.path.isfile(faa):
            return faa
    return None


def main():
    species_list    = list(snakemake.params.species)                       # noqa: F821
    busco_prefix    = snakemake.params.busco_prefix                        # noqa: F821
    min_occ         = float(snakemake.params.busco_min_occupancy)          # noqa: F821
    outdir          = snakemake.params.outdir                              # noqa: F821
    busco_dirs      = list(snakemake.input.busco_dirs)                     # noqa: F821
    gene_dir        = snakemake.output.gene_dir                            # noqa: F821
    occupancy_tsv   = snakemake.output.occupancy_tsv                       # noqa: F821

    os.makedirs(gene_dir, exist_ok=True)
    n_species = len(species_list)

    # Map species → busco_dir
    sp_to_dir = {sp: bd for sp, bd in zip(species_list, busco_dirs)}

    # Collect all single-copy complete gene IDs per species
    sp_to_genes = {}
    for sp, bd in sp_to_dir.items():
        genes = set()
        for subdir in glob.glob(os.path.join(bd, "run_*", "busco_sequences",
                                              "single_copy_busco_sequences")):
            for f in os.listdir(subdir):
                if f.endswith(".faa"):
                    genes.add(f[:-4])
        sp_to_genes[sp] = genes
        print(f"[extract_busco_aa] {sp}: {len(genes)} single-copy complete BUSCOs",
              flush=True)

    # Union of all gene IDs
    all_genes = set()
    for genes in sp_to_genes.values():
        all_genes.update(genes)

    min_count = max(1, round(min_occ * n_species))
    genes_written = []
    genes_skipped = []

    with open(occupancy_tsv, "w") as occ_fh:
        occ_fh.write("gene_id\tn_species\toccupancy_fraction\twritten\n")
        for gene_id in sorted(all_genes):
            present = [sp for sp in species_list if gene_id in sp_to_genes[sp]]
            n_present = len(present)
            occ_frac = n_present / n_species
            written = n_present >= min_count

            occ_fh.write(
                f"{gene_id}\t{n_present}\t{occ_frac:.4f}\t{'yes' if written else 'no'}\n"
            )

            if not written:
                genes_skipped.append(gene_id)
                continue

            # Write multi-FASTA: one record per species
            out_fasta = os.path.join(gene_dir, gene_id + ".faa")
            with open(out_fasta, "w") as ofh:
                for sp in present:
                    faa_path = _find_aa_fasta(sp_to_dir[sp], gene_id)
                    if faa_path is None:
                        continue
                    with open(faa_path) as ifh:
                        for line in ifh:
                            if line.startswith(">"):
                                ofh.write(f">{sp}\n")
                            else:
                                ofh.write(line)
            genes_written.append(gene_id)

    print(
        f"[extract_busco_aa] Genes written (occupancy ≥{min_occ:.0%}): {len(genes_written)}",
        flush=True,
    )
    print(
        f"[extract_busco_aa] Genes skipped (below threshold): {len(genes_skipped)}",
        flush=True,
    )
    if not genes_written:
        raise RuntimeError(
            "[extract_busco_aa] No BUSCO genes met the occupancy threshold. "
            "Lower busco_min_occupancy or check that BUSCO ran successfully."
        )


if "snakemake" in dir():
    main()
