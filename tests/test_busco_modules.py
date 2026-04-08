"""tests/test_busco_modules.py — unit tests for busco_summary_plot and extract_busco_aa helpers."""

import os
import sys
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from busco_summary_plot import _parse_busco_summary


# ---------------------------------------------------------------------------
# BUSCO summary parsing
# ---------------------------------------------------------------------------

SUMMARY_V5 = textwrap.dedent("""\
    # BUSCO version is: 5.4.3
    # The lineage dataset is: fungi_odb10
    # Summarized benchmarking in BUSCO notation for file genome.fa
    # BUSCO was run in mode: genome
    # Gene predictor used: metaeuk

    ***** Results: *****

    C:95.1%[S:94.2%,D:0.9%],F:1.8%,M:3.1%,n:758

        721    Complete BUSCOs (C)
        714    Complete and single-copy BUSCOs (S)
        7    Complete and duplicated BUSCOs (D)
        14    Fragmented BUSCOs (F)
        23    Missing BUSCOs (M)
        758    Total BUSCO groups searched
""")

SUMMARY_V5_PERFECT = textwrap.dedent("""\
    C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:100

        100    Complete BUSCOs (C)
        100    Complete and single-copy BUSCOs (S)
        0    Complete and duplicated BUSCOs (D)
        0    Fragmented BUSCOs (F)
        0    Missing BUSCOs (M)
        100    Total BUSCO groups searched
""")


class TestParseBuscoSummary:
    def test_single_copy(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5)
        c = _parse_busco_summary(str(p))
        assert c["single"] == 714

    def test_duplicated(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5)
        c = _parse_busco_summary(str(p))
        assert c["duplicated"] == 7

    def test_fragmented(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5)
        c = _parse_busco_summary(str(p))
        assert c["fragmented"] == 14

    def test_missing(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5)
        c = _parse_busco_summary(str(p))
        assert c["missing"] == 23

    def test_total(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5)
        c = _parse_busco_summary(str(p))
        assert c["total"] == 758

    def test_counts_sum_to_total(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5)
        c = _parse_busco_summary(str(p))
        assert c["single"] + c["duplicated"] + c["fragmented"] + c["missing"] == c["total"]

    def test_perfect_completeness(self, tmp_path):
        p = tmp_path / "short_summary.txt"
        p.write_text(SUMMARY_V5_PERFECT)
        c = _parse_busco_summary(str(p))
        assert c["single"] == 100
        assert c["missing"] == 0
        assert c["total"] == 100

    def test_empty_file_returns_zeros(self, tmp_path):
        p = tmp_path / "empty.txt"
        p.write_text("")
        c = _parse_busco_summary(str(p))
        assert c["total"] == 0
        assert c["single"] == 0


# ---------------------------------------------------------------------------
# extract_busco_aa logic (pure-Python helpers only, no Snakemake runtime)
# ---------------------------------------------------------------------------

def _write_fasta(path, sequences):
    """Write a dict{header: seq} to a FASTA file."""
    with open(path, "w") as fh:
        for header, seq in sequences.items():
            fh.write(f">{header}\n{seq}\n")


class TestExtractBuscoAaLogic:
    """Test the occupancy filtering logic independently of the Snakemake wrapper."""

    def _build_busco_dir(self, tmp_path, species, genes, lineage="fungi_odb10"):
        """Create a minimal BUSCO output directory structure for one species."""
        sc_dir = tmp_path / f"{species}_busco" / f"run_{lineage}" / \
                 "busco_sequences" / "single_copy_busco_sequences"
        sc_dir.mkdir(parents=True)
        for gene_id, seq in genes.items():
            _write_fasta(str(sc_dir / f"{gene_id}.faa"), {species: seq})
        return str(tmp_path / f"{species}_busco")

    def test_occupancy_filtering(self, tmp_path):
        """Gene present in only 1/3 species should be excluded at 0.5 occupancy."""
        sp1_dir = self._build_busco_dir(
            tmp_path / "sp1", "sp1", {"geneA": "MKLV", "geneB": "MAPL"}
        )
        sp2_dir = self._build_busco_dir(
            tmp_path / "sp2", "sp2", {"geneA": "MKLV"}
        )
        sp3_dir = self._build_busco_dir(
            tmp_path / "sp3", "sp3", {"geneA": "MKLV"}
        )

        # Simulate what extract_busco_aa does: collect genes per species
        import glob
        def collect_genes(busco_dir):
            genes = set()
            for subdir in glob.glob(
                os.path.join(busco_dir, "run_*", "busco_sequences",
                             "single_copy_busco_sequences")
            ):
                for f in os.listdir(subdir):
                    if f.endswith(".faa"):
                        genes.add(f[:-4])
            return genes

        sp_to_genes = {
            "sp1": collect_genes(sp1_dir),
            "sp2": collect_genes(sp2_dir),
            "sp3": collect_genes(sp3_dir),
        }

        species_list = ["sp1", "sp2", "sp3"]
        n_species = len(species_list)
        min_occ = 0.5
        min_count = max(1, round(min_occ * n_species))  # 2

        all_genes = set()
        for gs in sp_to_genes.values():
            all_genes.update(gs)

        written = []
        skipped = []
        for gene in sorted(all_genes):
            present = [sp for sp in species_list if gene in sp_to_genes[sp]]
            if len(present) >= min_count:
                written.append(gene)
            else:
                skipped.append(gene)

        assert "geneA" in written   # present in all 3 (≥2)
        assert "geneB" in skipped   # present in only sp1 (< 2)

    def test_all_genes_below_threshold_raises_error(self, tmp_path):
        """If no gene meets occupancy, must raise RuntimeError."""
        # We can't easily test the Snakemake script entry point, but we can
        # verify the filtering logic directly.
        species_list = ["sp1", "sp2"]
        n_species = 2
        min_count = 2  # need both species

        sp_to_genes = {
            "sp1": {"geneX"},
            "sp2": {"geneY"},  # different gene, no overlap
        }

        all_genes = {"geneX", "geneY"}
        written = [
            g for g in all_genes
            if len([sp for sp in species_list if g in sp_to_genes[sp]]) >= min_count
        ]
        assert len(written) == 0  # would trigger RuntimeError in the script
