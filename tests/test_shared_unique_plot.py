"""tests/test_shared_unique_plot.py — unit tests for shared_unique_plot.py helpers."""

import os
import sys
import tempfile
import textwrap

import numpy as np
import pytest

# Allow import without a Snakemake runtime
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from shared_unique_plot import (
    CLASS_ORDER,
    CLASS_PALETTE,
    LIGHT_PALETTE,
    OTHER_LABEL,
    _classify_te,
    _lighten,
    _parse_gff_name,
    _gff_coverage_and_families,
    _genome_size_from_prep,
    _empty_class_dict,
    _presence_absence_mode,
)


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------

class TestLighten:
    def test_black_becomes_grey(self):
        result = _lighten("#000000", factor=0.5)
        r = int(result[1:3], 16)
        assert r > 0

    def test_white_stays_white(self):
        assert _lighten("#FFFFFF", factor=0.5) == "#FFFFFF"

    def test_full_factor_gives_white(self):
        assert _lighten("#E32017", factor=1.0) == "#FFFFFF"

    def test_zero_factor_preserves(self):
        assert _lighten("#00782A", factor=0.0) == "#00782A"

    def test_lightened_palette_has_all_classes(self):
        assert set(LIGHT_PALETTE.keys()) == set(CLASS_PALETTE.keys())

    def test_lightened_is_actually_lighter(self):
        for cls, col in CLASS_PALETTE.items():
            light = LIGHT_PALETTE[cls]
            orig_sum = sum(int(col[i:i+2], 16) for i in (1, 3, 5))
            light_sum = sum(int(light[i:i+2], 16) for i in (1, 3, 5))
            assert light_sum >= orig_sum, f"Light palette not lighter for {cls}"


# ---------------------------------------------------------------------------
# TE classification
# ---------------------------------------------------------------------------

class TestClassifyTe:
    def test_ltr_variants(self):
        assert _classify_te("LTR/Gypsy") == "LTR"
        assert _classify_te("LTR/Copia") == "LTR"
        assert _classify_te("LTR") == "LTR"

    def test_line_variants(self):
        assert _classify_te("LINE/R2-Hero") == "LINE"
        assert _classify_te("LINE/RTE-BovB") == "LINE"

    def test_penelope_before_line(self):
        assert _classify_te("LINE/Penelope") == "Penelope"
        assert _classify_te("Penelope") == "Penelope"

    def test_sine(self):
        assert _classify_te("SINE/tRNA") == "SINE"

    def test_dna(self):
        assert _classify_te("DNA/TIR") == "DNA"
        assert _classify_te("DNA/CMC-EnSpm") == "DNA"

    def test_rolling_circle(self):
        assert _classify_te("RC/Helitron") == "Rolling Circle"
        assert _classify_te("RC") == "Rolling Circle"

    def test_other_categories(self):
        assert _classify_te("Simple_repeat") == "Other"
        assert _classify_te("Low_complexity") == "Other"
        assert _classify_te("Satellite") == "Other"
        assert _classify_te("rRNA") == "Other"

    def test_unclassified(self):
        assert _classify_te("Unknown") == "Unclassified"
        assert _classify_te("") == "Unclassified"
        assert _classify_te("ARTEFACT") == "Other"

    def test_whitespace_stripped(self):
        assert _classify_te("  LTR/Gypsy  ") == "LTR"


# ---------------------------------------------------------------------------
# GFF name parsing
# ---------------------------------------------------------------------------

class TestParseGffName:
    def test_typical(self):
        attr = "ID=GENOME1_RND-1_FAMILY-34_3;NAME=GENOME1_RND-1_FAMILY-34;TSTART=1"
        assert _parse_gff_name(attr) == "GENOME1_RND-1_FAMILY-34"

    def test_case_insensitive_key(self):
        assert _parse_gff_name("name=fam1") == "fam1"

    def test_missing_returns_none(self):
        assert _parse_gff_name("ID=foo;SCORE=100") is None

    def test_empty_string(self):
        assert _parse_gff_name("") is None

    def test_space_stripped(self):
        assert _parse_gff_name("NAME= fam1 ") == "fam1"


# ---------------------------------------------------------------------------
# GFF file parsing
# ---------------------------------------------------------------------------

GFF_CONTENT = textwrap.dedent("""\
    # comment line
    chr1\tEarl_Grey\tLTR/Gypsy\t101\t200\t.\t+\t.\tID=sp1_f1_1;NAME=sp1_f1
    chr1\tEarl_Grey\tLINE/R2-Hero\t301\t500\t.\t-\t.\tID=sp1_f2_1;NAME=sp1_f2
    chr1\tEarl_Grey\tLTR/Gypsy\t601\t700\t.\t+\t.\tID=sp1_f1_2;NAME=sp1_f1
    chr1\tEarl_Grey\tSimple_repeat\t801\t850\t.\t+\t.\tID=sp1_sr_1;NAME=sp1_sr
""")


class TestGffParsing:
    def test_family_class_mapping(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text(GFF_CONTENT)
        fc, hits = _gff_coverage_and_families(str(gff))
        assert fc["sp1_f1"] == "LTR"
        assert fc["sp1_f2"] == "LINE"
        assert fc["sp1_sr"] == "Other"

    def test_hit_count(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text(GFF_CONTENT)
        _, hits = _gff_coverage_and_families(str(gff))
        assert len(hits) == 4

    def test_bp_calculation(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text(GFF_CONTENT)
        _, hits = _gff_coverage_and_families(str(gff))
        # First hit: 101-200 = 100 bp
        name, cls, bp = hits[0]
        assert name == "sp1_f1"
        assert bp == 100

    def test_comment_and_empty_skipped(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text("# header\n\n" + GFF_CONTENT)
        _, hits = _gff_coverage_and_families(str(gff))
        assert len(hits) == 4


# ---------------------------------------------------------------------------
# Genome size from FASTA
# ---------------------------------------------------------------------------

class TestGenomeSizeFromPrep:
    def test_basic(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_text(">chr1\nACGTACGT\n>chr2\nACGT\n")
        assert _genome_size_from_prep(str(fa)) == 12

    def test_multiline_fasta(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_text(">chr1\nACGT\nACGT\n")
        assert _genome_size_from_prep(str(fa)) == 8


# ---------------------------------------------------------------------------
# Empty class dict
# ---------------------------------------------------------------------------

def test_empty_class_dict_has_all_classes():
    d = _empty_class_dict()
    assert set(d.keys()) == set(CLASS_ORDER)
    assert all(v == 0 for v in d.values())


# ---------------------------------------------------------------------------
# Presence/absence mode integration
# ---------------------------------------------------------------------------

PREP_CONTENT = ">chr1\n" + "A" * 10000 + "\n"

GFF_SP1 = textwrap.dedent("""\
    chr1\tEarl_Grey\tLTR/Gypsy\t101\t200\t.\t+\t.\tID=f1_1;NAME=fam_shared
    chr1\tEarl_Grey\tDNA/TIR\t301\t500\t.\t-\t.\tID=f2_1;NAME=fam_unique_sp1
""")

GFF_SP2 = textwrap.dedent("""\
    chr1\tEarl_Grey\tLTR/Gypsy\t101\t200\t.\t+\t.\tID=f1_1;NAME=fam_shared
    chr1\tEarl_Grey\tSINE/tRNA\t301\t400\t.\t-\t.\tID=f3_1;NAME=fam_unique_sp2
""")


class TestPresenceAbsenceMode:
    def _run(self, tmp_path):
        g1 = tmp_path / "sp1.gff"
        g2 = tmp_path / "sp2.gff"
        p1 = tmp_path / "sp1.prep"
        p2 = tmp_path / "sp2.prep"
        g1.write_text(GFF_SP1)
        g2.write_text(GFF_SP2)
        p1.write_text(PREP_CONTENT)
        p2.write_text(PREP_CONTENT)
        return _presence_absence_mode(
            ["sp1", "sp2"],
            [str(g1), str(g2)],
            [str(p1), str(p2)],
        )

    def test_shared_family_counted_in_both(self, tmp_path):
        fam_data, _ = self._run(tmp_path)
        assert fam_data["sp1"]["shared_by_class"]["LTR"] == 1
        assert fam_data["sp2"]["shared_by_class"]["LTR"] == 1

    def test_unique_family_in_correct_species(self, tmp_path):
        fam_data, _ = self._run(tmp_path)
        assert fam_data["sp1"]["unique_by_class"]["DNA"] == 1
        assert fam_data["sp2"]["unique_by_class"]["SINE"] == 1

    def test_unique_not_in_wrong_species(self, tmp_path):
        fam_data, _ = self._run(tmp_path)
        assert fam_data["sp2"]["unique_by_class"].get("DNA", 0) == 0

    def test_coverage_shared_correct(self, tmp_path):
        _, cov_data = self._run(tmp_path)
        # fam_shared: bp = 200 - 101 + 1 = 100 in both species
        assert cov_data["sp1"]["shared_bp_by_class"]["LTR"] == 100
        assert cov_data["sp2"]["shared_bp_by_class"]["LTR"] == 100

    def test_genome_size_recorded(self, tmp_path):
        _, cov_data = self._run(tmp_path)
        assert cov_data["sp1"]["genome_size"] == 10000

    def test_no_cross_species_pollution(self, tmp_path):
        fam_data, _ = self._run(tmp_path)
        # sp1 unique DNA should not appear as shared
        assert fam_data["sp1"]["shared_by_class"].get("DNA", 0) == 0


# ---------------------------------------------------------------------------
# CLASS_ORDER and palette completeness
# ---------------------------------------------------------------------------

def test_all_classes_have_palette():
    for cls in CLASS_ORDER:
        assert cls in CLASS_PALETTE, f"Missing palette entry for {cls}"
        assert cls in LIGHT_PALETTE, f"Missing light palette entry for {cls}"


def test_other_label_contains_keywords():
    assert "Simple" in OTHER_LABEL or "Sat" in OTHER_LABEL
