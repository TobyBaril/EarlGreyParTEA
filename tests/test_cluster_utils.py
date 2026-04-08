"""tests/test_cluster_utils.py — unit tests for cluster_utils.parse_clstr."""

import os
import sys
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from cluster_utils import parse_clstr


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_clstr(tmp_path, content, filename="test.clstr"):
    p = tmp_path / filename
    p.write_text(textwrap.dedent(content))
    return str(p)


# ---------------------------------------------------------------------------
# Basic parsing
# ---------------------------------------------------------------------------

SIMPLE_CLSTR = """\
>Cluster 0
0\t500nt, >sp1_RND-1_FAMILY-1#LTR/Gypsy... *
1\t480nt, >sp2_RND-1_FAMILY-1#LTR/Gypsy... at +/97.00%
>Cluster 1
0\t400nt, >sp1_RND-1_FAMILY-2#LINE/R2... *
>Cluster 2
0\t300nt, >REPMASKER_fungi_known_element#DNA... *
"""


class TestParseClstrBasic:
    def test_returns_four_values(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        result = parse_clstr(path, ["sp1", "sp2"])
        assert len(result) == 4

    def test_species_to_clusters_both_in_cluster0(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        s2c, _, _, _ = parse_clstr(path, ["sp1", "sp2"])
        assert 0 in s2c["sp1"]
        assert 0 in s2c["sp2"]

    def test_species_to_clusters_only_sp1_in_cluster1(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        s2c, _, _, _ = parse_clstr(path, ["sp1", "sp2"])
        assert 1 in s2c["sp1"]
        assert 1 not in s2c["sp2"]

    def test_existing_clusters_repmasker(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        _, existing, _, _ = parse_clstr(path, ["sp1", "sp2"])
        assert 2 in existing
        assert 0 not in existing
        assert 1 not in existing

    def test_cluster_to_species_shared(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        _, _, c2s, _ = parse_clstr(path, ["sp1", "sp2"])
        assert c2s[0] == frozenset({"sp1", "sp2"})
        assert c2s[1] == frozenset({"sp1"})

    def test_rep_name_stored(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        _, _, _, r2c = parse_clstr(path, ["sp1", "sp2"])
        # Representative of cluster 0 is sp1_RND-1_FAMILY-1#LTR/Gypsy
        # Stored with full name and bare name (no class)
        assert 0 in r2c.values()

    def test_bare_name_without_species_prefix_stored(self, tmp_path):
        path = _write_clstr(tmp_path, SIMPLE_CLSTR)
        _, _, _, r2c = parse_clstr(path, ["sp1", "sp2"])
        # bare_no_prefix for sp1_RND-1_FAMILY-1#LTR/Gypsy is RND-1_FAMILY-1
        assert "RND-1_FAMILY-1" in r2c
        assert r2c["RND-1_FAMILY-1"] == 0


# ---------------------------------------------------------------------------
# Species prefix substring safety
# ---------------------------------------------------------------------------

SUBSTRING_CLSTR = """\
>Cluster 0
0\t400nt, >Sp1_RND-1_FAM-A#DNA/TIR... *
>Cluster 1
0\t400nt, >Sp10_RND-1_FAM-B#DNA/TIR... *
"""


def test_species_prefix_longest_first(tmp_path):
    """Sp1 must not match sequences prefixed with Sp10."""
    path = _write_clstr(tmp_path, SUBSTRING_CLSTR)
    s2c, _, _, _ = parse_clstr(path, ["Sp1", "Sp10"])
    # Cluster 0 belongs to Sp1 only
    assert 0 in s2c["Sp1"]
    assert 0 not in s2c["Sp10"]
    # Cluster 1 belongs to Sp10 only
    assert 1 in s2c["Sp10"]
    assert 1 not in s2c["Sp1"]


# ---------------------------------------------------------------------------
# CUSTOM library prefix
# ---------------------------------------------------------------------------

CUSTOM_CLSTR = """\
>Cluster 0
0\t300nt, >CUSTOM_my_element#LTR/Copia... *
1\t290nt, >sp1_RND-1_FAM-1#LTR/Copia... at +/95.00%
"""


def test_custom_prefix_goes_to_existing(tmp_path):
    path = _write_clstr(tmp_path, CUSTOM_CLSTR)
    _, existing, _, _ = parse_clstr(path, ["sp1"])
    assert 0 in existing


# ---------------------------------------------------------------------------
# Empty/sentinel file
# ---------------------------------------------------------------------------

def test_empty_clstr(tmp_path):
    path = _write_clstr(tmp_path, "")
    s2c, existing, c2s, r2c = parse_clstr(path, ["sp1", "sp2"])
    assert s2c == {"sp1": set(), "sp2": set()}
    assert existing == set()
    assert c2s == {}
    assert r2c == {}


# ---------------------------------------------------------------------------
# Species with no sequences
# ---------------------------------------------------------------------------

NO_SP2_CLSTR = """\
>Cluster 0
0\t400nt, >sp1_RND-1_FAM-A#LTR/Gypsy... *
"""


def test_species_not_in_clstr_has_empty_set(tmp_path):
    path = _write_clstr(tmp_path, NO_SP2_CLSTR)
    s2c, _, _, _ = parse_clstr(path, ["sp1", "sp2"])
    assert s2c["sp2"] == set()
    assert len(s2c["sp1"]) == 1


# ---------------------------------------------------------------------------
# Multispecies cluster → frozenset
# ---------------------------------------------------------------------------

MULTI_CLSTR = """\
>Cluster 0
0\t500nt, >A_FAM1#LTR... *
1\t490nt, >B_FAM1#LTR... at +/98%
2\t480nt, >C_FAM1#LTR... at -/96%
"""


def test_three_species_cluster(tmp_path):
    path = _write_clstr(tmp_path, MULTI_CLSTR)
    _, _, c2s, _ = parse_clstr(path, ["A", "B", "C"])
    assert c2s[0] == frozenset({"A", "B", "C"})
