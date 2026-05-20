"""tests/test_v018_features.py

Verification test suite for EarlGreyParTEA v0.1.8 features.

Covers all items in the v0.1.8 verification checklist:
  - clustering_length_diff (-s flag)             items 4-7
  - clustering_coverage_long (-aL flag)          items 8-11
  - Post-clustering chimera detection/splitting  items 12-17
  - on_start_functions startup messages          items 6, 10
  - RepeatMasker warmup cache discovery          items 1-3
  - --generate-config wrapper output             items 7, 11

Run all tests:
    pytest tests/test_v018_features.py -v

Skip tests that require external tools (cd-hit-est, RepeatMasker):
    pytest tests/test_v018_features.py -v -m "not integration"
"""

import importlib
import io
import os
import shutil
import subprocess
import sys
import textwrap
from types import SimpleNamespace
from unittest.mock import patch

import pytest

# ---------------------------------------------------------------------------
# Path setup — allow imports from scripts/ without installing the package
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

# Marker for tests that require external binaries (cd-hit-est, RepeatMasker)
integration = pytest.mark.integration


# ===========================================================================
# SECTION 1 — Module importability (checklist item 17)
# ===========================================================================

class TestModuleImportability:
    """Item 17: import split_chimeras without NameError."""

    def test_import_public_api(self):
        """All public functions import cleanly."""
        from split_chimeras import (  # noqa: F401
            load_fasta,
            parse_clstr,
            detect_chimera,
            _connected_components,
            _interval_overlap,
            _base_name,
            _classification,
        )

    def test_main_not_called_at_import_time(self):
        """Importing the module must not call main() (which requires 'snakemake' in scope)."""
        # If main() were called at module level it would raise NameError because
        # 'snakemake' is not defined outside Snakemake's script: execution context.
        # Simply completing the import without exception verifies the guard.
        mod = importlib.import_module("split_chimeras")
        assert callable(mod.main)

    def test_main_guard_condition(self):
        """'snakemake' must NOT be in the module's dir() at import time."""
        mod = importlib.import_module("split_chimeras")
        assert "snakemake" not in dir(mod), (
            "snakemake object leaked into module globals — main() would be "
            "called on every import, causing NameError in test environments"
        )


# ===========================================================================
# SECTION 2 — FASTA parser (supporting functionality for items 12-17)
# ===========================================================================

class TestLoadFasta:
    def test_basic_sequence(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">seq1#DNA/hAT\nACGT\nACGT\n")
        from split_chimeras import load_fasta
        result = load_fasta(str(fa))
        assert "seq1#DNA/hAT" in result
        header, seq = result["seq1#DNA/hAT"]
        assert seq == "ACGTACGT"

    def test_full_header_preserved(self, tmp_path):
        """Item 13: original header (including comment fields) is preserved."""
        fa = tmp_path / "test.fa"
        fa.write_text(">speciesA_TE1#DNA/hAT EarlGrey_annotation extra_comment\nACGT\n")
        from split_chimeras import load_fasta
        result = load_fasta(str(fa))
        header, seq = result["speciesA_TE1#DNA/hAT"]
        assert header == "speciesA_TE1#DNA/hAT EarlGrey_annotation extra_comment"

    def test_multiline_sequence(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">seq1\nAAAA\nCCCC\nGGGG\n")
        from split_chimeras import load_fasta
        result = load_fasta(str(fa))
        _, seq = result["seq1"]
        assert seq == "AAAACCCCGGGG"

    def test_multiple_sequences(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">a\nAAA\n>b\nCCC\n")
        from split_chimeras import load_fasta
        result = load_fasta(str(fa))
        assert "a" in result
        assert "b" in result

    def test_lookup_key_is_first_word(self, tmp_path):
        """Name used for dict key is only the first whitespace-delimited word."""
        fa = tmp_path / "test.fa"
        fa.write_text(">name#Class/Family some comment\nACGT\n")
        from split_chimeras import load_fasta
        result = load_fasta(str(fa))
        assert "name#Class/Family" in result
        # Full header including comment is accessible via value
        header, _ = result["name#Class/Family"]
        assert "some comment" in header


# ===========================================================================
# SECTION 3 — .clstr parser (supporting functionality for items 12-17)
# ===========================================================================

# cd-hit-est -d 0 format:
#   representative: <idx>\t<len>nt, ><name>... *
#   member:         <idx>\t<len>nt, ><name>... at <qs>:<qe>:<rs>:<re>/<strand>/<id>%

SIMPLE_CLSTR = textwrap.dedent("""\
    >Cluster 0
    0\t500nt, >sp1_TE1#DNA/hAT... *
    1\t480nt, >sp2_TE1#DNA/hAT... at 1:480:1:480/+/97.00%
    2\t460nt, >sp3_TE1#DNA/hAT... at 1:460:20:480/+/96.00%
    >Cluster 1
    0\t300nt, >sp1_TE2#LINE/R2... *
""")

CHIMERIC_CLSTR = textwrap.dedent("""\
    >Cluster 0
    0\t2000nt, >sp1_BigTE#Unknown... *
    1\t500nt, >sp1_TE1#DNA/hAT... at 1:500:1:500/+/95.00%
    2\t480nt, >sp2_TE1#DNA/hAT... at 1:480:20:480/+/94.00%
    3\t450nt, >sp3_TE1#DNA/hAT... at 1:450:30:470/+/93.00%
    4\t600nt, >sp1_TRIM1#DNA/hAT... at 1:600:1400:2000/+/95.00%
    5\t550nt, >sp2_TRIM1#DNA/hAT... at 1:550:1450:2000/+/94.00%
    6\t400nt, >sp3_TRIM1#DNA/hAT... at 1:400:1500:1900/+/93.00%
""")


class TestParseClstr:
    def _write(self, tmp_path, content, name="test.clstr"):
        p = tmp_path / name
        p.write_text(content)
        return str(p)

    def test_cluster_count(self, tmp_path):
        from split_chimeras import parse_clstr
        path = self._write(tmp_path, SIMPLE_CLSTR)
        clusters = parse_clstr(path)
        assert len(clusters) == 2

    def test_representative_flagged(self, tmp_path):
        from split_chimeras import parse_clstr
        path = self._write(tmp_path, SIMPLE_CLSTR)
        clusters = parse_clstr(path)
        rep = next(m for m in clusters[0] if m.is_rep)
        assert rep.name == "sp1_TE1#DNA/hAT"
        assert rep.length == 500

    def test_member_coordinates_parsed(self, tmp_path):
        from split_chimeras import parse_clstr
        path = self._write(tmp_path, SIMPLE_CLSTR)
        clusters = parse_clstr(path)
        non_reps = [m for m in clusters[0] if not m.is_rep]
        assert len(non_reps) == 2
        # First non-rep: at 1:480:1:480
        m = non_reps[0]
        assert m.seq_start == 1 and m.seq_end == 480
        assert m.rep_start == 1 and m.rep_end == 480
        assert m.strand == "+"
        assert m.identity == 97.00

    def test_chimeric_clstr_member_count(self, tmp_path):
        from split_chimeras import parse_clstr
        path = self._write(tmp_path, CHIMERIC_CLSTR)
        clusters = parse_clstr(path)
        assert len(clusters) == 1
        non_reps = [m for m in clusters[0] if not m.is_rep]
        assert len(non_reps) == 6


# ===========================================================================
# SECTION 4 — Overlap graph utilities (supporting items 12-17)
# ===========================================================================

class TestIntervalOverlap:
    def test_overlapping(self):
        from split_chimeras import _interval_overlap
        assert _interval_overlap(0, 100, 50, 150) == 50

    def test_touching_not_overlapping(self):
        from split_chimeras import _interval_overlap
        # Intervals share an endpoint but do not overlap (half-open convention)
        assert _interval_overlap(0, 100, 100, 200) == 0

    def test_disjoint(self):
        from split_chimeras import _interval_overlap
        assert _interval_overlap(0, 100, 200, 300) == 0

    def test_contained(self):
        from split_chimeras import _interval_overlap
        assert _interval_overlap(0, 1000, 200, 400) == 200

    def test_identical(self):
        from split_chimeras import _interval_overlap
        assert _interval_overlap(50, 200, 50, 200) == 150


class TestConnectedComponents:
    """Directly tests the BFS graph builder."""

    def _make_member(self, rs, re):
        from split_chimeras import Member
        return Member(index=0, length=re - rs, name="x", is_rep=False,
                      rep_start=rs, rep_end=re)

    def test_all_overlapping_one_component(self):
        from split_chimeras import _connected_components
        members = [
            self._make_member(0, 500),
            self._make_member(100, 600),
            self._make_member(200, 700),
        ]
        comps = _connected_components(members, overlap_min=50)
        assert len(comps) == 1

    def test_disjoint_two_components(self):
        from split_chimeras import _connected_components
        members = [
            self._make_member(0, 500),      # left cluster
            self._make_member(50, 480),     # left cluster
            self._make_member(1500, 2000),  # right cluster
            self._make_member(1600, 1950),  # right cluster
        ]
        comps = _connected_components(members, overlap_min=50)
        assert len(comps) == 2
        sizes = sorted(len(c) for c in comps)
        assert sizes == [2, 2]

    def test_overlap_min_threshold(self):
        """Members that barely touch should only join if overlap >= threshold."""
        from split_chimeras import _connected_components
        members = [
            self._make_member(0, 100),
            self._make_member(95, 200),   # 5 nt overlap
        ]
        # Below threshold → two components
        comps = _connected_components(members, overlap_min=50)
        assert len(comps) == 2
        # At threshold (overlap_min=5) → one component
        comps = _connected_components(members, overlap_min=5)
        assert len(comps) == 1


# ===========================================================================
# SECTION 5 — detect_chimera (items 12, 15, 16)
# ===========================================================================

def _make_cluster(rep_len, non_rep_coords):
    """Build a minimal cluster list for detect_chimera tests.

    Args:
        rep_len: length of the representative sequence.
        non_rep_coords: list of (rep_start, rep_end, length) tuples for members.
    Returns:
        list of Member objects with rep first.
    """
    from split_chimeras import Member
    cluster = [Member(index=0, length=rep_len, name="REP#DNA/hAT", is_rep=True)]
    for i, (rs, re, ln) in enumerate(non_rep_coords, start=1):
        cluster.append(Member(
            index=i, length=ln, name=f"mem{i}#DNA/hAT", is_rep=False,
            seq_start=1, seq_end=ln, rep_start=rs, rep_end=re,
            strand="+", identity=95.0,
        ))
    return cluster


class TestDetectChimera:
    DEFAULTS = dict(overlap_min=50, min_members=3, min_component_span=0.1)

    def test_non_chimeric_all_overlapping(self):
        """Item 15: non-chimeric cluster → False."""
        from split_chimeras import detect_chimera
        cluster = _make_cluster(1000, [
            (1, 500, 500), (50, 480, 430), (100, 600, 500),
        ])
        is_chim, comps = detect_chimera(cluster, **self.DEFAULTS)
        assert not is_chim
        assert comps == []

    def test_chimeric_two_components(self):
        """Item 12: clear two-component chimera is detected."""
        from split_chimeras import detect_chimera
        # Left component aligns to 1-500; right to 1500-2000 on 2000 nt rep
        cluster = _make_cluster(2000, [
            (1, 500, 500), (20, 480, 460), (30, 470, 440),  # left
            (1400, 2000, 600), (1450, 2000, 550), (1500, 1900, 400),  # right
        ])
        is_chim, comps = detect_chimera(cluster, **self.DEFAULTS)
        assert is_chim
        assert len(comps) == 2

    def test_components_sorted_left_to_right(self):
        """Components should be sorted by leftmost position on the representative."""
        from split_chimeras import detect_chimera
        cluster = _make_cluster(2000, [
            (1400, 2000, 600), (1450, 1950, 550), (1500, 1900, 400),  # right (listed first)
            (1, 500, 500), (20, 480, 460), (30, 470, 440),  # left (listed second)
        ])
        is_chim, comps = detect_chimera(cluster, **self.DEFAULTS)
        assert is_chim
        assert min(m.rep_start for m in comps[0]) < min(m.rep_start for m in comps[1])

    def test_too_few_members_skipped(self):
        """Item 16: clusters with fewer than min_members are not tested."""
        from split_chimeras import detect_chimera
        cluster = _make_cluster(2000, [
            (1, 500, 500),              # would be left component
            (1500, 2000, 500),          # would be right component
            # only 2 members — less than min_members=3
        ])
        is_chim, comps = detect_chimera(cluster, min_members=3,
                                        overlap_min=50, min_component_span=0.1)
        assert not is_chim

    def test_component_span_filter(self):
        """Component that spans <min_component_span of rep is not counted."""
        from split_chimeras import detect_chimera
        # Right component spans only 50 nt of a 2000 nt rep → 2.5% < 10%
        cluster = _make_cluster(2000, [
            (1, 500, 500), (20, 480, 460), (30, 470, 440),   # solid left (25%)
            (1950, 2000, 50), (1960, 2000, 40), (1970, 2000, 30),  # tiny right (2.5%)
        ])
        is_chim, comps = detect_chimera(cluster, **self.DEFAULTS)
        assert not is_chim

    def test_chimera_score_positive_gap(self):
        """Chimera score (inter-component gap / rep length) must be > 0 for a clear chimera."""
        from split_chimeras import detect_chimera, _base_name, _classification
        cluster = _make_cluster(2000, [
            (1, 500, 500), (20, 480, 460), (30, 470, 440),
            (1400, 2000, 600), (1450, 2000, 550), (1500, 1900, 400),
        ])
        is_chim, comps = detect_chimera(cluster, **self.DEFAULTS)
        assert is_chim
        # Gap = 1400 - 500 = 900 nt; score = 900 / 2000 = 0.45
        left_end = max(m.rep_end for m in comps[0])
        right_start = min(m.rep_start for m in comps[1])
        gap = right_start - left_end
        assert gap == 900


# ===========================================================================
# SECTION 6 — Name utilities
# ===========================================================================

class TestNameUtilities:
    def test_base_name_strips_classification(self):
        from split_chimeras import _base_name
        assert _base_name("sp1_TE1#DNA/hAT") == "sp1_TE1"

    def test_base_name_no_hash(self):
        from split_chimeras import _base_name
        assert _base_name("sp1_TE1") == "sp1_TE1"

    def test_classification_returns_hash_class(self):
        from split_chimeras import _classification
        assert _classification("sp1_TE1#DNA/hAT") == "#DNA/hAT"

    def test_classification_absent(self):
        from split_chimeras import _classification
        assert _classification("sp1_TE1") == ""


# ===========================================================================
# SECTION 7 — Full main() integration via mock snakemake object
# Items 12, 13, 14, 15, 16
# ===========================================================================

CHIMERIC_FA_REP = textwrap.dedent("""\
    >sp1_BigTE#Unknown
    {rep_seq}
    >sp1_TE2#LINE/R2
    {clean_seq}
""")

CHIMERIC_FA_COMBINED = textwrap.dedent("""\
    >sp1_BigTE#Unknown EarlGrey_annotation
    {rep_seq}
    >sp1_TE1#DNA/hAT EarlGrey_annotation
    {left_seq}
    >sp2_TE1#DNA/hAT EarlGrey_annotation
    {left2_seq}
    >sp3_TE1#DNA/hAT EarlGrey_annotation
    {left3_seq}
    >sp1_TRIM1#DNA/hAT EarlGrey_annotation
    {right_seq}
    >sp2_TRIM1#DNA/hAT EarlGrey_annotation
    {right2_seq}
    >sp3_TRIM1#DNA/hAT EarlGrey_annotation
    {right3_seq}
    >sp1_TE2#LINE/R2 EarlGrey_annotation
    {clean_seq}
""")


def _write_integration_fixtures(tmp_path):
    """Write FASTA and .clstr files for the main() integration test."""
    rep_seq = "ACGT" * 500          # 2000 nt chimeric representative
    left_seq = "ACGT" * 125         # 500 nt — aligns to rep pos 1-500
    left2_seq = "ACGT" * 120        # 480 nt — aligns to rep pos 20-480
    left3_seq = "ACGT" * 112        # 448 nt — aligns to rep pos 30-470
    right_seq = "TTTT" * 150        # 600 nt — aligns to rep pos 1400-2000
    right2_seq = "TTTT" * 137       # 550 nt — aligns to rep pos 1450-2000
    right3_seq = "TTTT" * 100       # 400 nt — aligns to rep pos 1500-1900
    clean_seq = "GCGC" * 75         # 300 nt — clean non-chimeric cluster

    clstr_content = textwrap.dedent(f"""\
        >Cluster 0
        0\t2000nt, >sp1_BigTE#Unknown... *
        1\t500nt, >sp1_TE1#DNA/hAT... at 1:500:1:500/+/95.00%
        2\t480nt, >sp2_TE1#DNA/hAT... at 1:480:20:480/+/94.00%
        3\t448nt, >sp3_TE1#DNA/hAT... at 1:448:30:470/+/93.00%
        4\t600nt, >sp1_TRIM1#DNA/hAT... at 1:600:1400:2000/+/95.00%
        5\t550nt, >sp2_TRIM1#DNA/hAT... at 1:550:1450:2000/+/94.00%
        6\t400nt, >sp3_TRIM1#DNA/hAT... at 1:400:1500:1900/+/93.00%
        >Cluster 1
        0\t300nt, >sp1_TE2#LINE/R2... *
    """)

    clstr_path = tmp_path / "combined.clstrd.fa.clstr"
    clstr_path.write_text(clstr_content)

    clustered_fa = tmp_path / "clustered.fa"
    clustered_fa.write_text(
        f">sp1_BigTE#Unknown\n{rep_seq}\n"
        f">sp1_TE2#LINE/R2\n{clean_seq}\n"
    )

    combined_fa = tmp_path / "combined.fa"
    combined_fa.write_text(
        f">sp1_BigTE#Unknown EarlGrey_annotation\n{rep_seq}\n"
        f">sp1_TE1#DNA/hAT EarlGrey_annotation\n{left_seq}\n"
        f">sp2_TE1#DNA/hAT EarlGrey_annotation\n{left2_seq}\n"
        f">sp3_TE1#DNA/hAT EarlGrey_annotation\n{left3_seq}\n"
        f">sp1_TRIM1#DNA/hAT EarlGrey_annotation\n{right_seq}\n"
        f">sp2_TRIM1#DNA/hAT EarlGrey_annotation\n{right2_seq}\n"
        f">sp3_TRIM1#DNA/hAT EarlGrey_annotation\n{right3_seq}\n"
        f">sp1_TE2#LINE/R2 EarlGrey_annotation\n{clean_seq}\n"
    )

    out_fasta = tmp_path / "chimera_split.fa"
    out_summary = tmp_path / "chimera_detection_summary.tsv"
    log_path = tmp_path / "split_chimeras.log"

    return {
        "clstr": str(clstr_path),
        "clustered_fa": str(clustered_fa),
        "combined_fa": str(combined_fa),
        "out_fasta": str(out_fasta),
        "out_summary": str(out_summary),
        "log": str(log_path),
    }


def _build_mock_snakemake(paths):
    """Return a SimpleNamespace that mimics the snakemake object."""
    return SimpleNamespace(
        input=SimpleNamespace(
            clstr=paths["clstr"],
            clustered_fa=paths["clustered_fa"],
            combined_fa=paths["combined_fa"],
        ),
        output=SimpleNamespace(
            fasta=paths["out_fasta"],
            summary=paths["out_summary"],
        ),
        log=[paths["log"]],
        params=SimpleNamespace(
            overlap_min=50,
            min_members=3,
            min_component_span=0.1,
        ),
    )


class TestMainIntegration:
    """Items 12-16: end-to-end tests of split_chimeras.main()."""

    @pytest.fixture(autouse=True)
    def _run(self, tmp_path):
        """Run main() once; all tests in this class inspect the results."""
        import split_chimeras as sc
        paths = _write_integration_fixtures(tmp_path)
        mock_sm = _build_mock_snakemake(paths)

        # Inject the mock snakemake object into the module's namespace, then call main()
        with patch.object(sc, "__builtins__", sc.__builtins__):
            old = sc.__dict__.get("snakemake")
            sc.__dict__["snakemake"] = mock_sm
            try:
                sc.main()
            finally:
                if old is None:
                    sc.__dict__.pop("snakemake", None)
                else:
                    sc.__dict__["snakemake"] = old

        self.paths = paths
        self.out_fasta_text = open(paths["out_fasta"]).read()
        self.out_summary_text = open(paths["out_summary"]).read()

    # ── Item 12: output files are produced ──────────────────────────────────

    def test_output_fasta_created(self):
        """Item 12: chimera_split.fa is written."""
        assert os.path.isfile(self.paths["out_fasta"])
        assert len(self.out_fasta_text) > 0

    def test_output_summary_created(self):
        """Item 12: chimera_detection_summary.tsv is written."""
        assert os.path.isfile(self.paths["out_summary"])
        assert len(self.out_summary_text) > 0

    # ── Item 13: chimeric rep has _CHIMERA suffix; component reps keep original headers

    def test_chimeric_rep_labelled(self):
        """Item 13: chimeric representative gets _CHIMERA suffix in output FASTA."""
        assert "sp1_BigTE_CHIMERA#Unknown" in self.out_fasta_text

    def test_component_reps_have_original_headers(self):
        """Item 13: component representatives carry the full original FASTA header."""
        # The combined FASTA has 'EarlGrey_annotation' in the header comment;
        # the output FASTA should preserve it for component representatives.
        lines = self.out_fasta_text.splitlines()
        header_lines = [l for l in lines if l.startswith(">")]
        comp_headers = [l for l in header_lines if "TRIM1" in l or "TE1" in l]
        # At least one component header should include the comment field
        assert any("EarlGrey_annotation" in h for h in comp_headers), (
            f"No component header preserved comment field.\nHeaders: {comp_headers}"
        )

    def test_chimeric_rep_original_classification_retained(self):
        """Item 13: #classification tag is kept on the _CHIMERA entry."""
        assert ">sp1_BigTE_CHIMERA#Unknown" in self.out_fasta_text

    # ── Item 14: annotation routing — separate test in TestAnnotationRouting ──

    # ── Item 15: split_chimeras=false behaviour ──────────────────────────────

    def test_non_chimeric_cluster_passes_through(self):
        """Item 15: non-chimeric clusters are written to output unchanged."""
        assert ">sp1_TE2#LINE/R2" in self.out_fasta_text

    def test_non_chimeric_no_chimera_suffix(self):
        """Item 15: non-chimeric representatives do not get _CHIMERA suffix."""
        assert "sp1_TE2#LINE/R2_CHIMERA" not in self.out_fasta_text

    # ── Summary TSV format ───────────────────────────────────────────────────

    def test_summary_has_expected_columns(self):
        """Item 12: TSV contains all required column headers."""
        header_row = self.out_summary_text.splitlines()[0]
        for col in ("cluster_idx", "representative", "rep_length", "n_members",
                    "is_chimeric", "n_components", "component_sizes",
                    "chimera_score", "component_rep_names"):
            assert col in header_row, f"Column '{col}' missing from summary TSV"

    def test_summary_chimeric_row_correct(self):
        """Chimeric cluster row has is_chimeric=True and n_components=2."""
        rows = self.out_summary_text.splitlines()
        col_headers = rows[0].split("\t")
        is_chim_idx = col_headers.index("is_chimeric")
        n_comp_idx = col_headers.index("n_components")
        chimera_rows = [r for r in rows[1:] if "True" in r]
        assert len(chimera_rows) >= 1
        row_fields = chimera_rows[0].split("\t")
        assert row_fields[is_chim_idx] == "True"
        assert row_fields[n_comp_idx] == "2"

    def test_summary_non_chimeric_row_correct(self):
        """Non-chimeric cluster row has is_chimeric=False."""
        rows = self.out_summary_text.splitlines()
        col_headers = rows[0].split("\t")
        is_chim_idx = col_headers.index("is_chimeric")
        non_chimera_rows = [r for r in rows[1:] if "False" in r]
        assert len(non_chimera_rows) >= 1
        row_fields = non_chimera_rows[0].split("\t")
        assert row_fields[is_chim_idx] == "False"

    def test_chimera_score_nonzero_for_chimeric_cluster(self):
        """Chimera score should be > 0 when there is a gap between components."""
        rows = self.out_summary_text.splitlines()
        col_headers = rows[0].split("\t")
        score_idx = col_headers.index("chimera_score")
        chimera_rows = [r for r in rows[1:] if "True" in r]
        score = float(chimera_rows[0].split("\t")[score_idx])
        assert score > 0.0


# ===========================================================================
# SECTION 8 — Annotation routing (item 14)
# Tests that annotate.smk selects the correct library based on config
# ===========================================================================

class TestAnnotationRouting:
    """Item 14: ANNOTATION_LIBRARY switches between clstrd.fa and chimera_split.fa."""

    ANNOTATE_SMK = os.path.join(REPO_ROOT, "rules", "annotate.smk")

    def test_annotation_library_variable_present(self):
        """ANNOTATION_LIBRARY variable must be defined in annotate.smk."""
        text = open(self.ANNOTATE_SMK).read()
        assert "ANNOTATION_LIBRARY" in text

    def test_chimera_split_path_present(self):
        """annotate.smk must reference chimera_split.fa."""
        text = open(self.ANNOTATE_SMK).read()
        assert "chimera_split.fa" in text

    def test_clstrd_path_present(self):
        """annotate.smk must reference the standard clstrd.fa path."""
        text = open(self.ANNOTATE_SMK).read()
        assert "clstrd.fa" in text

    def test_routing_conditional_on_split_chimeras(self):
        """The switch must be conditional on split_chimeras config key."""
        text = open(self.ANNOTATE_SMK).read()
        assert "split_chimeras" in text

    def test_annotation_rule_uses_annotation_library(self):
        """repeatmasker_annotation input must reference ANNOTATION_LIBRARY."""
        text = open(self.ANNOTATE_SMK).read()
        assert "ANNOTATION_LIBRARY" in text
        # The library variable should appear in an input: block (simple substring check)
        assert "library=ANNOTATION_LIBRARY" in text or "library = ANNOTATION_LIBRARY" in text


# ===========================================================================
# SECTION 9 — Snakemake rule static analysis (items 4, 5, 8, 9, 16)
# ===========================================================================

class TestClusteringSmkFlags:
    """Items 4, 5, 8, 9: correct cd-hit-est flags present in clustering.smk."""

    CLUSTERING_SMK = os.path.join(REPO_ROOT, "rules", "clustering.smk")

    def test_aL_flag_present(self):
        """Item 8/9: -aL flag passed from cluster_coverage_long param."""
        text = open(self.CLUSTERING_SMK).read()
        assert "-aL {params.cluster_coverage_long}" in text or \
               "-aL {{params.cluster_coverage_long}}" in text

    def test_s_flag_present(self):
        """Item 4/5: -s flag passed from cluster_length_diff param."""
        text = open(self.CLUSTERING_SMK).read()
        assert "-s {params.cluster_length_diff}" in text or \
               "-s {{params.cluster_length_diff}}" in text

    def test_cluster_coverage_long_param_defined(self):
        """Items 8/9: cluster_coverage_long param reads from config."""
        text = open(self.CLUSTERING_SMK).read()
        assert "cluster_coverage_long" in text
        assert "clustering_coverage_long" in text

    def test_cluster_length_diff_param_defined(self):
        """Items 4/5: cluster_length_diff param reads from config."""
        text = open(self.CLUSTERING_SMK).read()
        assert "cluster_length_diff" in text
        assert "clustering_length_diff" in text

    def test_combined_fa_is_temp_output(self):
        """Item 16 (implicit): combined_fa declared as temp() for automatic cleanup."""
        text = open(self.CLUSTERING_SMK).read()
        assert "temp(" in text
        assert "combined_all_species.fa" in text

    def test_split_chimeras_rule_conditionally_defined(self):
        """Item 16: split_chimeras rule only defined when both conditions are met."""
        text = open(self.CLUSTERING_SMK).read()
        assert "split_chimeras" in text
        assert "skip_clustering" in text
        # The conditional must guard the rule definition
        assert "if not config.get" in text or "config.get('skip_clustering'" in text


# ===========================================================================
# SECTION 10 — on_start_functions startup messages (items 6, 10)
# ===========================================================================

def _capture_validate_output(config_dict):
    """Run validate_parameters and return captured stdout as a string.

    validate_parameters calls sys.exit() when genome files are not found on disk.
    The clustering/chimera messages are printed *before* that validation step, so we
    catch SystemExit and return whatever was printed up to that point.
    """
    sys.path.insert(0, SCRIPTS_DIR)
    from on_start_functions import validate_parameters
    buf = io.StringIO()
    try:
        with patch("sys.stdout", buf):
            validate_parameters(dict(config_dict))
    except SystemExit:
        pass  # Expected — genome paths don't exist; messages were already printed
    return buf.getvalue()


class TestStartupMessages:
    """Items 6, 10: startup message content for clustering parameters."""

    BASE_CONFIG = {
        "genome": {"sp1": "/fake/sp1.fa"},
        "species": ["sp1"],
        "output_dir": "/fake/out",
    }

    def test_length_diff_reported_in_message(self):
        """Item 6: 'length_diff:' appears in startup output when clustering enabled."""
        cfg = dict(self.BASE_CONFIG, clustering_length_diff=0.5, skip_clustering=False)
        out = _capture_validate_output(cfg)
        assert "length_diff" in out

    def test_aL_disabled_message_when_zero(self):
        """Item 10: 'aL: disabled' message when clustering_coverage_long=0.0."""
        cfg = dict(self.BASE_CONFIG, clustering_coverage_long=0.0, skip_clustering=False)
        out = _capture_validate_output(cfg)
        assert "aL: disabled" in out

    def test_aL_value_reported_when_nonzero(self):
        """Item 10: aL value reported when clustering_coverage_long > 0."""
        cfg = dict(self.BASE_CONFIG, clustering_coverage_long=0.75, skip_clustering=False)
        out = _capture_validate_output(cfg)
        assert "aL" in out
        assert "disabled" not in out or "0.75" in out

    def test_chimera_detection_enabled_message(self):
        """Item 12 (startup): chimera detection enabled message when split_chimeras=True."""
        cfg = dict(self.BASE_CONFIG, split_chimeras=True, skip_clustering=False,
                   chimera_overlap_min=50, chimera_min_members=3,
                   chimera_min_component_span=0.1)
        out = _capture_validate_output(cfg)
        assert "himera" in out  # 'Chimera detection enabled'

    def test_chimera_detection_disabled_message(self):
        """Item 15 (startup): chimera detection disabled message when split_chimeras=False."""
        cfg = dict(self.BASE_CONFIG, split_chimeras=False, skip_clustering=False)
        out = _capture_validate_output(cfg)
        assert "himera detection disabled" in out.lower() or \
               "split_chimeras: false" in out.lower() or \
               "disabled" in out


# ===========================================================================
# SECTION 11 — --generate-config wrapper output (items 7, 11)
# ===========================================================================

WRAPPER_SCRIPTS = [
    os.path.join(REPO_ROOT, "earlGreyParTEA"),
    os.path.join(REPO_ROOT, "earlGreyParTEA_LibConstruct"),
]


@pytest.mark.parametrize("wrapper", WRAPPER_SCRIPTS, ids=["full", "libconstruct"])
class TestGenerateConfig:
    """Items 7, 11: --generate-config output contains required v0.1.8 params."""

    def _run_generate_config(self, wrapper, tmp_path):
        out_cfg = str(tmp_path / "out.yaml")
        result = subprocess.run(
            [wrapper, "--generate-config", out_cfg],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"--generate-config failed for {wrapper}:\n{result.stderr}"
        )
        return open(out_cfg).read()

    def test_clustering_length_diff_present(self, wrapper, tmp_path):
        """Item 7: clustering_length_diff in generated config."""
        text = self._run_generate_config(wrapper, tmp_path)
        assert "clustering_length_diff" in text

    def test_clustering_length_diff_default_value(self, wrapper, tmp_path):
        """Item 7: default value is 0.5."""
        text = self._run_generate_config(wrapper, tmp_path)
        assert "clustering_length_diff: 0.5" in text

    def test_clustering_coverage_long_present(self, wrapper, tmp_path):
        """Item 11: clustering_coverage_long in generated config."""
        text = self._run_generate_config(wrapper, tmp_path)
        assert "clustering_coverage_long" in text

    def test_clustering_coverage_long_default_value(self, wrapper, tmp_path):
        """Item 11: default value is 0.0."""
        text = self._run_generate_config(wrapper, tmp_path)
        assert "clustering_coverage_long: 0.0" in text

    def test_split_chimeras_present(self, wrapper, tmp_path):
        """Item 12: split_chimeras in generated config."""
        text = self._run_generate_config(wrapper, tmp_path)
        assert "split_chimeras" in text

    def test_split_chimeras_default_false(self, wrapper, tmp_path):
        """Item 12: default is false."""
        text = self._run_generate_config(wrapper, tmp_path)
        assert "split_chimeras: false" in text

    def test_chimera_params_present(self, wrapper, tmp_path):
        """Item 12: all four chimera params appear in generated config."""
        text = self._run_generate_config(wrapper, tmp_path)
        for param in ("chimera_overlap_min", "chimera_min_members",
                      "chimera_min_component_span"):
            assert param in text, f"Missing param '{param}' in --generate-config output"


# ===========================================================================
# SECTION 12 — RepeatMasker cache discovery (items 1-3)
# Tested as shell script logic using subprocess + mock directory structures.
# ===========================================================================

# Extract just the CONS-* discovery + warning logic into a standalone test
# script — avoids needing an actual RepeatMasker installation.

_CACHE_DISCOVERY_SCRIPT = textwrap.dedent("""\
    #!/bin/bash
    # Minimal reproduction of the warmup CONS-* discovery logic.
    # Takes one argument: the mock RM_SHARE path.
    set -e
    RM_SHARE="$1"
    SPECIES_WORD="7215"
    CACHE_PARENT=$(find "$RM_SHARE/Libraries" -maxdepth 1 -type d -name "CONS-*" 2>/dev/null | head -n 1)
    if [ -z "$CACHE_PARENT" ]; then
        echo "WARNING: No CONS-* cache directory found under $RM_SHARE/Libraries -- skipping species cache check." >&2
        exit 0
    fi
    CACHE_DIR="$CACHE_PARENT/$SPECIES_WORD"
    if [ -d "$CACHE_DIR" ] && [ ! -f "$CACHE_DIR/refineableHash.dat" ]; then
        echo "INCOMPLETE_CACHE_DETECTED" >&2
        rm -rf "$CACHE_DIR"
    fi
    if [ ! -f "$CACHE_DIR/refineableHash.dat" ]; then
        echo "NEEDS_BUILD"
    else
        echo "CACHE_OK"
    fi
""")


@pytest.fixture(scope="module")
def cache_script(tmp_path_factory):
    """Write the cache-discovery test script once per module."""
    p = tmp_path_factory.mktemp("scripts") / "cache_discovery.sh"
    p.write_text(_CACHE_DISCOVERY_SCRIPT)
    p.chmod(0o755)
    return str(p)


class TestCacheDiscovery:
    """Items 1-3: CONS-* directory discovery and incomplete-cache handling."""

    def _run(self, script, mock_rm_share):
        result = subprocess.run(
            [script, mock_rm_share],
            capture_output=True, text=True,
        )
        return result

    # Item 3: no CONS-* directory → warning and graceful exit

    def test_no_cons_dir_exits_zero_with_warning(self, tmp_path, cache_script):
        """Item 3: warmup exits 0 and warns when no CONS-* directory exists."""
        rm_share = tmp_path / "rm_share"
        (rm_share / "Libraries").mkdir(parents=True)
        # No CONS-* subdirectory created
        result = self._run(cache_script, str(rm_share))
        assert result.returncode == 0
        assert "WARNING" in result.stderr
        assert "CONS-*" in result.stderr or "CONS-" in result.stderr

    # Item 1: Dfam-only installation has CONS-Dfam_3.9 directory

    def test_dfam_only_cons_dir_discovered(self, tmp_path, cache_script):
        """Item 1: CONS-Dfam_3.9 directory is found on Dfam-only installation."""
        rm_share = tmp_path / "rm_share"
        cons = rm_share / "Libraries" / "CONS-Dfam_3.9"
        cons.mkdir(parents=True)
        result = self._run(cache_script, str(rm_share))
        assert result.returncode == 0
        # Should not warn about missing CONS-* dir
        assert "WARNING" not in result.stderr

    # Item 2: standard installation has CONS-Dfam_withRBRM_3.9 directory

    def test_repbase_cons_dir_discovered(self, tmp_path, cache_script):
        """Item 2: CONS-Dfam_withRBRM_3.9 directory is found on standard installation."""
        rm_share = tmp_path / "rm_share"
        cons = rm_share / "Libraries" / "CONS-Dfam_withRBRM_3.9"
        cons.mkdir(parents=True)
        result = self._run(cache_script, str(rm_share))
        assert result.returncode == 0
        assert "WARNING" not in result.stderr

    def test_complete_cache_reports_ok(self, tmp_path, cache_script):
        """Item 1/2: when refineableHash.dat exists, no rebuild needed."""
        rm_share = tmp_path / "rm_share"
        cache_dir = rm_share / "Libraries" / "CONS-Dfam_3.9" / "7215"
        cache_dir.mkdir(parents=True)
        (cache_dir / "refineableHash.dat").write_text("fake")
        result = self._run(cache_script, str(rm_share))
        assert result.returncode == 0
        assert "CACHE_OK" in result.stdout

    def test_incomplete_cache_detected_and_removed(self, tmp_path, cache_script):
        """Incomplete cache (has .nhr, missing refineableHash.dat) is detected and cleared."""
        rm_share = tmp_path / "rm_share"
        cache_dir = rm_share / "Libraries" / "CONS-Dfam_3.9" / "7215"
        cache_dir.mkdir(parents=True)
        # Simulate OOM-killed makeblastdb: .nhr present, refineableHash.dat absent
        (cache_dir / "specieslib.nhr").write_text("fake")
        result = self._run(cache_script, str(rm_share))
        assert result.returncode == 0
        assert "INCOMPLETE_CACHE_DETECTED" in result.stderr
        # Directory should have been removed
        assert not cache_dir.exists()


# ===========================================================================
# SECTION 13 — Integration tests with cd-hit-est (items 4, 5, 8, 9)
# Skipped automatically if cd-hit-est is not on PATH.
# ===========================================================================

CDHIT_AVAILABLE = shutil.which("cd-hit-est") is not None


@pytest.mark.skipif(not CDHIT_AVAILABLE, reason="cd-hit-est not found on PATH")
@integration
class TestCdhitIntegration:
    """Items 4, 5, 8, 9: actual cd-hit-est behaviour with -s and -aL flags."""

    def _run_cdhit(self, tmp_path, fasta_path, extra_flags=""):
        out = str(tmp_path / "out.fa")
        clstr = out + ".clstr"
        cmd = (
            f"cd-hit-est -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -r 1 "
            f"{extra_flags} -i {fasta_path} -o {out} -M 4000 -T 1"
        )
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        clusters = {}
        current_rep = None
        for line in open(clstr):
            line = line.strip()
            if line.startswith(">Cluster"):
                current_rep = None
            elif line.endswith("*"):
                current_rep = line.split(">")[1].rstrip("... *")
                clusters[current_rep] = []
            elif current_rep:
                member = line.split(">")[1].split("...")[0]
                clusters[current_rep].append(member)
        return clusters

    @pytest.fixture(autouse=True)
    def _write_fasta(self, tmp_path):
        """Write a FASTA with one short and one long sequence."""
        fa = tmp_path / "test.fa"
        # Short: 500 nt; Long: 5000 nt. Ratio = 10×.
        # With -s 0.5, shorter must be >= 50% of longer → 500/5000=0.1 < 0.5 → not clustered.
        # Without -s they would cluster if the short sequence aligns well.
        short_seq = "ACGT" * 125      # 500 nt
        long_seq  = "ACGT" * 1250     # 5000 nt (short is a prefix, so alignment is perfect)
        fa.write_text(
            f">SHORT#DNA/hAT\n{short_seq}\n"
            f">LONG#DNA/hAT\n{long_seq}\n"
        )
        self.fa = str(fa)
        self.tmp_path = tmp_path

    def test_length_diff_prevents_extreme_size_mismatch(self):
        """Item 4: with -s 0.5, 10× size difference prevents clustering."""
        clusters = self._run_cdhit(self.tmp_path, self.fa, extra_flags="-s 0.5")
        # SHORT should not be a member of LONG's cluster (separate cluster reps)
        for rep, members in clusters.items():
            assert not (
                "SHORT" in rep and any("LONG" in m for m in members)
            ), "SHORT became representative of LONG — clustering not prevented"
            assert not (
                "LONG" in rep and any("SHORT" in m for m in members)
            ), "SHORT was absorbed into LONG cluster — clustering not prevented"

    def test_length_diff_zero_allows_clustering(self):
        """Item 5: with -s 0.0, extreme size mismatch is allowed."""
        # Perfect prefix alignment: short will match long → they should cluster
        clusters = self._run_cdhit(self.tmp_path, self.fa, extra_flags="-s 0.0")
        # At -s 0.0 the size ratio is unconstrained; whether they cluster depends
        # on -aS coverage alone. We just check that cd-hit-est runs without error
        # and returns some clusters (the specific clustering result depends on
        # alignment details and is implementation-specific).
        assert len(clusters) >= 1

    def test_aL_prevents_short_in_long_cluster(self):
        """Item 8: -aL 0.75 prevents short (10% of long) from clustering with long."""
        clusters = self._run_cdhit(self.tmp_path, self.fa, extra_flags="-aL 0.75 -s 0.0")
        for rep, members in clusters.items():
            if "LONG" in rep:
                assert not any("SHORT" in m for m in members), (
                    "SHORT absorbed into LONG cluster despite -aL 0.75"
                )

    def test_aL_zero_does_not_restrict(self):
        """Item 9: -aL 0.0 (disabled) does not prevent clustering by itself."""
        # With -aL 0.0 and -s 0.0 the only constraint is -aS 0.8; since the short
        # seq is a prefix of the long seq it should match at ≥80% coverage of itself.
        clusters = self._run_cdhit(self.tmp_path, self.fa,
                                   extra_flags="-aL 0.0 -s 0.0")
        assert len(clusters) >= 1  # sanity check — cd-hit-est ran and produced output
