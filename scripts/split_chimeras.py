#!/usr/bin/env python3
"""
split_chimeras.py — Post-clustering chimera detection and cluster splitting.

A chimeric cluster representative is one where non-representative members align
to mutually exclusive, non-overlapping windows of the representative. This
indicates the representative spans sequence from ≥2 distinct TE families that
were joined during consensus building.

Algorithm
---------
1. Parse the cd-hit-est .clstr file; for each cluster extract (rep_start,
   rep_end) alignment coordinates for every non-representative member.
2. Build an overlap graph: two members share an edge if their alignment windows
   on the representative overlap by ≥ chimera_overlap_min nucleotides.
3. Find connected components (BFS). ≥2 components → chimeric candidate.
4. Validate: each component must span ≥ chimera_min_component_span fraction of
   the representative's length (guards against trivially partial alignments).
5. For each confirmed chimeric cluster:
   - Retain the original representative, relabelled with a _CHIMERA suffix.
   - For each component, take the longest member as the new representative.
     Its sequence is retrieved from the pre-clustering combined FASTA using
     its original, unmodified header.

Outputs
-------
- <prefix>.chimera_split.fa        Modified FASTA (chimeric reps replaced)
- <prefix>.chimera_detection.tsv   Per-cluster summary
"""

import re
import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, List, Tuple


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class Member:
    index: int
    length: int
    name: str
    is_rep: bool
    seq_start: Optional[int] = None
    seq_end: Optional[int] = None
    rep_start: Optional[int] = None
    rep_end: Optional[int] = None
    strand: Optional[str] = None
    identity: Optional[float] = None


# ── FASTA / .clstr parsers ────────────────────────────────────────────────────

def load_fasta(path: str) -> dict:
    """Return {name: (original_header_line, sequence)} dict.
    `name` is the first word after `>` (used for lookups and de-duplication).
    `original_header_line` is the full text after `>` (preserved on output).
    """
    seqs: dict = {}
    name: Optional[str] = None
    header: Optional[str] = None
    buf: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = (header, ''.join(buf))
                header = line[1:]           # full header text, excluding '>'
                name = header.split()[0]    # first word only (for lookups)
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = (header, ''.join(buf))
    return seqs


def parse_clstr(path: str) -> List[List[Member]]:
    """Return list of clusters; each cluster is a list of Member objects."""
    # With -d 0 cd-hit writes the full name (up to first space) then '...'
    rep_re = re.compile(
        r'^(\d+)\t(\d+)nt,\s+>(\S+)\.\.\.\s+\*'
    )
    mem_re = re.compile(
        r'^(\d+)\t(\d+)nt,\s+>(\S+)\.\.\.'
        r'\s+at\s+(\d+):(\d+):(\d+):(\d+)/([+-])/(\d+\.\d+)%'
    )
    clusters: List[List[Member]] = []
    current: List[Member] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>Cluster'):
                if current:
                    clusters.append(current)
                current = []
                continue
            m = rep_re.match(line)
            if m:
                current.append(Member(
                    index=int(m.group(1)),
                    length=int(m.group(2)),
                    name=m.group(3),
                    is_rep=True,
                ))
                continue
            m = mem_re.match(line)
            if m:
                current.append(Member(
                    index=int(m.group(1)),
                    length=int(m.group(2)),
                    name=m.group(3),
                    is_rep=False,
                    seq_start=int(m.group(4)),
                    seq_end=int(m.group(5)),
                    rep_start=int(m.group(6)),
                    rep_end=int(m.group(7)),
                    strand=m.group(8),
                    identity=float(m.group(9)),
                ))
    if current:
        clusters.append(current)
    return clusters


# ── Overlap graph utilities ───────────────────────────────────────────────────

def _interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _connected_components(members: List[Member], overlap_min: int) -> List[List[int]]:
    """
    Return connected components as lists of indices into `members`.
    Two members are adjacent when their representative alignment windows
    overlap by >= overlap_min nucleotides.
    """
    n = len(members)
    adj: dict = defaultdict(set)
    for i in range(n):
        for j in range(i + 1, n):
            ov = _interval_overlap(
                members[i].rep_start, members[i].rep_end,
                members[j].rep_start, members[j].rep_end,
            )
            if ov >= overlap_min:
                adj[i].add(j)
                adj[j].add(i)

    visited = [False] * n
    components: List[List[int]] = []
    for start in range(n):
        if visited[start]:
            continue
        component: List[int] = []
        queue = [start]
        visited[start] = True
        while queue:
            node = queue.pop()
            component.append(node)
            for neighbour in adj[node]:
                if not visited[neighbour]:
                    visited[neighbour] = True
                    queue.append(neighbour)
        components.append(component)
    return components


# ── Chimera detection ─────────────────────────────────────────────────────────

def detect_chimera(
    cluster: List[Member],
    overlap_min: int,
    min_members: int,
    min_component_span: float,
) -> Tuple[bool, List[List[Member]]]:
    """
    Returns (is_chimeric, validated_components).
    Components are sorted by leftmost alignment position on the representative.
    """
    non_reps = [m for m in cluster if not m.is_rep]
    if len(non_reps) < min_members:
        return False, []

    rep = next((m for m in cluster if m.is_rep), None)
    if rep is None:
        return False, []
    rep_len = rep.length

    raw_components = _connected_components(non_reps, overlap_min)
    if len(raw_components) < 2:
        return False, []

    # Validate: each component must span a meaningful fraction of the rep
    valid = []
    for indices in raw_components:
        comp = [non_reps[i] for i in indices]
        span = max(m.rep_end for m in comp) - min(m.rep_start for m in comp)
        if span / rep_len >= min_component_span:
            valid.append(comp)

    if len(valid) < 2:
        return False, []

    # Sort components by their leftmost position on the representative
    valid.sort(key=lambda c: min(m.rep_start for m in c))
    return True, valid


# ── Name utilities ────────────────────────────────────────────────────────────

def _base_name(name: str) -> str:
    """Strip #classification suffix."""
    return name.split('#')[0]


def _classification(name: str) -> str:
    """Return '#classification' or '' if absent."""
    parts = name.split('#', 1)
    return '#' + parts[1] if len(parts) > 1 else ''


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    # Snakemake injects `snakemake` into globals when using the script: directive
    clstr_path       = snakemake.input.clstr         # noqa: F821
    clustered_fa     = snakemake.input.clustered_fa  # noqa: F821
    combined_fa      = snakemake.input.combined_fa   # noqa: F821
    out_fasta        = snakemake.output.fasta        # noqa: F821
    out_summary      = snakemake.output.summary      # noqa: F821
    log_path         = snakemake.log[0]              # noqa: F821
    overlap_min      = int(snakemake.params.overlap_min)        # noqa: F821
    min_members      = int(snakemake.params.min_members)        # noqa: F821
    min_comp_span    = float(snakemake.params.min_component_span)  # noqa: F821

    logging.basicConfig(
        filename=log_path, level=logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s',
    )
    log = logging.getLogger()

    log.info("Loading clustered representative FASTA (%s)", clustered_fa)
    rep_seqs = load_fasta(clustered_fa)   # {name: (header, seq)}

    log.info("Loading pre-clustering combined FASTA (%s)", combined_fa)
    all_seqs = load_fasta(combined_fa)    # {name: (header, seq)}

    log.info("Parsing .clstr file (%s)", clstr_path)
    clusters = parse_clstr(clstr_path)
    log.info("  %d clusters loaded", len(clusters))

    # Output containers (insertion-ordered dict preserves cluster order)
    out_seqs: dict = {}     # name → (header, sequence)
    summary_rows: list = []

    chimeric_count = 0
    component_reps_added = 0

    for cluster_idx, cluster in enumerate(clusters):
        rep = next((m for m in cluster if m.is_rep), None)
        if rep is None:
            log.warning("Cluster %d has no representative — skipping", cluster_idx)
            continue

        if rep.name not in rep_seqs:
            log.warning(
                "Cluster %d: representative %r not found in clustered FASTA",
                cluster_idx, rep.name,
            )
            continue

        is_chimeric, components = detect_chimera(
            cluster, overlap_min, min_members, min_comp_span,
        )

        if not is_chimeric:
            _header, rep_seq = rep_seqs[rep.name]
            out_seqs[rep.name] = (f">{rep.name}", rep_seq)
            summary_rows.append({
                'cluster_idx':        cluster_idx,
                'representative':     rep.name,
                'rep_length':         rep.length,
                'n_members':          len(cluster) - 1,
                'is_chimeric':        False,
                'n_components':       1,
                'component_sizes':    str(len(cluster) - 1),
                'chimera_score':      0.0,
                'component_rep_names': rep.name,
            })
            continue

        # ── Chimeric cluster ─────────────────────────────────────────────
        chimeric_count += 1
        n_non_reps = len(cluster) - 1
        log.info(
            "Cluster %d: chimeric representative %r (%d nt, %d members, %d components)",
            cluster_idx, rep.name, rep.length, n_non_reps, len(components),
        )

        # Retain chimeric rep with _CHIMERA label for traceability
        _header, rep_seq = rep_seqs[rep.name]
        chimera_name = f"{_base_name(rep.name)}_CHIMERA{_classification(rep.name)}"
        out_seqs[chimera_name] = (f">{chimera_name}", rep_seq)

        # Chimera score: largest inter-component gap / rep length
        gaps: List[int] = []
        for i in range(len(components) - 1):
            left_end   = max(m.rep_end   for m in components[i])
            right_start = min(m.rep_start for m in components[i + 1])
            gaps.append(max(0, right_start - left_end))
        chimera_score = round(max(gaps) / rep.length, 4) if gaps else 0.0

        comp_rep_names: List[str] = []
        for comp_num, comp in enumerate(components, start=1):
            longest = max(comp, key=lambda m: m.length)
            comp_seq = all_seqs.get(longest.name)
            if comp_seq is None:
                log.warning(
                    "  Component %d: longest member %r not found in combined FASTA"
                    " — component skipped",
                    comp_num, longest.name,
                )
                continue

            # Preserve the original sequence name and header from the combined FASTA
            # so that TE classification tags and species prefixes are retained.
            comp_name = longest.name
            if comp_name in out_seqs:
                log.warning(
                    "  Component %d: name %r already present in output (collision) "
                    "— component skipped",
                    comp_num, comp_name,
                )
                continue
            _orig_header, comp_seq = all_seqs[longest.name]
            # Write with the full original header line (including any comment fields)
            out_seqs[comp_name] = (f">{_orig_header}", comp_seq)
            comp_rep_names.append(comp_name)
            component_reps_added += 1

            coord_range = (
                f"{min(m.rep_start for m in comp)}"
                f"-{max(m.rep_end for m in comp)}"
            )
            log.info(
                "  Component %d: %d members, rep coverage %s, "
                "new representative %r (%d nt)",
                comp_num, len(comp), coord_range, comp_name, longest.length,
            )

        summary_rows.append({
            'cluster_idx':        cluster_idx,
            'representative':     rep.name,
            'rep_length':         rep.length,
            'n_members':          n_non_reps,
            'is_chimeric':        True,
            'n_components':       len(components),
            'component_sizes':    ';'.join(str(len(c)) for c in components),
            'chimera_score':      chimera_score,
            'component_rep_names': ';'.join(comp_rep_names),
        })

    # ── Write output FASTA ────────────────────────────────────────────────────
    with open(out_fasta, 'w') as fh:
        for header, seq in out_seqs.values():
            fh.write(header + '\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + '\n')

    # ── Write summary TSV ─────────────────────────────────────────────────────
    cols = [
        'cluster_idx', 'representative', 'rep_length', 'n_members',
        'is_chimeric', 'n_components', 'component_sizes',
        'chimera_score', 'component_rep_names',
    ]
    with open(out_summary, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for row in summary_rows:
            fh.write('\t'.join(str(row[c]) for c in cols) + '\n')

    log.info(
        "Done. %d chimeric clusters detected; %d component representatives written.",
        chimeric_count, component_reps_added,
    )
    log.info("Output FASTA  : %s", out_fasta)
    log.info("Summary TSV   : %s", out_summary)


# Guard: only execute when run as a Snakemake script (snakemake object is
# injected into globals by Snakemake before executing the script file).
# When imported directly (e.g. for testing), main() is NOT called.
if 'snakemake' in dir():
    main()
