"""shared_unique_plot.py — shared/unique TE content plots for EarlGreyParTEA.

Generates stacked bar charts showing the number of TE families and the genome
coverage (bp and %) that are shared across species vs unique to a single species.

Bars are coloured by TE classification using the EarlGrey palette.  Shared
families/coverage are drawn in the full class colour; unique families/coverage
are drawn in a lightened version of the same colour.  This makes it easy to
distinguish both the shared/unique split AND the TE class composition at a
glance.

Two detection modes
-------------------
cluster (full pipeline mode)
    Parses the cd-hit-est .clstr file to determine which species contributed
    sequences to each cluster.  A cluster is 'shared' if it contains sequences
    from ≥2 species.  TE class is read from the representative sequence header
    (which carries a '#class/subclass' suffix) and confirmed/updated from each
    species' GFF column 3.

presence_absence (annotate pipeline mode)
    No .clstr file is available.  The NAME= attribute and TE class (column 3)
    from each species' filteredRepeats.gff are used.  A family is 'shared' if
    the same NAME= value appears in ≥2 species' GFF files.  Sequence-level
    divergence cannot be accounted for; homologous families annotated under
    different names will appear as unique in each species.

GFF format produced by EarlGrey
--------------------------------
Column 3 (0-based col[2]) — raw TE type, e.g. "LINE/R2-Hero", "LTR/Gypsy".
Column 9 — attributes including NAME=<family_name>.  Example:
    ID=GENOME1_RND-1_FAMILY-34_3;NAME=GENOME1_RND-1_FAMILY-34;TSTART=1;…
GFF coordinates are 1-based, end-inclusive → bp = end - start + 1.

Outputs
-------
shared_unique_families.{pdf,tsv}
    Family-count stacked bars (by TE class, shared vs unique) + data table.
shared_unique_coverage.{pdf,tsv}
    Genome-coverage stacked bars + data table.
shared_unique_families_phylo.pdf   (only when tree_path is set)
shared_unique_coverage_phylo.pdf
    Same plots with species in phylogenetic order and a cladogram sidebar.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_utils import parse_clstr as _parse_clstr_full

# ---------------------------------------------------------------------------
# EarlGrey TE classification colour palette
# Matches the autoPie / R colour definitions used throughout EarlGrey.
# ---------------------------------------------------------------------------
CLASS_PALETTE = {
    "LTR":            "#00782A",
    "LINE":           "#0098D4",
    "SINE":           "#9B0056",
    "DNA":            "#E32017",
    "Rolling Circle": "#EE7C0E",
    "Penelope":       "#7156A5",
    "Other":          "#F3A9BB",
    "Unclassified":   "#A0A5A9",
}

# Stacking order for bar segments — loosely by typical TE abundance.
CLASS_ORDER = [
    "LTR", "LINE", "SINE", "DNA",
    "Rolling Circle", "Penelope", "Other", "Unclassified",
]

# Human-readable label for "Other" used in legend and TSV headers.
OTHER_LABEL = "Other (Simple repeat / Satellite / Low-complexity / RNA)"


def _lighten(hex_colour, factor=0.55):
    """Return a lightened hex colour by blending with white at *factor* (0–1)."""
    r = int(hex_colour[1:3], 16)
    g = int(hex_colour[3:5], 16)
    b = int(hex_colour[5:7], 16)
    r = round(r + (255 - r) * factor)
    g = round(g + (255 - g) * factor)
    b = round(b + (255 - b) * factor)
    return f"#{r:02X}{g:02X}{b:02X}"


# Pre-computed lightened palette for unique segments.
LIGHT_PALETTE = {cls: _lighten(col) for cls, col in CLASS_PALETTE.items()}


def _classify_te(te_type):
    """Map a raw GFF column-3 TE type to the top-level EarlGrey category.

    Checks are ordered so that subtypes that share a common prefix (e.g.
    LINE/Penelope vs LINE) are resolved correctly.
    """
    t = te_type.strip()
    if t.startswith("LTR"):
        return "LTR"
    # Penelope before generic LINE to avoid mis-classifying LINE/Penelope.
    if t.startswith("LINE/Penelope") or t.lower().startswith("penelope"):
        return "Penelope"
    if t.startswith("LINE"):
        return "LINE"
    if t.startswith("SINE"):
        return "SINE"
    if t.startswith("DNA"):
        return "DNA"
    if (t.startswith("RC") or "helitron" in t.lower()
            or t.lower().startswith("rolling")):
        return "Rolling Circle"
    if any(t.startswith(s) for s in (
        "Simple_repeat", "Low_complexity", "Satellite",
        "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA",
        "Microsatellite", "ARTEFACT",
    )):
        return "Other"
    return "Unclassified"


# ---------------------------------------------------------------------------
# GFF / FASTA parsing helpers
# ---------------------------------------------------------------------------

def _parse_gff_name(attr_string):
    """Extract the NAME= value from a GFF attribute string.

    Handles semi-colon-separated fields.  Returns None if NAME= is absent.
    """
    for part in attr_string.split(";"):
        part = part.strip()
        if part.upper().startswith("NAME="):
            return part[5:].strip()
    return None


def _gff_coverage_and_families(gff_path):
    """Parse a filteredRepeats.gff and return class-annotated family data.

    Returns
    -------
    family_class : dict{ family_name -> top_class }
        Mapping of each unique NAME= value to its (last-seen) EarlGrey class.
    hits : list of (name, top_class, bp, is_nested)
        Per-annotation-hit tuples.  is_nested is True when the GFF attribute
        field contains ``NESTED=FULLY_NESTED``, matching the EarlGrey convention
        used in highLevelCount.txt (e.g. "DNA-nested").
    """
    family_class = {}
    hits = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue
            name = _parse_gff_name(parts[8])
            if name is None:
                continue
            top_class = _classify_te(parts[2])
            bp = end - start + 1
            is_nested = "NESTED=FULLY_NESTED" in parts[8].upper()
            family_class[name] = top_class
            hits.append((name, top_class, bp, is_nested))
    return family_class, hits


def _genome_size_from_prep(prep_path):
    """Count total non-header bp in a FASTA file."""
    total = 0
    with open(prep_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def _empty_class_dict():
    """Return a zeroed dict keyed by CLASS_ORDER."""
    return {cls: 0 for cls in CLASS_ORDER}


# ---------------------------------------------------------------------------
# Mode A: cluster-based
# ---------------------------------------------------------------------------

def _cluster_mode(species_list, gff_paths, prep_paths, clstr_file):
    """Return (fam_data, cov_data) using cluster membership from .clstr file.

    fam_data : dict[ species -> { 'shared_by_class': {cls: int},
                                   'unique_by_class': {cls: int} } ]
    cov_data : dict[ species -> { 'shared_bp_by_class': {cls: int},
                                   'unique_bp_by_class': {cls: int},
                                   'genome_size': int } ]
    """
    _, existing_clusters, cluster_to_species, rep_name_to_cluster = \
        _parse_clstr_full(clstr_file, species_list)

    # Derive cluster → TE class from representative sequence names.
    # Keys in rep_name_to_cluster that include '#' carry the class suffix.
    cluster_class = {}
    for stored_name, cid in rep_name_to_cluster.items():
        if "#" in stored_name:
            raw_class = stored_name.split("#", 1)[1]
            cluster_class[cid] = _classify_te(raw_class)

    # Coverage: parse each species' GFF, match hits to clusters, accumulate.
    # GFF data is also used to confirm/update cluster_class (GFF is authoritative).
    # NOTE: annotate.smk applies `toupper($9)` to the GFF attributes, so NAME=
    # values are always uppercase.  rep_name_to_cluster keys are lowercase (from
    # the cluster file), so we build a case-folded alias map for lookup.
    rep_name_lower = {k.lower(): v for k, v in rep_name_to_cluster.items()}

    cov_data = {}
    for sp, gff_path, prep_path in zip(species_list, gff_paths, prep_paths):
        genome_size = _genome_size_from_prep(prep_path)
        _, hits = _gff_coverage_and_families(gff_path)
        shared_bp        = _empty_class_dict()
        unique_bp        = _empty_class_dict()
        shared_nested_bp = _empty_class_dict()
        unique_nested_bp = _empty_class_dict()
        for name, top_class, bp, is_nested in hits:
            cid = rep_name_to_cluster.get(name)
            if cid is None:
                cid = rep_name_lower.get(name.lower())
            if cid is None:
                bare = name.split("#")[0]
                cid = rep_name_to_cluster.get(bare)
            if cid is None:
                cid = rep_name_lower.get(bare.lower())
            if cid is not None:
                cluster_class[cid] = top_class   # GFF overrides header-derived class
                sp_set = cluster_to_species.get(cid, frozenset())
                if is_nested:
                    if len(sp_set) == 1:
                        unique_nested_bp[top_class] += bp
                    else:
                        shared_nested_bp[top_class] += bp
                else:
                    if len(sp_set) == 1:
                        unique_bp[top_class] += bp
                    else:
                        shared_bp[top_class] += bp
            # Hits with no matching cluster are rare (simple repeats not in
            # the combined library) and are silently omitted.
        cov_data[sp] = {
            "shared_bp_by_class":        shared_bp,
            "unique_bp_by_class":         unique_bp,
            "shared_nested_bp_by_class": shared_nested_bp,
            "unique_nested_bp_by_class": unique_nested_bp,
            "genome_size": genome_size,
        }

    # Family counts: one count per cluster, attributed to each contributing species.
    fam_data = {sp: {
        "shared_by_class": _empty_class_dict(),
        "unique_by_class": _empty_class_dict(),
    } for sp in species_list}
    for cid, sp_set in cluster_to_species.items():
        top_class = cluster_class.get(cid, "Unclassified")
        if len(sp_set) == 1:
            sp = next(iter(sp_set))
            if sp in fam_data:
                fam_data[sp]["unique_by_class"][top_class] += 1
        else:
            for sp in sp_set:
                if sp in fam_data:
                    fam_data[sp]["shared_by_class"][top_class] += 1

    return fam_data, cov_data


# ---------------------------------------------------------------------------
# Mode B: presence/absence-based
# ---------------------------------------------------------------------------

def _presence_absence_mode(species_list, gff_paths, prep_paths):
    """Return (fam_data, cov_data) using cross-species presence/absence.

    fam_data / cov_data have the same structure as _cluster_mode output.
    """
    family_species = {}   # name -> set[species]
    family_class   = {}   # name -> top_class (last-seen)
    family_hits    = {}   # name -> {species -> total_bp}
    genome_sizes   = {}

    family_nested_hits = {}  # name -> {species -> total_nested_bp}

    for sp, gff_path, prep_path in zip(species_list, gff_paths, prep_paths):
        genome_sizes[sp] = _genome_size_from_prep(prep_path)
        fc, hits = _gff_coverage_and_families(gff_path)
        for name, cls in fc.items():
            family_species.setdefault(name, set()).add(sp)
            family_class[name] = cls
        for name, top_class, bp, is_nested in hits:
            if is_nested:
                family_nested_hits.setdefault(name, {}).setdefault(sp, 0)
                family_nested_hits[name][sp] += bp
            else:
                family_hits.setdefault(name, {}).setdefault(sp, 0)
                family_hits[name][sp] += bp

    # Family counts (nested hits don't change family membership)
    fam_data = {sp: {
        "shared_by_class": _empty_class_dict(),
        "unique_by_class": _empty_class_dict(),
    } for sp in species_list}
    for name, sp_set in family_species.items():
        cls = family_class.get(name, "Unclassified")
        for sp in sp_set:
            if sp not in fam_data:
                continue
            if len(sp_set) == 1:
                fam_data[sp]["unique_by_class"][cls] += 1
            else:
                fam_data[sp]["shared_by_class"][cls] += 1

    # Coverage — non-nested bp (used for plots)
    cov_data = {sp: {
        "shared_bp_by_class":        _empty_class_dict(),
        "unique_bp_by_class":         _empty_class_dict(),
        "shared_nested_bp_by_class": _empty_class_dict(),
        "unique_nested_bp_by_class": _empty_class_dict(),
        "genome_size": genome_sizes[sp],
    } for sp in species_list}
    for name, sp_bp_map in family_hits.items():
        sp_set = family_species.get(name, set())
        cls = family_class.get(name, "Unclassified")
        is_shared = len(sp_set) > 1
        for sp, bp in sp_bp_map.items():
            if sp not in cov_data:
                continue
            if is_shared:
                cov_data[sp]["shared_bp_by_class"][cls] += bp
            else:
                cov_data[sp]["unique_bp_by_class"][cls] += bp
    # Nested bp — recorded separately in TSV, excluded from plots
    for name, sp_bp_map in family_nested_hits.items():
        sp_set = family_species.get(name, set())
        cls = family_class.get(name, "Unclassified")
        is_shared = len(sp_set) > 1
        for sp, bp in sp_bp_map.items():
            if sp not in cov_data:
                continue
            if is_shared:
                cov_data[sp]["shared_nested_bp_by_class"][cls] += bp
            else:
                cov_data[sp]["unique_nested_bp_by_class"][cls] += bp

    return fam_data, cov_data


# ---------------------------------------------------------------------------
# TSV writers
# ---------------------------------------------------------------------------

def _write_family_tsv(path, species_order, fam_data, method):
    col_s = [f"shared_{cls.replace(' ', '_')}_families" for cls in CLASS_ORDER]
    col_u = [f"unique_{cls.replace(' ', '_')}_families" for cls in CLASS_ORDER]
    with open(path, "w") as fh:
        if method == "presence_absence":
            fh.write(
                "# NOTE: presence/absence classification — sequence-level "
                "divergence not accounted for. Homologous families with "
                "different names will appear as unique in each species.\n"
            )
        fh.write(
            "species\tshared_families\tunique_families\ttotal_families\t"
            + "\t".join(col_s) + "\t" + "\t".join(col_u) + "\tmethod\n"
        )
        for sp in species_order:
            d = fam_data[sp]
            shared = sum(d["shared_by_class"].values())
            unique = sum(d["unique_by_class"].values())
            sc = "\t".join(str(d["shared_by_class"].get(cls, 0)) for cls in CLASS_ORDER)
            uc = "\t".join(str(d["unique_by_class"].get(cls, 0)) for cls in CLASS_ORDER)
            fh.write(
                f"{sp}\t{shared}\t{unique}\t{shared + unique}\t{sc}\t{uc}\t{method}\n"
            )


def _write_coverage_tsv(path, species_order, cov_data, method):
    # Column name helpers — follow highLevelCount.txt convention: "{class}-nested"
    def _col(prefix, cls, nested=False):
        safe = cls.replace(' ', '_')
        return f"{prefix}_{safe}-nested_bp" if nested else f"{prefix}_{safe}_bp"

    col_s  = [_col("shared", cls)          for cls in CLASS_ORDER]
    col_u  = [_col("unique", cls)          for cls in CLASS_ORDER]
    col_sn = [_col("shared", cls, True)    for cls in CLASS_ORDER]
    col_un = [_col("unique", cls, True)    for cls in CLASS_ORDER]
    with open(path, "w") as fh:
        if method == "presence_absence":
            fh.write(
                "# NOTE: presence/absence classification — sequence-level "
                "divergence not accounted for.\n"
            )
        fh.write(
            "species\tshared_bp\tunique_bp\ttotal_te_bp\tgenome_size_bp\t"
            "shared_pct\tunique_pct\t"
            + "\t".join(col_s) + "\t" + "\t".join(col_u) + "\t"
            + "\t".join(col_sn) + "\t" + "\t".join(col_un) + "\tmethod\n"
        )
        for sp in species_order:
            d = cov_data[sp]
            shared_bp = sum(d["shared_bp_by_class"].values())
            unique_bp = sum(d["unique_bp_by_class"].values())
            gs = d["genome_size"]
            spct = 100.0 * shared_bp / gs if gs > 0 else 0.0
            upct = 100.0 * unique_bp / gs if gs > 0 else 0.0
            sc  = "\t".join(str(d["shared_bp_by_class"].get(cls, 0))        for cls in CLASS_ORDER)
            uc  = "\t".join(str(d["unique_bp_by_class"].get(cls, 0))         for cls in CLASS_ORDER)
            scn = "\t".join(str(d["shared_nested_bp_by_class"].get(cls, 0)) for cls in CLASS_ORDER)
            ucn = "\t".join(str(d["unique_nested_bp_by_class"].get(cls, 0)) for cls in CLASS_ORDER)
            fh.write(
                f"{sp}\t{shared_bp}\t{unique_bp}\t{shared_bp + unique_bp}\t{gs}\t"
                f"{spct:.4f}\t{upct:.4f}\t{sc}\t{uc}\t{scn}\t{ucn}\t{method}\n"
            )


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _legend_label(cls):
    return OTHER_LABEL if cls == "Other" else cls


def _stacked_bar(ax, species_order, shared_by_class, unique_by_class,
                 xlabel, title, value_fmt="{:.0f}"):
    """Draw a TE-class-coloured stacked horizontal bar chart onto *ax*.

    Each bar has two halves:
      • Left  (solid class colour, no hatch)   — shared (≥2 species).
      • Right (lightened colour + o hatch)  — unique (1 species).
    A thin vertical black line is drawn at the shared→unique boundary.
    Within each half, segments are stacked in CLASS_ORDER.

    Parameters
    ----------
    shared_by_class, unique_by_class : dict{ class -> ndarray of shape (n,) }
        Per-class values indexed to match *species_order*.
    """
    n = len(species_order)
    y = np.arange(n)
    bar_h = 0.55

    classes_present = set()

    # Shared segments — solid class colours, no hatch
    left = np.zeros(n)
    for cls in CLASS_ORDER:
        vals = shared_by_class.get(cls, np.zeros(n))
        if np.any(vals > 0):
            ax.barh(y, vals, left=left, height=bar_h,
                    color=CLASS_PALETTE[cls], edgecolor="none")
            classes_present.add(cls)
        left = left + vals

    shared_totals = left.copy()

    # Precompute unique totals to know where to draw dividing lines
    unique_totals = np.zeros(n)
    for cls in CLASS_ORDER:
        unique_totals += unique_by_class.get(cls, np.zeros(n))

    # Dividing line between shared and unique portions of each bar
    for i in range(n):
        if shared_totals[i] > 0 and unique_totals[i] > 0:
            ax.plot(
                [shared_totals[i], shared_totals[i]],
                [i - bar_h / 2, i + bar_h / 2],
                color="#eeeeee", lw=1.5, zorder=5, solid_capstyle="butt",
            )

    # Unique segments — lightened colour + diagonal hatch for clear contrast
    left_u = left.copy()
    for cls in CLASS_ORDER:
        vals = unique_by_class.get(cls, np.zeros(n))
        if np.any(vals > 0):
            ax.barh(y, vals, left=left_u, height=bar_h,
                    color=LIGHT_PALETTE[cls], hatch="o",
                    edgecolor="#eeeeee", linewidth=0)
            classes_present.add(cls)
        left_u = left_u + vals

    # Value labels at the right end of each bar (total = shared + unique)
    totals = left_u
    for i, total in enumerate(totals):
        if total > 0:
            ax.text(
                total + total * 0.01, i, value_fmt.format(total),
                va="center", ha="left", fontsize=8, color="#444444",
            )

    ax.set_yticks(y)
    ax.set_yticklabels(species_order, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    xlim = ax.get_xlim()
    ax.set_xlim(0, xlim[1] * 1.15)

    # Legend outside axes (right side) so it never overlaps bars
    shade_patches = [
        mpatches.Patch(facecolor="#555555", edgecolor="none",
                       label="Shared (≥2 species)"),
        mpatches.Patch(facecolor="#BBBBBB", hatch="o", edgecolor="#555555",
                       linewidth=0, label="Unique (1 species)"),
    ]
    class_patches = [
        mpatches.Patch(facecolor=CLASS_PALETTE[cls], edgecolor="none",
                       label=_legend_label(cls))
        for cls in CLASS_ORDER
        if cls in classes_present
    ]
    ax.legend(
        handles=shade_patches + class_patches,
        frameon=False, fontsize=7.5,
        loc="upper left", bbox_to_anchor=(1.01, 1),
    )


# ---------------------------------------------------------------------------
# Cladogram helper
# ---------------------------------------------------------------------------

def _draw_cladogram(ax, species_order, tree_path):
    """Draw a unit-branch cladogram on *ax* with tips matching *species_order*.

    Uses Bio.Phylo to read the newick tree; draws horizontal and vertical
    connecting lines manually so no plt.show() or new figure is created.
    The y-axis of *ax* must span 0..len(species_order)-1.
    """
    try:
        from Bio import Phylo
    except ImportError:
        ax.axis("off")
        ax.text(0.5, 0.5, "biopython\nnot installed",
                ha="center", va="center", transform=ax.transAxes, fontsize=7)
        return

    with open(tree_path) as fh:
        tree = Phylo.read(fh, "newick")
    tree.ladderize(reverse=True)

    name_to_y = {sp: i for i, sp in enumerate(species_order)}
    if not any(t.name in name_to_y for t in tree.get_terminals()):
        ax.axis("off")
        return

    def _assign_depth(clade, depth=0):
        clade._depth = depth
        for child in clade.clades:
            _assign_depth(child, depth + 1)

    _assign_depth(tree.root)
    max_depth = max(t._depth for t in tree.get_terminals())

    def _draw_clade(clade):
        if clade.is_terminal():
            y = name_to_y.get(clade.name)
            if y is None:
                return
            ax.plot([clade._depth, max_depth], [y, y],
                    color="#444444", lw=1.2, solid_capstyle="round")
        else:
            child_ys = []
            for child in clade.clades:
                _draw_clade(child)
                if child.is_terminal():
                    cy = name_to_y.get(child.name)
                else:
                    c_tips = [t for t in child.get_terminals() if t.name in name_to_y]
                    cy = np.mean([name_to_y[t.name] for t in c_tips]) if c_tips else None
                if cy is not None:
                    child_ys.append(cy)
                    ax.plot([clade._depth, child._depth], [cy, cy],
                            color="#444444", lw=1.2, solid_capstyle="round")
            if len(child_ys) >= 2:
                ax.plot([clade._depth, clade._depth],
                        [min(child_ys), max(child_ys)],
                        color="#444444", lw=1.2, solid_capstyle="round")
                # Bootstrap value: FastTree outputs SH-like local support (0–1).
                # Display as a percentage integer at the internal node.
                # Skip root (no confidence) and trivial 100% nodes to reduce clutter.
                conf = getattr(clade, "confidence", None)
                if conf is not None and clade is not tree.root:
                    pct = int(round(conf * 100))
                    node_y = np.mean(child_ys)
                    ax.text(clade._depth + 0.05, node_y,
                            str(pct),
                            ha="left", va="center",
                            fontsize=6, color="#666666")

    _draw_clade(tree.root)
    ax.set_xlim(-0.5, max_depth + 0.5)
    ax.set_ylim(-0.5, len(species_order) - 0.5)
    ax.axis("off")


def _make_plot(species_order, shared_by_class, unique_by_class,
               xlabel, title, out_path,
               value_fmt="{:.0f}", tree_path=None):
    """Save a class-coloured stacked-bar plot (with optional cladogram sidebar)."""
    fig_h = max(4, len(species_order) * 0.65 + 1.8)
    if tree_path:
        fig, (ax_tree, ax_bar) = plt.subplots(
            1, 2, figsize=(12, fig_h),
            gridspec_kw={"width_ratios": [1, 3]},
        )
        _draw_cladogram(ax_tree, species_order, tree_path)
        ax_bar.set_yticks(range(len(species_order)))
        ax_bar.set_yticklabels(species_order, fontsize=9)
        _stacked_bar(ax_bar, species_order, shared_by_class, unique_by_class,
                     xlabel, title, value_fmt)
    else:
        fig, ax_bar = plt.subplots(figsize=(10, fig_h))
        _stacked_bar(ax_bar, species_order, shared_by_class, unique_by_class,
                     xlabel, title, value_fmt)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Phylo-ordered species list helper
# ---------------------------------------------------------------------------

def _phylo_order(species_list, tree_path):
    """Return species in phylogenetic tip order (ladderized).

    Species absent from the tree are appended at the end in their original order.
    """
    try:
        from Bio import Phylo
    except ImportError:
        return list(species_list)
    with open(tree_path) as fh:
        tree = Phylo.read(fh, "newick")
    tree.ladderize(reverse=True)
    tip_names = [t.name for t in tree.get_terminals()]
    ordered   = [t for t in tip_names if t in set(species_list)]
    remaining = [s for s in species_list if s not in set(ordered)]
    return ordered + remaining


# ---------------------------------------------------------------------------
# Array-building helpers
# ---------------------------------------------------------------------------

def _class_arrays(species_order, data_dict, key_prefix):
    """Build per-class numpy arrays from fam_data or cov_data.

    *key_prefix* should be 'shared' or 'unique' for fam_data (resolves to
    'shared_by_class') and 'shared_bp' or 'unique_bp' for cov_data (resolves
    to 'shared_bp_by_class').
    """
    key = key_prefix + "_by_class"
    return {
        cls: np.array(
            [data_dict[sp][key].get(cls, 0) for sp in species_order],
            dtype=float,
        )
        for cls in CLASS_ORDER
    }


# ---------------------------------------------------------------------------
# Entry point (called as a Snakemake script: directive)
# ---------------------------------------------------------------------------

def main():
    species_list   = list(snakemake.params.species)                       # noqa: F821
    detection_mode = snakemake.params.detection_mode                      # noqa: F821
    has_phylo_tree = bool(snakemake.params.get("has_phylo_tree", False))  # noqa: F821

    gff_paths  = list(snakemake.input.gffs)   # noqa: F821
    prep_paths = list(snakemake.input.preps)  # noqa: F821

    fam_pdf = snakemake.output.fam_pdf  # noqa: F821
    fam_tsv = snakemake.output.fam_tsv  # noqa: F821
    cov_pdf = snakemake.output.cov_pdf  # noqa: F821
    cov_tsv = snakemake.output.cov_tsv  # noqa: F821

    fam_phylo_pdf = getattr(snakemake.output, "fam_phylo_pdf", None)  # noqa: F821
    cov_phylo_pdf = getattr(snakemake.output, "cov_phylo_pdf", None)  # noqa: F821
    tree_path     = getattr(snakemake.input,  "tree",         None)   # noqa: F821

    os.makedirs(os.path.dirname(os.path.abspath(fam_pdf)), exist_ok=True)

    # ---- Compute per-class shared/unique data ----
    if detection_mode == "cluster":
        clstr_file = snakemake.input.clstr   # noqa: F821
        fam_data, cov_data = _cluster_mode(
            species_list, gff_paths, prep_paths, clstr_file
        )
        method_label = "cluster"
        print("[shared_unique] Using cluster-based detection (full mode)", flush=True)
    else:
        fam_data, cov_data = _presence_absence_mode(
            species_list, gff_paths, prep_paths
        )
        method_label = "presence_absence"
        print(
            "[shared_unique] Using presence/absence detection (annotate mode). "
            "Homologous families under different names appear as unique.",
            flush=True,
        )

    species_order = sorted(species_list)
    genome_sizes  = np.array(
        [cov_data[sp]["genome_size"] for sp in species_order], dtype=float
    )

    # Per-class arrays for alphabetical plots
    shared_fam_cls = _class_arrays(species_order, fam_data, "shared")
    unique_fam_cls = _class_arrays(species_order, fam_data, "unique")
    shared_cov_cls = _class_arrays(species_order, cov_data, "shared_bp")
    unique_cov_cls = _class_arrays(species_order, cov_data, "unique_bp")

    # Convert bp coverage → % genome
    def _to_pct(d):
        return {
            cls: np.where(genome_sizes > 0, 100.0 * vals / genome_sizes, 0.0)
            for cls, vals in d.items()
        }

    shared_cov_pct = _to_pct(shared_cov_cls)
    unique_cov_pct = _to_pct(unique_cov_cls)

    # ---- Write TSVs ----
    _write_family_tsv(fam_tsv, species_order, fam_data, method_label)
    _write_coverage_tsv(cov_tsv, species_order, cov_data, method_label)
    print(f"[shared_unique] Saved family TSV   → {fam_tsv}", flush=True)
    print(f"[shared_unique] Saved coverage TSV → {cov_tsv}", flush=True)

    # ---- Standard (alphabetical) plots ----
    _make_plot(
        species_order, shared_fam_cls, unique_fam_cls,
        "TE families", "Shared vs Unique TE Families per Species", fam_pdf,
    )
    _make_plot(
        species_order, shared_cov_pct, unique_cov_pct,
        "% genome", "Shared vs Unique TE Coverage per Species",
        cov_pdf, value_fmt="{:.2f}%",
    )
    print(f"[shared_unique] Saved family plot   → {fam_pdf}", flush=True)
    print(f"[shared_unique] Saved coverage plot → {cov_pdf}", flush=True)

    # ---- Phylo-ordered plots ----
    if has_phylo_tree and tree_path and fam_phylo_pdf and cov_phylo_pdf:
        phylo_order = _phylo_order(species_list, tree_path)
        gs_p = np.array(
            [cov_data[sp]["genome_size"] for sp in phylo_order], dtype=float
        )

        sf_p = _class_arrays(phylo_order, fam_data, "shared")
        uf_p = _class_arrays(phylo_order, fam_data, "unique")
        sc_p_bp = _class_arrays(phylo_order, cov_data, "shared_bp")
        uc_p_bp = _class_arrays(phylo_order, cov_data, "unique_bp")
        sc_p = {
            cls: np.where(gs_p > 0, 100.0 * vals / gs_p, 0.0)
            for cls, vals in sc_p_bp.items()
        }
        uc_p = {
            cls: np.where(gs_p > 0, 100.0 * vals / gs_p, 0.0)
            for cls, vals in uc_p_bp.items()
        }

        _make_plot(
            phylo_order, sf_p, uf_p,
            "TE families",
            "Shared vs Unique TE Families per Species (Phylogenetic Order)",
            fam_phylo_pdf, tree_path=tree_path,
        )
        _make_plot(
            phylo_order, sc_p, uc_p,
            "% genome",
            "Shared vs Unique TE Coverage per Species (Phylogenetic Order)",
            cov_phylo_pdf, value_fmt="{:.2f}%", tree_path=tree_path,
        )
        print(f"[shared_unique] Saved phylo family plot   → {fam_phylo_pdf}", flush=True)
        print(f"[shared_unique] Saved phylo coverage plot → {cov_phylo_pdf}", flush=True)


# Only execute when run as a Snakemake script (snakemake object is injected).
if "snakemake" in dir():
    main()
