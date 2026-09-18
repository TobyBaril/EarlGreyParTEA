"""busco_summary_plot.py — BUSCO completeness summary plots for EarlGreyParTEA.

Reads BUSCO short_summary.txt files (one per species) and generates a
horizontal stacked bar chart showing the completeness classification
for each genome, using the standard BUSCO colour scheme.

Called as a Snakemake script (script: directive).

snakemake.input.summaries : list of paths to BUSCO short_summary.txt files
snakemake.params.species  : list of species names (same order as summaries)
snakemake.output.pdf      : output PDF path
snakemake.output.tsv      : output TSV path
"""

import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Standard BUSCO colour scheme
BUSCO_COLOURS = {
    "Complete (single-copy)": "#1E90FF",   # dodger blue
    "Complete (duplicated)":  "#FF7F00",   # orange
    "Fragmented":             "#FFDD00",   # yellow
    "Missing":                "#D83232",   # red
}
BUSCO_ORDER = [
    "Complete (single-copy)",
    "Complete (duplicated)",
    "Fragmented",
    "Missing",
]
BUSCO_ABBREV = {
    "Complete (single-copy)": "S",
    "Complete (duplicated)":  "D",
    "Fragmented":             "F",
    "Missing":                "M",
}


def _parse_busco_summary(path):
    """Parse a BUSCO short_summary.txt and return a dict of raw counts.

    Handles both BUSCO v4/v5 short_summary format.
    Returns dict with keys: single, duplicated, fragmented, missing, total.
    """
    counts = {"single": 0, "duplicated": 0, "fragmented": 0, "missing": 0, "total": 0}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            # v5 format: "C:xx%[S:xx%,D:xx%],F:xx%,M:xx%,n:NN"
            m = re.search(r"C:[\d.]+%\[S:[\d.]+%,D:[\d.]+%\],F:[\d.]+%,M:[\d.]+%,n:(\d+)", line)
            if m:
                counts["total"] = int(m.group(1))
            # Numeric counts lines
            m = re.match(r"\s*(\d+)\s+Complete\s+and\s+single.copy", line, re.I)
            if m:
                counts["single"] = int(m.group(1))
            m = re.match(r"\s*(\d+)\s+Complete\s+and\s+duplicated", line, re.I)
            if m:
                counts["duplicated"] = int(m.group(1))
            m = re.match(r"\s*(\d+)\s+Fragmented", line, re.I)
            if m:
                counts["fragmented"] = int(m.group(1))
            m = re.match(r"\s*(\d+)\s+Missing", line, re.I)
            if m:
                counts["missing"] = int(m.group(1))
            m = re.match(r"\s*(\d+)\s+Total BUSCO groups", line, re.I)
            if m:
                counts["total"] = int(m.group(1))
    return counts


def _write_tsv(path, species_order, all_counts):
    with open(path, "w") as fh:
        fh.write(
            "species\tcomplete_single\tcomplete_duplicated\tfragmented\t"
            "missing\ttotal\tcomplete_pct\tfragmented_pct\tmissing_pct\n"
        )
        for sp in species_order:
            c = all_counts[sp]
            total = c["total"] if c["total"] > 0 else 1
            comp = c["single"] + c["duplicated"]
            fh.write(
                f"{sp}\t{c['single']}\t{c['duplicated']}\t{c['fragmented']}\t"
                f"{c['missing']}\t{c['total']}\t"
                f"{100.0 * comp / total:.2f}\t"
                f"{100.0 * c['fragmented'] / total:.2f}\t"
                f"{100.0 * c['missing'] / total:.2f}\n"
            )


def _read_tsv(tsv_path):
    """Read an existing busco_completeness.tsv into an all_counts dict."""
    all_counts = {}
    with open(tsv_path) as fh:
        header = fh.readline()  # discard
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            sp = parts[0]
            all_counts[sp] = {
                "single":     int(parts[1]),
                "duplicated": int(parts[2]),
                "fragmented": int(parts[3]),
                "missing":    int(parts[4]),
                "total":      int(parts[5]),
            }
    return all_counts


def _phylo_order(species_list, tree_path):
    """Return species in phylogenetic tip order (ladderized, decreasing)."""
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


def _plot(species_order, all_counts, out_pdf, title="BUSCO Completeness Assessment",
          tree_path=None):
    """Render and save the horizontal stacked bar chart, with optional cladogram sidebar."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from shared_unique_plot import _draw_cladogram

    n = len(species_order)
    y = np.arange(n)
    bar_h = 0.6
    fig_h = max(4, n * 0.65 + 1.8)

    val_map = {
        "Complete (single-copy)": "single",
        "Complete (duplicated)":  "duplicated",
        "Fragmented":             "fragmented",
        "Missing":                "missing",
    }

    if tree_path:
        fig, (ax_tree, ax_bar) = plt.subplots(
            1, 2, figsize=(13, fig_h),
            gridspec_kw={"width_ratios": [1, 3]},
            sharey=True,
        )
        ax_tree.set_ylim(-0.5, n - 0.5)
        _draw_cladogram(ax_tree, species_order, tree_path)
    else:
        fig, ax_bar = plt.subplots(figsize=(10, fig_h))

    left = np.zeros(n)
    for label in BUSCO_ORDER:
        key = val_map[label]
        vals = np.array([all_counts[sp][key] for sp in species_order], dtype=float)
        totals = np.array(
            [all_counts[sp]["total"] if all_counts[sp]["total"] > 0 else 1
             for sp in species_order], dtype=float
        )
        pcts = 100.0 * vals / totals
        ax_bar.barh(y, pcts, left=left, height=bar_h,
                    color=BUSCO_COLOURS[label], label=label)
        for i, (pct, lft) in enumerate(zip(pcts, left)):
            if pct >= 3:
                ax_bar.text(lft + pct / 2, i, f"{BUSCO_ABBREV[label]}:{pct:.1f}%",
                            va="center", ha="center", fontsize=7, color="white",
                            fontweight="bold")
        left = left + pcts

    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(species_order, fontsize=9)
    ax_bar.tick_params(labelleft=True)   # restore labels suppressed by sharey=True
    ax_bar.set_xlabel("% BUSCO groups", fontsize=10)
    ax_bar.set_title(title, fontsize=11, fontweight="bold")
    ax_bar.set_xlim(0, 100)
    ax_bar.set_ylim(-0.5, n - 0.5) # force exact match with tree axis if present
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.legend(frameon=False, fontsize=8,
                  bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    # Phylo-reorder mode: input is existing TSV + tree, output is phylo-ordered PDF only.
    tree_path = getattr(snakemake.input, "tree", None)   # noqa: F821
    if tree_path:
        in_tsv       = snakemake.input.tsv               # noqa: F821
        out_pdf      = snakemake.output.pdf              # noqa: F821
        species_list = list(snakemake.params.species)    # noqa: F821
        os.makedirs(os.path.dirname(os.path.abspath(out_pdf)), exist_ok=True)
        all_counts    = _read_tsv(in_tsv)
        species_order = _phylo_order(species_list, tree_path)
        _plot(species_order, all_counts, out_pdf,
              title="BUSCO Completeness Assessment (phylogenetic order)",
              tree_path=tree_path)
        print(f"[busco_summary] Saved phylo plot → {out_pdf}", flush=True)
        return

    # Normal mode: parse short_summary files, write TSV, plot alphabetically.
    summaries    = list(snakemake.input.summaries)   # noqa: F821
    species_list = list(snakemake.params.species)    # noqa: F821
    out_pdf      = snakemake.output.pdf              # noqa: F821
    out_tsv      = snakemake.output.tsv              # noqa: F821

    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)), exist_ok=True)

    all_counts = {}
    for sp, path in zip(species_list, summaries):
        all_counts[sp] = _parse_busco_summary(path)
        print(f"[busco_summary] Parsed {sp}: {all_counts[sp]}", flush=True)

    species_order = sorted(species_list)
    _write_tsv(out_tsv, species_order, all_counts)
    _plot(species_order, all_counts, out_pdf)
    print(f"[busco_summary] Saved plot → {out_pdf}", flush=True)
    print(f"[busco_summary] Saved TSV  → {out_tsv}", flush=True)


if "snakemake" in dir():
    main()
