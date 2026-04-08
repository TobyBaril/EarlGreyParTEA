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


def main():
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

    n = len(species_order)
    y = np.arange(n)
    bar_h = 0.6

    val_map = {
        "Complete (single-copy)": "single",
        "Complete (duplicated)":  "duplicated",
        "Fragmented":             "fragmented",
        "Missing":                "missing",
    }

    fig, ax = plt.subplots(figsize=(10, max(4, n * 0.65 + 1.8)))

    left = np.zeros(n)
    for label in BUSCO_ORDER:
        key = val_map[label]
        vals = np.array([all_counts[sp][key] for sp in species_order], dtype=float)
        totals = np.array(
            [all_counts[sp]["total"] if all_counts[sp]["total"] > 0 else 1
             for sp in species_order], dtype=float
        )
        pcts = 100.0 * vals / totals
        ax.barh(y, pcts, left=left, height=bar_h,
                color=BUSCO_COLOURS[label], label=label)
        # Label inside segment if wide enough
        for i, (pct, lft) in enumerate(zip(pcts, left)):
            if pct >= 3:
                ax.text(lft + pct / 2, i, f"{BUSCO_ABBREV[label]}:{pct:.1f}%",
                        va="center", ha="center", fontsize=7, color="white",
                        fontweight="bold")
        left = left + pcts

    ax.set_yticks(y)
    ax.set_yticklabels(species_order, fontsize=9)
    ax.set_xlabel("% BUSCO groups", fontsize=10)
    ax.set_title("BUSCO Completeness Assessment", fontsize=11, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[busco_summary] Saved plot → {out_pdf}", flush=True)
    print(f"[busco_summary] Saved TSV  → {out_tsv}", flush=True)


if "snakemake" in dir():
    main()
