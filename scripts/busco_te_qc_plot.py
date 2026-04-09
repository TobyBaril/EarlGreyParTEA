"""busco_te_qc_plot.py — BUSCO completeness vs TE content scatter plot.

Visualises the relationship between genome assembly quality (BUSCO completeness
%) and repeat content (total TE coverage %) across all species, using TE class
breakdown as colour (same EarlGrey palette as shared_unique_plot.py).

The scatter plot has:
  x-axis: BUSCO completeness (%)  = (complete single-copy + duplicated) / total
  y-axis: total TE coverage (% genome)
  point colour: dominant TE class (the class with the most coverage)
  point size: genome size (proportional, with legend)
  point label: species name

Called as a Snakemake script (script: directive).

snakemake.input.busco_summaries : list of BUSCO short_summary.txt paths
snakemake.input.cov_tsv         : shared_unique_coverage.tsv (from shared_unique_plot)
snakemake.params.species        : list[str]
snakemake.output.pdf            : output scatter plot PDF
snakemake.output.tsv            : combined QC table TSV
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shared_unique_plot import CLASS_PALETTE, CLASS_ORDER
from busco_summary_plot import _parse_busco_summary


def _load_coverage_tsv(path, species_list):
    """Read shared_unique_coverage.tsv and return per-species class coverage %.

    Returns dict{ species -> { 'total_pct': float, 'dominant_class': str,
                                'class_pcts': {cls: float} } }
    """
    sp_set = set(species_list)
    results = {}
    cls_col_prefix_shared = "shared_"
    cls_col_prefix_unique = "unique_"

    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    if not lines:
        return results

    header = lines[0].rstrip("\n").split("\t")
    for line in lines[1:]:
        parts = line.rstrip("\n").split("\t")
        row = dict(zip(header, parts))
        sp = row.get("species", "")
        if sp not in sp_set:
            continue
        gs = float(row.get("genome_size_bp", 1) or 1)
        # Sum shared + unique bp per class
        class_pcts = {}
        for cls in CLASS_ORDER:
            safe_cls = cls.replace(" ", "_")
            s_bp = float(row.get(f"shared_{safe_cls}_bp", 0) or 0)
            u_bp = float(row.get(f"unique_{safe_cls}_bp", 0) or 0)
            class_pcts[cls] = 100.0 * (s_bp + u_bp) / gs if gs > 0 else 0.0
        total_pct = sum(class_pcts.values())
        dominant = max(class_pcts, key=class_pcts.get) if class_pcts else "Unclassified"
        results[sp] = {
            "total_pct":    total_pct,
            "dominant_class": dominant,
            "class_pcts":   class_pcts,
            "genome_size":  gs,
        }
    return results


def main():
    busco_summaries = list(snakemake.input.busco_summaries)  # noqa: F821
    cov_tsv         = snakemake.input.cov_tsv                # noqa: F821
    species_list    = list(snakemake.params.species)         # noqa: F821
    out_pdf         = snakemake.output.pdf                   # noqa: F821
    out_tsv         = snakemake.output.tsv                   # noqa: F821

    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)), exist_ok=True)

    # Parse BUSCO summaries
    busco_data = {}
    for sp, path in zip(species_list, busco_summaries):
        c = _parse_busco_summary(path)
        total = c["total"] if c["total"] > 0 else 1
        comp_pct = 100.0 * (c["single"] + c["duplicated"]) / total
        busco_data[sp] = {"completeness_pct": comp_pct, "raw": c}

    # Load coverage data
    cov_data = _load_coverage_tsv(cov_tsv, species_list)

    # Write combined TSV
    with open(out_tsv, "w") as fh:
        fh.write(
            "species\tbusco_completeness_pct\ttotal_te_pct\t"
            "dominant_te_class\tgenome_size_bp\n"
        )
        for sp in sorted(species_list):
            bc = busco_data.get(sp, {}).get("completeness_pct", float("nan"))
            cd = cov_data.get(sp, {})
            te_pct = cd.get("total_pct", float("nan"))
            dom    = cd.get("dominant_class", "Unclassified")
            gs     = cd.get("genome_size", 0)
            fh.write(f"{sp}\t{bc:.4f}\t{te_pct:.4f}\t{dom}\t{gs}\n")

    # Build scatter data
    sp_plot = [sp for sp in species_list
               if sp in busco_data and sp in cov_data]
    if not sp_plot:
        # Write empty plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes)
        fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
        plt.close(fig)
        return

    x = np.array([busco_data[sp]["completeness_pct"] for sp in sp_plot])
    y = np.array([cov_data[sp]["total_pct"] for sp in sp_plot])
    colours = [CLASS_PALETTE.get(cov_data[sp]["dominant_class"], "#A0A5A9")
               for sp in sp_plot]
    genome_sizes = np.array([cov_data[sp]["genome_size"] for sp in sp_plot])

    # Normalise genome size to point area (min 50, max 400)
    gs_norm = genome_sizes / genome_sizes.max() if genome_sizes.max() > 0 else genome_sizes
    sizes = 50 + 350 * gs_norm

    fig, ax = plt.subplots(figsize=(9, 7))
    sc = ax.scatter(x, y, c=colours, s=sizes, alpha=0.85, edgecolors="#333333",
                    linewidths=0.6, zorder=3)

    # Species labels
    for xi, yi, sp in zip(x, y, sp_plot):
        ax.annotate(
            sp, (xi, yi),
            xytext=(4, 4), textcoords="offset points",
            fontsize=7, color="#333333",
        )

    ax.set_xlabel("BUSCO completeness (%)", fontsize=11)
    ax.set_ylabel("Total TE coverage (% genome)", fontsize=11)
    ax.set_title("Genome Quality vs TE Content", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, linestyle="--", alpha=0.4, zorder=0)

    # Legend: TE classes present
    present_classes = set(cov_data[sp]["dominant_class"] for sp in sp_plot)
    class_patches = [
        mpatches.Patch(facecolor=CLASS_PALETTE[cls], edgecolor="none", label=cls)
        for cls in CLASS_ORDER if cls in present_classes
    ]
    ax.legend(handles=class_patches, title="Dominant TE class",
              frameon=False, fontsize=8, title_fontsize=8,
              bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[busco_te_qc] Saved scatter plot → {out_pdf}", flush=True)
    print(f"[busco_te_qc] Saved TSV          → {out_tsv}", flush=True)


if "snakemake" in dir():
    main()
