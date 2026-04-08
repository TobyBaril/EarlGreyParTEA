"""Saturation plot for EarlGreyParTEA.

Visualises how the number of unique TE families in the pangenome library
accumulates as additional genomes are included in the pipeline, analogous to
pangenome open/closed accumulation curves.

Method
------
The cd-hit-est .clstr file produced during cluster_all_species is parsed to
attribute every cluster (= unique TE family) to its contributing source(s).
Because sequence headers are prefixed with the source name during clustering
({species}_, REPMASKER_{repspec}_, or CUSTOM_), attribution requires no
re-clustering.

N random permutations of genome addition order are run and the cumulative
unique-cluster count is recorded at each step.  Mean ± 95th-percentile CI
across permutations forms the saturation curve.

Fallback (skip_clustering=True)
--------------------------------
When skip_clustering is True in the config, cd-hit-est is not run and the
.clstr file is an empty sentinel.  The script detects this (file size == 0)
and falls back to counting raw sequences per per-genome strained FASTA.
Cross-genome deduplication cannot be applied in this mode.

Outputs
-------
saturation_plot.pdf  : saturation curve (mean + shaded 95% CI)
saturation_data.tsv  : n_genomes | mean_unique_families | ci_lower_95 | ci_upper_95
"""

import csv
import os
import sys
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# cluster_utils lives in the same scripts/ directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_utils import parse_clstr as _parse_clstr_full


def parse_clstr(clstr_file, species_list):
    """Thin wrapper around cluster_utils.parse_clstr for backward compatibility.

    Returns only the two values that saturation_plot.py historically used.
    """
    species_to_clusters, existing_clusters, _, _ = _parse_clstr_full(
        clstr_file, species_list
    )
    return species_to_clusters, existing_clusters


def count_sequences_fasta(fasta_path):
    """Count FASTA sequences (header lines) in a file."""
    count = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


# ---------------------------------------------------------------------------
# Permutation analysis
# ---------------------------------------------------------------------------

def run_permutations(species_to_clusters, existing_clusters, species_list,
                     n_permutations, seed=42):
    """Permutation-based accumulation analysis using cluster attribution.

    Parameters
    ----------
    species_to_clusters : dict[str, set[int]]
    existing_clusters : set[int]
    species_list : list[str]
    n_permutations : int
    seed : int

    Returns
    -------
    results : ndarray, shape (n_permutations, n_species + 1)
        results[p, 0]  = baseline (existing library cluster count)
        results[p, k]  = cumulative unique families after k genomes added
    """
    rng = random.Random(seed)
    n = len(species_list)
    results = np.zeros((n_permutations, n + 1), dtype=int)
    baseline = len(existing_clusters)

    for p in range(n_permutations):
        order = list(species_list)
        rng.shuffle(order)
        discovered = set(existing_clusters)
        results[p, 0] = baseline
        for k, sp in enumerate(order):
            discovered |= species_to_clusters[sp]
            results[p, k + 1] = len(discovered)

    return results


def run_permutations_fallback(species_to_counts, species_list,
                               n_permutations, seed=42):
    """Fallback permutation when skip_clustering=True.

    No cluster membership is available so each sequence is treated as a
    unique family.  Cross-genome redundancy is NOT accounted for.

    Returns
    -------
    results : ndarray, shape (n_permutations, n_species + 1)
    """
    rng = random.Random(seed)
    n = len(species_list)
    results = np.zeros((n_permutations, n + 1), dtype=int)

    for p in range(n_permutations):
        order = list(species_list)
        rng.shuffle(order)
        cumulative = 0
        results[p, 0] = 0
        for k, sp in enumerate(order):
            cumulative += species_to_counts[sp]
            results[p, k + 1] = cumulative

    return results


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_stats(results):
    """Return mean and 95th-percentile CI (2.5 / 97.5) across permutations."""
    means = np.mean(results, axis=0)
    ci_lower = np.percentile(results, 2.5, axis=0)
    ci_upper = np.percentile(results, 97.5, axis=0)
    return means, ci_lower, ci_upper


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_table(table_path, x_values, means, ci_lower, ci_upper):
    """Write TSV with saturation statistics."""
    with open(table_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["n_genomes", "mean_unique_families",
                         "ci_lower_95", "ci_upper_95"])
        for x, m, lo, hi in zip(x_values, means, ci_lower, ci_upper):
            writer.writerow([x, f"{m:.2f}", f"{lo:.2f}", f"{hi:.2f}"])


def make_plot(x_values, means, ci_lower, ci_upper, n_permutations,
              show_baseline, repspec, has_custom, is_fallback, plot_path):
    """Render and save the saturation curve as a PDF."""
    fig, ax = plt.subplots(figsize=(7, 5))

    color = "#1f77b4"
    ax.fill_between(x_values, ci_lower, ci_upper,
                    alpha=0.2, color=color, label="95% CI")
    ax.plot(x_values, means, color=color, linewidth=2,
            marker="o", markersize=5,
            label=f"Mean ({n_permutations} permutations)")

    if show_baseline:
        baseline_val = means[0]
        ax.axhline(baseline_val, linestyle="--", color="grey",
                   linewidth=1, alpha=0.8)
        lib_parts = []
        if repspec:
            lib_parts.append(f"RepeatMasker ({repspec})")
        if has_custom:
            lib_parts.append("custom library")
        annotation = (
            f"Existing library baseline: {int(round(baseline_val))} families"
            + (f"\n({', '.join(lib_parts)})" if lib_parts else "")
        )
        ax.annotate(
            annotation,
            xy=(x_values[0], baseline_val),
            xytext=(0.5, 0.08),
            textcoords="axes fraction",
            fontsize=8,
            color="grey",
            ha="center",
            arrowprops=dict(arrowstyle="->", color="grey", lw=0.8),
        )

    if is_fallback:
        ax.text(
            0.98, 0.02,
            "Warning: skip_clustering=True\nCross-genome deduplication not applied",
            transform=ax.transAxes, fontsize=7, color="darkorange",
            ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                      edgecolor="darkorange", alpha=0.85),
        )

    ax.set_xlabel("Genomes added", fontsize=12)
    ax.set_ylabel("Cumulative unique TE families", fontsize=12)
    ax.set_title("TE Family Library Saturation", fontsize=13, fontweight="bold")
    ax.set_xticks(x_values)
    ax.set_xlim(left=x_values[0] - 0.3, right=x_values[-1] + 0.3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    fig.tight_layout()
    fig.savefig(plot_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point (called as a Snakemake script)
# ---------------------------------------------------------------------------

def main():
    clstr_file = snakemake.input.clstr                       # noqa: F821
    strained_files = list(snakemake.input.strained)          # noqa: F821
    species_list = list(snakemake.params.species)            # noqa: F821
    repspec = snakemake.params.repspec or ""                 # noqa: F821
    has_custom = bool(snakemake.params.has_custom)           # noqa: F821
    n_permutations = int(snakemake.params.permutations)      # noqa: F821
    plot_path = snakemake.output.plot                        # noqa: F821
    table_path = snakemake.output.table                      # noqa: F821

    os.makedirs(os.path.dirname(os.path.abspath(plot_path)), exist_ok=True)

    # Detect skip_clustering fallback: empty sentinel .clstr file
    is_fallback = os.path.getsize(clstr_file) == 0

    if is_fallback:
        print(
            "[saturation_plot] WARNING: .clstr file is empty "
            "(skip_clustering=True). Falling back to raw sequence counts "
            "per genome. Cross-genome redundancy is NOT accounted for.",
            flush=True,
        )
        species_to_counts = {
            sp: count_sequences_fasta(fasta)
            for sp, fasta in zip(species_list, strained_files)
        }
        results = run_permutations_fallback(
            species_to_counts, species_list, n_permutations
        )
        show_baseline = False
    else:
        species_to_clusters, existing_clusters = parse_clstr(
            clstr_file, species_list
        )
        results = run_permutations(
            species_to_clusters, existing_clusters,
            species_list, n_permutations,
        )
        show_baseline = bool(existing_clusters)

    means, ci_lower, ci_upper = compute_stats(results)

    # Build x-axis: include x=0 only when there is a meaningful existing baseline
    if show_baseline:
        x_values = list(range(len(species_list) + 1))
        plot_means, plot_ci_lower, plot_ci_upper = means, ci_lower, ci_upper
    else:
        x_values = list(range(1, len(species_list) + 1))
        plot_means = means[1:]
        plot_ci_lower = ci_lower[1:]
        plot_ci_upper = ci_upper[1:]

    write_table(table_path, x_values, plot_means, plot_ci_lower, plot_ci_upper)
    make_plot(
        x_values, plot_means, plot_ci_lower, plot_ci_upper,
        n_permutations=n_permutations,
        show_baseline=show_baseline,
        repspec=repspec,
        has_custom=has_custom,
        is_fallback=is_fallback,
        plot_path=plot_path,
    )

    print(f"[saturation_plot] Saved plot  → {plot_path}", flush=True)
    print(f"[saturation_plot] Saved table → {table_path}", flush=True)


main()
