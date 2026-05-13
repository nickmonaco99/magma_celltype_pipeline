#!/usr/bin/env python3
"""
05_make_plots.py — Generate publication-style plots from annotated MAGMA outputs.

Reads from <snapshot>/annotated/ (produced by 04_annotate_results.py) and writes
figures to <snapshot>/figures/.

Generates three plots:
  - celltype_manhattan.png        : all 461 Siletti clusters, -log10(P) vs
                                    cluster ID, colored by supercluster, with
                                    Bonferroni threshold line.
  - celltype_forest_top20.png     : top 20 clusters with β ± SE error bars,
                                    P-values annotated, Bonferroni hits bolded.
  - genes_manhattan_annotated.png : gene-level Manhattan with GWS-significant
                                    gene symbols labeled.

The original figures in <snapshot>/figures/ (manhattan.png, qqplot.png) are
NOT modified. These are additional, annotated plots.

Usage:
  python pipeline/05_make_plots.py                       # latest snapshot
  python pipeline/05_make_plots.py results/snapshots/X   # explicit snapshot

Idempotent: re-running overwrites the three plots, updates manifest.
"""

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ----------------------------------------------------------------------------
# Path resolution (identical to 04_annotate_results.py for consistency)
# ----------------------------------------------------------------------------

def find_project_root():
    here = Path(__file__).resolve().parent
    if (here.parent / "pipeline" / "00_config.yaml").exists():
        return here.parent
    if (Path.cwd() / "pipeline" / "00_config.yaml").exists():
        return Path.cwd()
    raise RuntimeError("Cannot find project root.")


def find_latest_snapshot(snapshots_dir):
    snaps = sorted(d for d in snapshots_dir.iterdir() if d.is_dir())
    if not snaps:
        raise RuntimeError(f"No snapshots in {snapshots_dir}")
    return snaps[-1]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# Siletti 2023 supercluster ordering (matches Duncan 2025 Fig 3 legend).
# The 461 L2 clusters group hierarchically into these 31 superclusters, in
# this anatomical/developmental order:
#   1     : Miscellaneous (mix of non-neurons + rare neurons)
#   2-11  : Non-neurons (glia, vascular, ependymal)
#   12-21 : Cortical/hippocampal/striatal excitatory
#   22    : Splatter (subcortical neuromodulatory: HY, MB, HB, mammillary)
#   23-25 : Cortical inhibitory interneurons
#   26-28 : Cerebellar / rhombic lip
#   29-31 : Thalamic, mammillary, midbrain
SILETTI_SUPERCLUSTER_ORDER = [
    "Miscellaneous",
    "Microglia",
    "Vascular",
    "Fibroblast",
    "Oligodendrocyte precursor",
    "Committed oligodendrocyte precursor",
    "Oligodendrocyte",
    "Bergmann glia",
    "Astrocyte",
    "Ependymal",
    "Choroid plexus",
    "Deep-layer near-projecting",
    "Deep-layer corticothalamic and 6b",
    "Hippocampal CA1-3",
    "Upper-layer intratelencephalic",
    "Deep-layer intratelencephalic",
    "Amygdala excitatory",
    "Hippocampal CA4",
    "Hippocampal dentate gyrus",
    "Medium spiny neuron",
    "Eccentric medium spiny neuron",
    "Splatter",
    "MGE interneuron",
    "LAMP5-LHX6 and Chandelier",
    "CGE interneuron",
    "Upper rhombic lip",
    "Cerebellar inhibitory",
    "Lower rhombic lip",
    "Mammillary body",
    "Thalamic excitatory",
    "Midbrain-derived inhibitory",
]


def build_siletti_palette():
    """
    Build a 31-color palette matching Duncan Fig 3 style:
      - #1   : distinctive brown for Miscellaneous
      - #2-11: grayscale gradient (dark to light) for non-neurons
      - #12-31: rainbow (red -> orange -> yellow -> green -> cyan -> blue -> purple -> pink)
                for neuronal superclusters
    Returns a list of RGBA tuples in Siletti supercluster ID order (1..31).
    """
    colors = []
    # #1 Miscellaneous
    colors.append((0.55, 0.30, 0.10, 1.0))
    # #2-11 Non-neurons: grayscale
    for c in plt.cm.gray(np.linspace(0.20, 0.85, 10)):
        colors.append(tuple(c))
    # #12-31 Neurons: rainbow (20 colors)
    for c in plt.cm.gist_rainbow(np.linspace(0.0, 0.92, 20)):
        colors.append(tuple(c))
    return colors


def get_supercluster_palette(superclusters):
    """
    Returns dict {supercluster_name: {'color': rgba, 'order': int (1-indexed)}}
    in Siletti's published ordering. Unknown superclusters get assigned a
    fallback gray and appended after #31.
    """
    palette = build_siletti_palette()
    out = {}
    for i, name in enumerate(SILETTI_SUPERCLUSTER_ORDER):
        out[name] = {"color": palette[i], "order": i + 1}
    unique_in_data = set(superclusters.dropna().unique())
    unknown = sorted(unique_in_data - set(SILETTI_SUPERCLUSTER_ORDER))
    for j, name in enumerate(unknown):
        out[name] = {"color": (0.5, 0.5, 0.5, 1.0), "order": 32 + j}
    return out


def chrom_to_numeric(chr_series):
    """Convert CHR column to numeric; X -> 23."""
    return pd.to_numeric(
        chr_series.astype(str).replace({"X": "23", "Y": "24", "MT": "25"}),
        errors="coerce",
    )


# ----------------------------------------------------------------------------
# Plot 1: Cell-type Manhattan
# ----------------------------------------------------------------------------

def plot_celltype_manhattan(ct_df, bonferroni, gwas_name, out_path):
    """
    Cell-type Manhattan plot mirroring Duncan 2025 Fig 2 / Fig 3 style:
      - X-axis: Siletti L2 cluster ID (0..460), which already groups clusters
                by supercluster in anatomical/developmental order.
      - Y-axis: -log10(P), with a y-min of 0 and y-max floored at 5 so the
                Bonferroni line is visible even when no clusters reach it.
      - Colors: Siletti's published supercluster ordering (1..31), mapped to
                a curated rainbow (browns/grays -> warm -> cool -> pink),
                NOT alphabetical.
      - Visual weighting: clusters with P>0.05 are faint (alpha ~0.2),
                P<0.05 medium, P<0.001 fully saturated. Bigger markers for
                stronger signals. This makes the eye focus on the signal,
                not the noise.
      - Bonferroni line in red dashed, suggestive line (P<1e-3) in light gray.
      - Legend below the plot, 4 columns, with Siletti supercluster IDs.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    df = ct_df.copy().sort_values("cluster_id").reset_index(drop=True)
    df["minus_log10_P"] = -np.log10(df["P"])

    palette = get_supercluster_palette(df["supercluster"])

    # Y-axis range: always show 0 -> max(actual_max+0.7, 5.0) for breathing room
    y_max_data = df["minus_log10_P"].max()
    y_max_plot = max(y_max_data + 0.7, 5.0)

    # Per-point alpha + size based on significance tier
    # P > 0.05            : faint (noise zone)
    # 0.001 < P < 0.05    : moderate
    # 1e-5 < P < 0.001    : strong
    # P < 1e-5            : top-tier
    def _alpha(p):
        if p < 1e-5: return 1.00
        if p < 1e-3: return 0.90
        if p < 0.05: return 0.55
        return 0.22

    def _size(p):
        if p < 1e-5: return 50
        if p < 1e-3: return 32
        if p < 0.05: return 20
        return 14

    df["_alpha"] = df["P"].apply(_alpha)
    df["_size"] = df["P"].apply(_size)

    # Build per-point RGBA colors (apply per-point alpha to supercluster color)
    def _rgba(row):
        info = palette.get(row["supercluster"])
        if info is None:
            return (0.5, 0.5, 0.5, row["_alpha"])
        c = info["color"]
        return (c[0], c[1], c[2], row["_alpha"])

    point_colors = df.apply(_rgba, axis=1).tolist()

    fig, ax = plt.subplots(figsize=(14, 6.8))

    # Noise-zone shading (P > 0.05) — very subtle
    ax.axhspan(0, -np.log10(0.05), color="black", alpha=0.025, zorder=0)

    # Suggestive line (P < 1e-3) — light gray reference
    sug_p = 1e-3
    if -np.log10(sug_p) < y_max_plot:
        ax.axhline(
            -np.log10(sug_p), color="gray", linestyle=":", linewidth=0.8,
            alpha=0.6, zorder=1,
        )
        ax.text(
            462, -np.log10(sug_p), f" Suggestive (P < {sug_p:.0e})",
            fontsize=7.5, color="gray", va="center", ha="left", zorder=1,
        )

    # Bonferroni line
    ax.axhline(
        -np.log10(bonferroni), color="#C0392B", linestyle="--", linewidth=1.1,
        zorder=1,
    )
    ax.text(
        462, -np.log10(bonferroni),
        f" Bonferroni (P < {bonferroni:.2e})",
        fontsize=8, color="#C0392B", va="center", ha="left",
        fontweight="bold", zorder=1,
    )

    # Single scatter call with per-point color/alpha/size
    ax.scatter(
        df["cluster_id"], df["minus_log10_P"],
        c=point_colors, s=df["_size"],
        edgecolors="white", linewidths=0.3, zorder=2,
    )

    # Annotate clusters passing Bonferroni
    sig = df[df["P"] < bonferroni].sort_values("P")
    for _, row in sig.iterrows():
        ax.annotate(
            row["cluster_name"],
            (row["cluster_id"], -np.log10(row["P"])),
            xytext=(7, 7), textcoords="offset points",
            fontsize=10, fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
            zorder=5,
        )

    # If fewer than 3 sig, also label top 3 by P (sub-Bonferroni context)
    if len(sig) < 3:
        top3 = df.nsmallest(3, "P")
        for _, row in top3.iterrows():
            if row["P"] >= bonferroni:
                ax.annotate(
                    row["cluster_name"],
                    (row["cluster_id"], -np.log10(row["P"])),
                    xytext=(7, 7), textcoords="offset points",
                    fontsize=8.5, color="#444444", style="italic",
                    zorder=5,
                )

    # Axes
    ax.set_xlabel("Siletti cluster ID", fontsize=11)
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=11)
    ax.set_title(
        f"Cell-type analysis — {gwas_name}\n"
        f"461 Siletti L2 clusters  •  Bonferroni P < {bonferroni:.2e}",
        fontsize=12, pad=12,
    )
    ax.set_xlim(-8, 480)
    ax.set_ylim(-0.1, y_max_plot)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.3, zorder=0)
    ax.tick_params(axis="both", which="major", labelsize=10)

    # Legend below plot: Siletti supercluster ID + name, in numeric order,
    # only for superclusters actually present in the data.
    present = [
        sc for sc in SILETTI_SUPERCLUSTER_ORDER
        if sc in df["supercluster"].dropna().unique()
    ]
    extra = [
        sc for sc in df["supercluster"].dropna().unique()
        if sc not in SILETTI_SUPERCLUSTER_ORDER
    ]
    legend_handles = []
    for sc in present + extra:
        info = palette[sc]
        legend_handles.append(
            Line2D(
                [0], [0], marker="o", color="w",
                markerfacecolor=info["color"], markersize=8,
                markeredgecolor="white", markeredgewidth=0.4,
                label=f"{info['order']:>2}  {sc}",
            )
        )

    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=4,
        fontsize=8,
        framealpha=0.95,
        title="Superclusters (Siletti 2023 / Duncan 2025 numbering)",
        title_fontsize=9,
        handletextpad=0.4,
        columnspacing=1.0,
        borderaxespad=0.5,
    )

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Plot 2: Forest plot of top 20
# ----------------------------------------------------------------------------

def plot_celltype_forest(ct_df, bonferroni, gwas_name, out_path, top_n=20):
    """
    Horizontal forest plot of top N clusters by P, with β ± SE error bars.
    Bonferroni-significant rows have bold labels.
    """
    top = ct_df.nsmallest(top_n, "P").reset_index(drop=True)
    # Reverse so smallest P (most significant) appears at TOP of plot
    top = top.iloc[::-1].reset_index(drop=True)

    palette = get_supercluster_palette(top["supercluster"])
    colors = [
        palette[sc]["color"] if sc in palette else (0.5, 0.5, 0.5, 1.0)
        for sc in top["supercluster"]
    ]

    fig, ax = plt.subplots(figsize=(11, max(6, 0.4 * top_n)))
    y_pos = np.arange(len(top))

    # Error bars + points
    ax.errorbar(
        top["BETA"], y_pos,
        xerr=top["SE"],
        fmt="o", markersize=7,
        ecolor="gray", capsize=3, capthick=0.8, elinewidth=0.8,
        markerfacecolor="none", markeredgecolor="gray", markeredgewidth=0.8,
        zorder=1,
    )
    # Colored markers on top
    ax.scatter(top["BETA"], y_pos, c=colors, s=60, zorder=2,
               edgecolors="black", linewidths=0.6)

    ax.axvline(0, color="gray", linewidth=0.5, zorder=0)

    # Y-axis labels: cluster_name + supercluster (with Siletti ID prefix)
    labels = []
    for _, row in top.iterrows():
        label = row["cluster_name"]
        if pd.notna(row.get("supercluster")) and row["supercluster"] in palette:
            sc_num = palette[row["supercluster"]]["order"]
            label = f"{label}  [#{sc_num} {row['supercluster']}]"
        elif pd.notna(row.get("supercluster")):
            label = f"{label}  [{row['supercluster']}]"
        labels.append(label)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)

    # Bold significant labels
    for tick, (_, row) in zip(ax.get_yticklabels(), top.iterrows()):
        if row["P"] < bonferroni:
            tick.set_fontweight("bold")
            tick.set_color("red")

    # P-value annotations to the right
    xmax = ax.get_xlim()[1]
    xpos = xmax + 0.02 * (xmax - ax.get_xlim()[0])
    for i, (_, row) in enumerate(top.iterrows()):
        marker = " *" if row["P"] < bonferroni else ""
        ax.text(
            xpos, i, f"P = {row['P']:.2e}{marker}",
            fontsize=8, va="center", family="monospace",
        )

    ax.set_xlabel(r"$\beta$ (effect size on gene Z-stat)")
    ax.set_title(
        f"Top {top_n} cell-type associations — {gwas_name}\n"
        f"Bold red = passes Bonferroni (P < {bonferroni:.2e})"
    )
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Plot 3: Annotated gene-level Manhattan
# ----------------------------------------------------------------------------

def plot_gene_manhattan(genes_df, gwas_name, out_path, gws_threshold=5e-8):
    """
    Gene-level Manhattan with GWS-significant gene symbols labeled.
    Chromosomes alternate color.
    """
    df = genes_df.copy()
    df["CHR_num"] = chrom_to_numeric(df["CHR"])
    df = df.dropna(subset=["CHR_num", "P"])
    df = df.sort_values(["CHR_num", "START"]).reset_index(drop=True)
    df["x_pos"] = np.arange(len(df))
    df["minus_log10_P"] = -np.log10(df["P"])

    fig, ax = plt.subplots(figsize=(14, 5.5))

    # Plot each chromosome with alternating color
    chrom_midpoints = {}
    chroms_present = sorted(df["CHR_num"].unique())
    for i, chrom in enumerate(chroms_present):
        mask = df["CHR_num"] == chrom
        color = "#3674B5" if i % 2 == 0 else "#F18F01"
        ax.scatter(
            df.loc[mask, "x_pos"], df.loc[mask, "minus_log10_P"],
            color=color, s=8, alpha=0.7, edgecolors="none",
        )
        chrom_midpoints[chrom] = df.loc[mask, "x_pos"].median()

    # GWS threshold line
    ax.axhline(
        -np.log10(gws_threshold), color="red", linestyle="--", linewidth=1.0,
        label=f"GWS (P < {gws_threshold:.0e})",
    )

    # Annotate genes passing GWS
    gws = df[df["P"] < gws_threshold].copy()
    if len(gws) > 0:
        # Avoid label overlap: when multiple genes are in same chromosome region,
        # only label the strongest signal per ~500kb window
        gws = gws.sort_values("P")
        labeled_chr_pos = []  # (chrom, position) pairs already labeled
        for _, row in gws.iterrows():
            chrom, start = row["CHR_num"], row["START"]
            too_close = any(
                c == chrom and abs(p - start) < 5e5
                for c, p in labeled_chr_pos
            )
            if too_close:
                continue
            symbol = row.get("gene_symbol", None)
            if pd.isna(symbol) or symbol is None or str(symbol).strip() == "":
                symbol = f"E:{row['GENE']}"
            ax.annotate(
                symbol,
                (row["x_pos"], row["minus_log10_P"]),
                xytext=(4, 6), textcoords="offset points",
                fontsize=8.5, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
            )
            labeled_chr_pos.append((chrom, start))

    # Chromosome ticks
    ax.set_xticks(list(chrom_midpoints.values()))
    ax.set_xticklabels([
        "X" if c == 23 else "Y" if c == 24 else str(int(c))
        for c in chrom_midpoints.keys()
    ])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title(
        f"Gene-level Manhattan — {gwas_name}\n"
        f"{(df['P'] < gws_threshold).sum()} genes pass GWS (P < {gws_threshold:.0e})"
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Snapshot updates
# ----------------------------------------------------------------------------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def regenerate_checksums(snapshot_dir):
    out = snapshot_dir / "checksums.sha256"
    lines = []
    for f in sorted(snapshot_dir.rglob("*")):
        if f.is_file() and f.name != "checksums.sha256":
            lines.append(f"{sha256_file(f)}  {f.relative_to(snapshot_dir)}")
    out.write_text("\n".join(lines) + "\n")


def update_manifest_plots(snapshot_dir, plot_info):
    path = snapshot_dir / "manifest.json"
    with open(path) as f:
        manifest = json.load(f)
    manifest["plots"] = plot_info
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("snapshot_path", nargs="?", default=None)
    parser.add_argument("--gws", type=float, default=5e-8,
                        help="Genome-wide significance threshold (default: 5e-8)")
    args = parser.parse_args()

    project_root = find_project_root()

    if args.snapshot_path:
        snap_dir = Path(args.snapshot_path).resolve()
    else:
        snap_dir = find_latest_snapshot(project_root / "results" / "snapshots")
    print(f"[plots] Snapshot: {snap_dir.relative_to(project_root)}")

    ct_path = snap_dir / "annotated" / "celltype_annotated.csv"
    g_path = snap_dir / "annotated" / "genes_annotated.csv"
    manifest_path = snap_dir / "manifest.json"

    for p in (ct_path, g_path, manifest_path):
        if not p.exists():
            sys.exit(f"[plots] ERROR: required input not found: {p}\n"
                     f"        Did you run 04_annotate_results.py first?")

    with open(manifest_path) as f:
        manifest = json.load(f)
    gwas_name = manifest.get("gwas", {}).get("name", "GWAS")
    celltype_bonf = manifest.get("celltype_bonferroni", 0.05 / 461)

    figures_dir = snap_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    ct = pd.read_csv(ct_path)
    g = pd.read_csv(g_path)
    print(f"[plots] Loaded {len(ct)} clusters, {len(g)} genes")

    print(f"[plots] Plot 1/3: celltype Manhattan...")
    p1 = figures_dir / "celltype_manhattan.png"
    plot_celltype_manhattan(ct, celltype_bonf, gwas_name, p1)
    print(f"        wrote {p1.relative_to(snap_dir)}")

    print(f"[plots] Plot 2/3: celltype forest (top 20)...")
    p2 = figures_dir / "celltype_forest_top20.png"
    plot_celltype_forest(ct, celltype_bonf, gwas_name, p2, top_n=20)
    print(f"        wrote {p2.relative_to(snap_dir)}")

    print(f"[plots] Plot 3/3: annotated gene Manhattan...")
    p3 = figures_dir / "genes_manhattan_annotated.png"
    plot_gene_manhattan(g, gwas_name, p3, gws_threshold=args.gws)
    print(f"        wrote {p3.relative_to(snap_dir)}")

    plot_info = {
        "plotted_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "gws_threshold": args.gws,
        "celltype_bonferroni": celltype_bonf,
        "files": {
            "celltype_manhattan": str(p1.relative_to(snap_dir)),
            "celltype_forest_top20": str(p2.relative_to(snap_dir)),
            "genes_manhattan_annotated": str(p3.relative_to(snap_dir)),
        },
    }
    update_manifest_plots(snap_dir, plot_info)
    print(f"[plots] Updated manifest.json (added 'plots' key)")

    print(f"[plots] Regenerating checksums.sha256...")
    regenerate_checksums(snap_dir)
    print(f"[plots] Done.")


if __name__ == "__main__":
    main()
