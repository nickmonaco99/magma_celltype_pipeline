#!/usr/bin/env python3
"""
04_annotate_results.py — Annotate MAGMA outputs with biological metadata.

Reads a snapshot directory (default: latest under results/snapshots/) and
writes annotated CSV tables to <snapshot>/annotated/:

  - celltype_annotated.csv : 461 clusters joined with Siletti supercluster,
                             cell class, neurotransmitter, anatomical regions,
                             and top marker genes.
  - genes_annotated.csv    : all genes joined with HGNC symbols from
                             NCBI37.3.gene.loc (column 6, if present).

The original MAGMA outputs in <snapshot>/gene_level/ and <snapshot>/cell_type/
are NOT modified. Annotations are derived, written to a new subdirectory.

Inputs (paths read from pipeline/00_config.yaml):
  - <snapshot>/gene_level/genes.out
  - <snapshot>/cell_type/cell_type.gsa.out
  - auxfiles/NCBI37.3.gene.loc
  - auxfiles/science.add7046_table_s3.xlsx        (Siletti 2023 Table S3)

Outputs:
  - <snapshot>/annotated/celltype_annotated.csv
  - <snapshot>/annotated/genes_annotated.csv
  - <snapshot>/manifest.json                       (adds top-level 'annotation' key)
  - <snapshot>/checksums.sha256                    (regenerated for full snapshot)

Usage:
  python pipeline/04_annotate_results.py                       # latest snapshot
  python pipeline/04_annotate_results.py results/snapshots/X   # explicit snapshot

Idempotent: safe to re-run; overwrites previous annotated/ contents.
"""

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


# ----------------------------------------------------------------------------
# Path resolution
# ----------------------------------------------------------------------------

def find_project_root():
    """Walk up from this script's location to find pipeline/00_config.yaml."""
    here = Path(__file__).resolve().parent
    candidate = here.parent / "pipeline" / "00_config.yaml"
    if candidate.exists():
        return here.parent
    # Fallback: maybe script is being run from project root
    if (Path.cwd() / "pipeline" / "00_config.yaml").exists():
        return Path.cwd()
    raise RuntimeError(
        "Cannot find project root. Expected pipeline/00_config.yaml relative "
        "to this script or current working directory."
    )


def load_config(project_root):
    with open(project_root / "pipeline" / "00_config.yaml") as f:
        return yaml.safe_load(f)


def find_latest_snapshot(snapshots_dir):
    if not snapshots_dir.exists():
        raise RuntimeError(f"No snapshots directory at {snapshots_dir}")
    snapshots = sorted(d for d in snapshots_dir.iterdir() if d.is_dir())
    if not snapshots:
        raise RuntimeError(f"No snapshots found under {snapshots_dir}")
    return snapshots[-1]


# ----------------------------------------------------------------------------
# Parsers
# ----------------------------------------------------------------------------

def parse_siletti_annotations(xlsx_path):
    """
    Load Siletti 2023 Table S3 and return a clean DataFrame.

    The Cluster ID column is float (0.0, 1.0, ...) in the source. We cast to
    int and add a `cluster_label` column matching MAGMA's output format
    ("Cluster0", "Cluster1", ...).
    """
    df = pd.read_excel(xlsx_path)
    df = df.dropna(subset=["Cluster ID"]).copy()
    df["Cluster ID"] = df["Cluster ID"].astype(int)
    df["cluster_label"] = "Cluster" + df["Cluster ID"].astype(str)
    rename_map = {
        "Cluster ID": "cluster_id",
        "Cluster name": "cluster_name",
        "Supercluster": "supercluster",
        "Class auto-annotation": "class_auto",
        "Neurotransmitter auto-annotation": "neurotransmitter",
        "Top three regions": "top_regions",
        "Top three dissections": "top_dissections",
        "Top Enriched Genes": "top_genes",
        "Number of cells": "n_cells",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    return df


def parse_magma_celltype(path):
    """Parse MAGMA cell-type output (.gsa.out): whitespace-sep, '#' comments."""
    df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
    df["VARIABLE"] = df["VARIABLE"].astype(str).str.strip()
    return df


def parse_magma_genes(path):
    """Parse MAGMA gene-level output (genes.out)."""
    df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
    df["GENE"] = df["GENE"].astype(str)
    return df


def parse_gene_loc(path):
    """
    Parse MAGMA gene location file.

    Standard MAGMA format: gene_id  chr  start  stop  strand  [gene_symbol]
    Column 6 (gene symbol) is optional in MAGMA's spec but present in the
    canonical NCBI37.3.gene.loc distributed with MAGMA aux files.
    """
    # Peek at first non-empty line to detect column count
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                ncols = len(line.split())
                break
        else:
            raise ValueError(f"Empty gene.loc file: {path}")

    if ncols >= 6:
        names = ["entrez_id", "chr", "start", "stop", "strand", "gene_symbol"]
        usecols = [0, 1, 2, 3, 4, 5]
    elif ncols == 5:
        names = ["entrez_id", "chr", "start", "stop", "strand"]
        usecols = [0, 1, 2, 3, 4]
    else:
        raise ValueError(f"Unexpected gene.loc column count: {ncols}")

    df = pd.read_csv(
        path, sep=r"\s+", header=None, names=names, usecols=usecols,
        dtype={"entrez_id": str}, engine="python",
    )
    if "gene_symbol" not in df.columns:
        df["gene_symbol"] = None
    return df[["entrez_id", "gene_symbol"]]


# ----------------------------------------------------------------------------
# Merging
# ----------------------------------------------------------------------------

def annotate_celltype(gsa_df, siletti_df):
    merged = gsa_df.merge(
        siletti_df, left_on="VARIABLE", right_on="cluster_label", how="left"
    )
    cols = [
        "VARIABLE", "cluster_id", "cluster_name", "supercluster", "class_auto",
        "neurotransmitter", "top_regions", "top_dissections", "top_genes",
        "n_cells", "NGENES", "BETA", "BETA_STD", "SE", "P",
    ]
    cols = [c for c in cols if c in merged.columns]
    return merged[cols].sort_values("P").reset_index(drop=True)


def annotate_genes(genes_df, gene_loc_df):
    merged = genes_df.merge(
        gene_loc_df, left_on="GENE", right_on="entrez_id", how="left"
    )
    cols_first = [
        "GENE", "gene_symbol", "CHR", "START", "STOP",
        "NSNPS", "NPARAM", "N", "ZSTAT", "P",
    ]
    cols_first = [c for c in cols_first if c in merged.columns]
    other = [c for c in merged.columns if c not in cols_first and c != "entrez_id"]
    return merged[cols_first + other].sort_values("P").reset_index(drop=True)


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
            rel = f.relative_to(snapshot_dir)
            lines.append(f"{sha256_file(f)}  {rel}")
    out.write_text("\n".join(lines) + "\n")


def update_manifest(snapshot_dir, annotation_info):
    path = snapshot_dir / "manifest.json"
    with open(path) as f:
        manifest = json.load(f)
    manifest["annotation"] = annotation_info
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "snapshot_path", nargs="?", default=None,
        help="Path to snapshot directory (default: latest under results/snapshots/)",
    )
    args = parser.parse_args()

    project_root = find_project_root()
    config = load_config(project_root)

    if args.snapshot_path:
        snap_dir = Path(args.snapshot_path).resolve()
    else:
        snap_dir = find_latest_snapshot(project_root / "results" / "snapshots")

    print(f"[annotate] Snapshot: {snap_dir.relative_to(project_root)}")

    siletti_path = project_root / config["paths"]["siletti_annotations"]
    gene_loc_path = project_root / config["paths"]["gene_loc"]
    gsa_path = snap_dir / "cell_type" / "cell_type.gsa.out"
    genes_path = snap_dir / "gene_level" / "genes.out"

    for p in (siletti_path, gene_loc_path, gsa_path, genes_path):
        if not p.exists():
            sys.exit(f"[annotate] ERROR: required input not found: {p}")

    out_dir = snap_dir / "annotated"
    out_dir.mkdir(exist_ok=True)

    # Cell-type annotation
    print(f"[annotate] Reading cell-type MAGMA output...")
    gsa = parse_magma_celltype(gsa_path)
    print(f"[annotate] Reading Siletti Table S3...")
    siletti = parse_siletti_annotations(siletti_path)
    print(f"[annotate] Merging cell-type results with Siletti annotations...")
    ct_annotated = annotate_celltype(gsa, siletti)
    ct_out = out_dir / "celltype_annotated.csv"
    ct_annotated.to_csv(ct_out, index=False)
    n_unmatched_ct = int(ct_annotated.get("cluster_id", pd.Series([])).isna().sum())
    print(f"[annotate]   wrote {ct_out.relative_to(snap_dir)} "
          f"({len(ct_annotated)} rows, {n_unmatched_ct} unmatched clusters)")

    # Gene-level annotation
    print(f"[annotate] Reading gene-level MAGMA output...")
    genes = parse_magma_genes(genes_path)
    print(f"[annotate] Reading gene location file...")
    gene_loc = parse_gene_loc(gene_loc_path)
    has_symbols = gene_loc["gene_symbol"].notna().any()
    print(f"[annotate] Gene symbols present in gene.loc: {has_symbols}")
    print(f"[annotate] Merging gene-level results with gene symbols...")
    g_annotated = annotate_genes(genes, gene_loc)
    g_out = out_dir / "genes_annotated.csv"
    g_annotated.to_csv(g_out, index=False)
    n_unmatched_g = (
        int(g_annotated["gene_symbol"].isna().sum()) if has_symbols else None
    )
    print(f"[annotate]   wrote {g_out.relative_to(snap_dir)} "
          f"({len(g_annotated)} rows"
          f"{f', {n_unmatched_g} unmatched symbols' if has_symbols else ''})")

    # Manifest update
    annotation_info = {
        "annotated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "inputs": {
            "siletti_annotations": str(siletti_path.relative_to(project_root)),
            "gene_loc": str(gene_loc_path.relative_to(project_root)),
        },
        "outputs": {
            "celltype_annotated": str(ct_out.relative_to(snap_dir)),
            "genes_annotated": str(g_out.relative_to(snap_dir)),
        },
        "gene_symbol_mapping_available": bool(has_symbols),
        "n_clusters_total": int(len(ct_annotated)),
        "n_clusters_unmatched": n_unmatched_ct,
        "n_genes_total": int(len(g_annotated)),
        "n_genes_unmatched_symbol": n_unmatched_g,
    }
    update_manifest(snap_dir, annotation_info)
    print(f"[annotate] Updated manifest.json (added 'annotation' key)")

    # Checksum regeneration
    print(f"[annotate] Regenerating checksums.sha256...")
    regenerate_checksums(snap_dir)
    print(f"[annotate] Done.")


if __name__ == "__main__":
    main()
