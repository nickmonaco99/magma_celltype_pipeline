#!/usr/bin/env python3
"""
snapshot_results.py — Create a timestamped, immutable snapshot of a completed
MAGMA cell-type pipeline run.

Usage:
    python snapshot_results.py <gwas_key>

Example:
    python snapshot_results.py insomnia_ukb

What it captures:
    - The exact 00_config.yaml that produced these results
    - The .genes.out (gene-level p-values, ~1.5 MB)
    - The .gsa.out (cell-type p-values, the actual result, ~30 KB)
    - All MAGMA logs and our run.log
    - The Manhattan and QQ figures
    - A manifest.json with provenance + top hits parsed out
    - SHA-256 checksums for every file

What it does NOT capture (re-generatable from raw data + config):
    - Prepared sumstats (.no_heading, ~300 MB)
    - SNP location file
    - .genes.annot (SNP→gene mapping, ~60 MB)
    - .genes.raw (gene-level raw stats, ~11 MB)

Output: results/snapshots/<gwas_key>_v<N>_<YYYY-MM-DD>/
"""
import sys
import json
import yaml
import shutil
import hashlib
import pathlib
import platform
import subprocess
from datetime import datetime, timezone


def sha256_file(path: pathlib.Path, chunk_size: int = 65536) -> str:
    """Compute SHA-256 of a file."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def get_git_commit(project_root: pathlib.Path) -> str:
    """Return current git commit hash, or 'not-git-tracked' if no repo."""
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            cwd=project_root, capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return 'not-git-tracked'


def parse_gsa_top_hits(gsa_path: pathlib.Path, top_n: int = 20) -> list[dict]:
    """Extract top N hits from a MAGMA .gsa.out file."""
    rows = []
    in_data = False
    with open(gsa_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[0] == 'VARIABLE':
                cols = parts
                in_data = True
                continue
            if in_data and len(parts) == len(cols):
                row = dict(zip(cols, parts))
                rows.append(row)
    # Sort by P (ascending) and return top N
    rows.sort(key=lambda r: float(r['P']))
    return rows[:top_n]


def parse_genes_out_summary(genes_out: pathlib.Path,
                             bonferroni: float = 2.57e-6,
                             gws_strict: float = 5e-8) -> dict:
    """Count significant genes at multiple thresholds."""
    n_total = 0
    n_bonf = 0
    n_gws = 0
    with open(genes_out) as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) < 9:
                continue
            n_total += 1
            try:
                p = float(parts[-1])
                if p < bonferroni: n_bonf += 1
                if p < gws_strict: n_gws += 1
            except ValueError:
                continue
    return {
        'total_genes_tested':       n_total,
        'bonferroni_threshold':     bonferroni,
        'genes_passing_bonferroni': n_bonf,
        'gws_strict_threshold':     gws_strict,
        'genes_passing_gws_strict': n_gws,
    }


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    gwas_key = sys.argv[1]
    pipeline_dir = pathlib.Path(__file__).resolve().parent
    project_root = pipeline_dir.parent

    # ---- Load config ----
    config_path = pipeline_dir / '00_config.yaml'
    with open(config_path) as f:
        config = yaml.safe_load(f)

    if gwas_key not in config['gwas']:
        print(f"ERROR: '{gwas_key}' not in 00_config.yaml gwas: section")
        sys.exit(1)

    gwas = config['gwas'][gwas_key]
    folder = gwas['output_folder']
    results_dir = project_root / config['paths']['results_dir'] / folder

    if not results_dir.exists():
        print(f"ERROR: results dir not found: {results_dir}")
        sys.exit(1)

    # ---- Determine snapshot version ----
    snapshots_root = project_root / config['paths']['results_dir'] / 'snapshots'
    snapshots_root.mkdir(parents=True, exist_ok=True)

    today = datetime.now(timezone.utc).strftime('%Y-%m-%d')
    existing = sorted(snapshots_root.glob(f"{gwas_key}_v*_*"))
    version = 1 + max(
        (int(p.name.split('_v')[1].split('_')[0]) for p in existing),
        default=0
    )
    snapshot_name = f"{gwas_key}_v{version}_{today}"
    snapshot_dir = snapshots_root / snapshot_name
    snapshot_dir.mkdir(parents=True)

    print(f"Creating snapshot: {snapshot_dir.relative_to(project_root)}")

    # ---- Set up subdirs ----
    (snapshot_dir / 'gene_level').mkdir()
    (snapshot_dir / 'cell_type').mkdir()
    (snapshot_dir / 'figures').mkdir()

    # ---- Copy files ----
    file_specs = [
        # (source, dest_within_snapshot)
        (config_path,                                                              '00_config.yaml'),
        (results_dir / 'run.log',                                                  'run.log'),
        (results_dir / f"{gwas_key}.no_heading.step2.genes.out",                  'gene_level/genes.out'),
        (results_dir / f"{gwas_key}.no_heading.step2.log",                        'gene_level/step2.log'),
        (results_dir / f"Siletti_l2_conti-spe_{folder}.gsa.out",                  'cell_type/cell_type.gsa.out'),
        (results_dir / f"Siletti_l2_conti-spe_{folder}.gsa.log",                  'cell_type/cell_type.gsa.log'),
        (project_root / 'figures' / f'01_manhattan_{folder}.png',                 'figures/manhattan.png'),
        (project_root / 'figures' / f'01_qqplot_{folder}.png',                    'figures/qqplot.png'),
    ]

    copied = []
    missing = []
    for src, dest_rel in file_specs:
        dest = snapshot_dir / dest_rel
        if src.exists():
            shutil.copy2(src, dest)
            copied.append((src, dest))
        else:
            missing.append(src)

    print(f"  Copied {len(copied)} files; {len(missing)} expected files missing")
    if missing:
        for m in missing:
            print(f"    MISSING: {m.relative_to(project_root)}")

    # ---- Compute checksums ----
    checksums = []
    for f in sorted(snapshot_dir.rglob('*')):
        if f.is_file() and f.name != 'checksums.sha256':
            rel = f.relative_to(snapshot_dir)
            digest = sha256_file(f)
            checksums.append((digest, str(rel)))

    checksum_file = snapshot_dir / 'checksums.sha256'
    with open(checksum_file, 'w') as f:
        for digest, name in checksums:
            f.write(f"{digest}  {name}\n")
    print(f"  Wrote {len(checksums)} checksums to checksums.sha256")

    # ---- Parse results for manifest ----
    gsa_file = snapshot_dir / 'cell_type' / 'cell_type.gsa.out'
    top_celltypes = parse_gsa_top_hits(gsa_file, top_n=20) if gsa_file.exists() else []

    genes_file = snapshot_dir / 'gene_level' / 'genes.out'
    gene_summary = parse_genes_out_summary(genes_file) if genes_file.exists() else {}

    bonf_celltype = 0.05 / 461
    sig_celltypes = [r for r in top_celltypes if float(r['P']) < bonf_celltype]

    # ---- Build manifest ----
    manifest = {
        'snapshot_id':         snapshot_name,
        'snapshot_created':    datetime.now(timezone.utc).isoformat(),
        'snapshot_host':       f"{platform.system()} {platform.machine()}",
        'git_commit':          get_git_commit(project_root),
        'gwas': {
            'key':             gwas_key,
            'name':            gwas['name'],
            'citation':        gwas['citation'],
            'n_cases':         gwas.get('n_cases'),
            'n_controls':      gwas.get('n_controls'),
            'n_effective':     gwas.get('n_effective'),
            'ancestry':        gwas.get('ancestry'),
            'build':           gwas.get('build'),
        },
        'magma': {
            'required_version': config['magma']['required_version'],
            'window_upstream_kb':   config['magma']['gene_window_upstream_kb'],
            'window_downstream_kb': config['magma']['gene_window_downstream_kb'],
            'model_direction':      config['magma']['model_direction'],
        },
        'qc': {
            'maf_min':            config['qc']['maf_min'],
            'info_min':           config['qc']['info_min'],
            'min_g1000_overlap':  config['qc']['min_g1000_overlap_frac'],
        },
        'gene_level_summary':     gene_summary,
        'celltype_bonferroni':    bonf_celltype,
        'celltype_significant':   {
            'count': len(sig_celltypes),
            'list':  [{'cluster': r['VARIABLE'], 'beta': float(r['BETA']),
                       'p': float(r['P'])} for r in sig_celltypes],
        },
        'celltype_top20':         [
            {'cluster': r['VARIABLE'], 'beta': float(r['BETA']),
             'beta_std': float(r['BETA_STD']), 'p': float(r['P']),
             'ngenes': int(r['NGENES'])}
            for r in top_celltypes
        ],
        'files_copied':           [str(p[1].relative_to(snapshot_dir)) for p in copied],
        'files_missing':          [str(m.relative_to(project_root)) for m in missing],
    }

    manifest_file = snapshot_dir / 'manifest.json'
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"  Wrote manifest.json")

    # ---- Summary ----
    total_size = sum(f.stat().st_size for f in snapshot_dir.rglob('*') if f.is_file())
    print()
    print("=" * 70)
    print(f"✓ Snapshot created: {snapshot_name}")
    print("=" * 70)
    print(f"Location:   {snapshot_dir}")
    print(f"Total size: {total_size / 1e6:.1f} MB")
    print(f"Files:      {len(checksums)}")
    print()
    if sig_celltypes:
        print(f"Bonferroni-significant cell types (P < {bonf_celltype:.2e}):")
        for r in sig_celltypes:
            print(f"  {r['VARIABLE']:12s}  β={float(r['BETA']):.2f}  P={float(r['P']):.2e}")
    else:
        print(f"No cell types passed Bonferroni (threshold {bonf_celltype:.2e})")

    if gene_summary:
        print(f"\nGene-level: {gene_summary['genes_passing_bonferroni']} of "
              f"{gene_summary['total_genes_tested']} genes pass Bonferroni "
              f"({gene_summary['bonferroni_threshold']:.2e})")
    print()
    print("Next steps:")
    print(f"  1. Verify integrity: cd {snapshot_dir.relative_to(project_root)} && shasum -a 256 -c checksums.sha256")
    print(f"  2. Upload to team Google Drive: drag the folder")
    print(f"  3. Commit to git: git add results/snapshots/{snapshot_name}")


if __name__ == '__main__':
    main()
