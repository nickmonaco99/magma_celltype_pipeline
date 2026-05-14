#!/usr/bin/env python3
"""
Build a self-contained team handoff folder for cross-trait cell-type analysis.

Pulls from:
  - results/snapshots/<trait>_v*_*/annotated/celltype_annotated.csv  (our 5 traits)
  - Akshita BMI data/celltype_results_BMI_2018.csv                   (BMI, different schema)
  - Armaan Snapshots/<trait>_finngen_v1_*/annotated/celltype_annotated.csv  (Armaan's 6 FinnGen)
  - results/<TRAIT>/conditional/independent_clusters.txt             (3 traits with conditional)
  - results/<TRAIT>/Siletti_l2_conti-spe_<TRAIT>.gsa.out              (raw MAGMA marginal)

Produces:
  team_handoff_<date>/
    README.md
    per_trait/             one CSV per trait, uniform schema
    cross_trait/           long, wide, and overlap tables
    metadata/              cluster annotations, analysis constants

Run from the project root:
  python pipeline/08_build_team_handoff.py
"""
import argparse
import sys
from pathlib import Path
from datetime import date
import pandas as pd
import numpy as np

# ============================================================================
# Constants
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
N_CLUSTERS = 461
BONF = 0.05 / N_CLUSTERS  # 1.085e-4

# Target uniform schema for per-trait CSVs
CANONICAL_COLS = [
    'trait_key', 'trait_display', 'trait_category',
    'VARIABLE', 'cluster_id_numeric',
    'cluster_name', 'supercluster', 'class_auto', 'neurotransmitter',
    'top_regions', 'top_dissections', 'top_genes',
    'NGENES', 'BETA', 'BETA_STD', 'SE', 'P',
    'log10_neg_P',
    'bonferroni_significant', 'rank_within_trait',
    'independent', 'independent_rank',
]

# Akshita -> canonical column rename map
AKSHIKA_RENAME = {
    'Cluster name': 'cluster_name',
    'Supercluster': 'supercluster',
    'Class auto-annotation': 'class_auto',
    'Top three regions': 'top_regions',
}

# ============================================================================
# Trait registry
# ============================================================================
TRAITS = [
    # 5 traits run through our pipeline (canonical snapshot schema)
    {
        'key': 'scz_pgc3_2022', 'display': 'SCZ',
        'category': 'psychiatric / brain-direct',
        'gwas': 'PGC3 wave3 EUR public v3 (Trubetskoy 2022)',
        'snapshot': 'results/snapshots/scz_pgc3_2022_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': 'results/SCZ_PGC3/Siletti_l2_conti-spe_SCZ_PGC3.gsa.out',
        'independent': 'results/SCZ_PGC3/conditional/independent_clusters.txt',
        'schema': 'canonical',
    },
    {
        'key': 'chronotype_loh2018_ukb', 'display': 'Chronotype',
        'category': 'sleep / brain-direct',
        'gwas': 'UKB chronotype (Loh 2018)',
        'snapshot': 'results/snapshots/chronotype_loh2018_ukb_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': 'results/Chronotype_Loh2018/Siletti_l2_conti-spe_Chronotype_Loh2018.gsa.out',
        'independent': 'results/Chronotype_Loh2018/conditional/independent_clusters.txt',
        'schema': 'canonical',
    },
    {
        'key': 'insomnia_ukb', 'display': 'Insomnia',
        'category': 'sleep / brain-direct',
        'gwas': 'UKB insomnia (UKB-only release)',
        'snapshot': 'results/snapshots/insomnia_ukb_v2_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': 'results/Insomnia_UKB/Siletti_l2_conti-spe_Insomnia_UKB.gsa.out',
        'independent': None,  # only 1 sig cluster, conditional not run
        'schema': 'canonical',
    },
    {
        'key': 't2d_xue2018', 'display': 'T2D',
        'category': 'metabolic / body-mediated',
        'gwas': 'T2D (Xue 2018)',
        'snapshot': 'results/snapshots/t2d_xue2018_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': 'results/T2D_Xue2018/Siletti_l2_conti-spe_T2D_Xue2018.gsa.out',
        'independent': None,  # 0 sig clusters
        'schema': 'canonical',
    },
    {
        'key': 'ibs_eijsbouts2021', 'display': 'IBS',
        'category': 'gastrointestinal / body-mediated',
        'gwas': 'IBS (Eijsbouts 2021)',
        'snapshot': 'results/snapshots/ibs_eijsbouts2021_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': 'results/IBS_Eijsbouts2021/Siletti_l2_conti-spe_IBS_Eijsbouts2021.gsa.out',
        'independent': None,  # 0 sig clusters
        'schema': 'canonical',
    },
    # BMI: Akshita ran the GWAS+gene-level (.genes.raw shared), we ran the marginal cell-type call
    # ourselves with our spec matrix, so BMI is consistent with the rest of our pipeline.
    {
        'key': 'bmi_yengo2018', 'display': 'BMI',
        'category': 'metabolic / brain-direct',
        'gwas': 'BMI (Yengo 2018)',
        'snapshot': 'Akshita BMI data/celltype_results_BMI_2018.csv',
        'gsa_out': 'results/BMI_Yengo2018/Siletti_l2_conti-spe_BMI_Yengo2018.gsa.out',
        'independent': 'results/BMI_Yengo2018/conditional/independent_clusters.txt',
        'schema': 'akshika',
    },
    # Armaan's 6 FinnGen traits (canonical schema)
    {
        'key': 'menorrhagia_finngen', 'display': 'Menorrhagia',
        'category': 'reproductive / body-mediated',
        'gwas': 'FinnGen menorrhagia',
        'snapshot': 'Armaan Snapshots/menorrhagia_finngen_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': None, 'independent': None, 'schema': 'canonical',
    },
    {
        'key': 'menst_excessive_finngen', 'display': 'Excessive menstruation',
        'category': 'reproductive / body-mediated',
        'gwas': 'FinnGen excessive menstruation',
        'snapshot': 'Armaan Snapshots/menst_excessive_finngen_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': None, 'independent': None, 'schema': 'canonical',
    },
    {
        'key': 'menst_irregular_finngen', 'display': 'Irregular menstruation',
        'category': 'reproductive / body-mediated',
        'gwas': 'FinnGen irregular menstruation',
        'snapshot': 'Armaan Snapshots/menst_irregular_finngen_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': None, 'independent': None, 'schema': 'canonical',
    },
    {
        'key': 'migraine_finngen', 'display': 'Migraine',
        'category': 'neurological / brain-direct',
        'gwas': 'FinnGen migraine',
        'snapshot': 'Armaan Snapshots/migraine_finngen_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': None, 'independent': None, 'schema': 'canonical',
    },
    {
        'key': 'sleep_apnea_finngen', 'display': 'Sleep apnea',
        'category': 'sleep / mixed',
        'gwas': 'FinnGen sleep apnea',
        'snapshot': 'Armaan Snapshots/sleep_apnea_finngen_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': None, 'independent': None, 'schema': 'canonical',
    },
    {
        'key': 'sleep_combined_finngen', 'display': 'Sleep combined',
        'category': 'sleep / mixed',
        'gwas': 'FinnGen sleep combined',
        'snapshot': 'Armaan Snapshots/sleep_combined_finngen_v1_2026-05-13/annotated/celltype_annotated.csv',
        'gsa_out': None, 'independent': None, 'schema': 'canonical',
    },
]

# ============================================================================
# Helpers
# ============================================================================
def normalize_to_canonical(df, schema):
    """Rename and reshape a snapshot dataframe to canonical schema."""
    if schema == 'akshika':
        df = df.rename(columns=AKSHIKA_RENAME).copy()
    # Ensure columns we expect to populate exist
    for c in ['cluster_name', 'supercluster', 'class_auto', 'neurotransmitter',
              'top_regions', 'top_dissections', 'top_genes']:
        if c not in df.columns:
            df[c] = pd.NA
    return df

def load_independent_clusters(path):
    """Return ordered list of cluster names (e.g. ['Cluster239', 'Cluster280', ...])."""
    p = PROJECT_ROOT / path
    if not p.exists():
        return None
    with open(p) as f:
        return [line.strip() for line in f if line.strip()]

def cluster_id_from_name(name):
    """'Cluster239' -> 239 (int). Returns None on failure."""
    if not isinstance(name, str) or not name.startswith('Cluster'):
        return None
    try:
        return int(name.replace('Cluster', ''))
    except ValueError:
        return None

def build_per_trait_df(trait_cfg):
    """Load + normalize one trait. Returns a dataframe with CANONICAL_COLS, or None."""
    snap = PROJECT_ROOT / trait_cfg['snapshot']
    if not snap.exists():
        print(f"  [skip] {trait_cfg['key']}: snapshot not found at {snap}")
        return None
    df = pd.read_csv(snap)
    df = normalize_to_canonical(df, trait_cfg['schema'])

    # Trait identifiers
    df['trait_key'] = trait_cfg['key']
    df['trait_display'] = trait_cfg['display']
    df['trait_category'] = trait_cfg['category']

    # Numeric cluster ID for plotting
    df['cluster_id_numeric'] = df['VARIABLE'].apply(cluster_id_from_name)

    # Derived: -log10(P) and Bonferroni
    df['log10_neg_P'] = -np.log10(df['P'].astype(float))
    df['bonferroni_significant'] = df['P'].astype(float) < BONF

    # Rank within trait by P
    df = df.sort_values('P').reset_index(drop=True)
    df['rank_within_trait'] = df.index + 1

    # Independent flag (if conditional analysis was run)
    if trait_cfg.get('independent'):
        order = load_independent_clusters(trait_cfg['independent'])
        if order:
            indep_set = set(order)
            rank_map = {c: i + 1 for i, c in enumerate(order)}
            df['independent'] = df['VARIABLE'].isin(indep_set)
            df['independent_rank'] = df['VARIABLE'].map(rank_map)
        else:
            df['independent'] = False
            df['independent_rank'] = pd.NA
    else:
        df['independent'] = pd.NA  # not analyzed (vs explicitly False)
        df['independent_rank'] = pd.NA

    # Reorder to canonical columns; preserve any extras after
    keep = [c for c in CANONICAL_COLS if c in df.columns]
    extras = [c for c in df.columns if c not in CANONICAL_COLS]
    return df[keep + extras]

# ============================================================================
# Main build
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=str,
                        default=str(PROJECT_ROOT / f'team_handoff_{date.today().isoformat()}'),
                        help='Output directory (default: team_handoff_YYYY-MM-DD/ in project root)')
    args = parser.parse_args()

    out = Path(args.out_dir)
    (out / 'per_trait').mkdir(parents=True, exist_ok=True)
    (out / 'cross_trait').mkdir(parents=True, exist_ok=True)
    (out / 'metadata').mkdir(parents=True, exist_ok=True)

    print(f"Building team handoff at: {out}\n")

    # ---- Per-trait CSVs ----
    print("[1/4] Per-trait CSVs")
    per_trait_dfs = {}
    for i, t in enumerate(TRAITS, 1):
        df = build_per_trait_df(t)
        if df is None:
            continue
        out_path = out / 'per_trait' / f"{i:02d}_{t['key']}.csv"
        df.to_csv(out_path, index=False)
        per_trait_dfs[t['key']] = df
        n_sig = df['bonferroni_significant'].sum()
        n_indep = df['independent'].sum() if df['independent'].dtype != object else \
                  (df['independent'] == True).sum()
        print(f"  [{i:02d}] {t['key']}: {len(df)} clusters, {n_sig} Bonf-sig, "
              f"{n_indep if not pd.isna(n_indep) else 'n/a'} independent")

    if not per_trait_dfs:
        print("ERROR: no traits could be loaded. Check that snapshots exist.")
        sys.exit(1)

    # ---- Cross-trait long format ----
    print("\n[2/4] Cross-trait long-format table")
    long_df = pd.concat(per_trait_dfs.values(), ignore_index=True)
    long_path = out / 'cross_trait' / 'all_traits_long.csv'
    long_df.to_csv(long_path, index=False)
    print(f"  Wrote {len(long_df)} rows across {len(per_trait_dfs)} traits -> {long_path}")

    # ---- Cross-trait wide format (-log10(P) per cluster per trait) ----
    print("\n[3/4] Cross-trait wide-format table")
    wide_pieces = []
    for k, df in per_trait_dfs.items():
        wide_pieces.append(df[['VARIABLE', 'log10_neg_P']].rename(
            columns={'log10_neg_P': f'neglog10P__{k}'}))
    # Merge on VARIABLE
    wide = wide_pieces[0]
    for piece in wide_pieces[1:]:
        wide = wide.merge(piece, on='VARIABLE', how='outer')
    wide['cluster_id_numeric'] = wide['VARIABLE'].apply(cluster_id_from_name)
    wide = wide.sort_values('cluster_id_numeric').reset_index(drop=True)
    cols = ['VARIABLE', 'cluster_id_numeric'] + [c for c in wide.columns
                                                  if c.startswith('neglog10P__')]
    wide_path = out / 'cross_trait' / 'all_traits_wide_neglog10P.csv'
    wide[cols].to_csv(wide_path, index=False)
    print(f"  Wrote {len(wide)} clusters x {len(per_trait_dfs)} traits -> {wide_path}")

    # ---- Independent cross-trait overlap ----
    print("\n[4/4] Independent cluster overlap")
    indep_sets = {}
    for t in TRAITS:
        if t.get('independent'):
            order = load_independent_clusters(t['independent'])
            if order:
                indep_sets[t['key']] = order

    if indep_sets:
        # Build per-trait independent CSV with biology
        per_trait_indep_rows = []
        for trait_key, clusters in indep_sets.items():
            tdf = per_trait_dfs[trait_key]
            for idx, c in enumerate(clusters, 1):
                row = tdf[tdf['VARIABLE'] == c].iloc[0]
                per_trait_indep_rows.append({
                    'trait_key': trait_key,
                    'trait_display': row['trait_display'],
                    'independent_rank': idx,
                    'VARIABLE': c,
                    'cluster_id_numeric': cluster_id_from_name(c),
                    'cluster_name': row.get('cluster_name'),
                    'supercluster': row.get('supercluster'),
                    'class_auto': row.get('class_auto'),
                    'top_regions': row.get('top_regions'),
                    'P': row['P'],
                    'log10_neg_P': row['log10_neg_P'],
                })
        indep_df = pd.DataFrame(per_trait_indep_rows)
        indep_path = out / 'cross_trait' / 'independent_clusters_per_trait.csv'
        indep_df.to_csv(indep_path, index=False)
        print(f"  Wrote {len(indep_df)} (trait, cluster) rows -> {indep_path}")

        # Cross-trait overlap: clusters appearing as independent in 2+ traits
        from collections import defaultdict
        cluster_to_traits = defaultdict(list)
        for trait_key, clusters in indep_sets.items():
            for c in clusters:
                cluster_to_traits[c].append(trait_key)

        overlap_rows = []
        for c, traits in cluster_to_traits.items():
            if len(traits) >= 2:
                # Get biology from any trait's snapshot
                any_trait = traits[0]
                tdf = per_trait_dfs[any_trait]
                row = tdf[tdf['VARIABLE'] == c].iloc[0]
                overlap_rows.append({
                    'VARIABLE': c,
                    'cluster_id_numeric': cluster_id_from_name(c),
                    'n_traits_independent': len(traits),
                    'traits_independent': ', '.join(sorted(traits)),
                    'cluster_name': row.get('cluster_name'),
                    'supercluster': row.get('supercluster'),
                    'class_auto': row.get('class_auto'),
                    'top_regions': row.get('top_regions'),
                })
        if overlap_rows:
            overlap_df = pd.DataFrame(overlap_rows).sort_values(
                ['n_traits_independent', 'cluster_id_numeric'],
                ascending=[False, True])
            overlap_path = out / 'cross_trait' / 'independent_cross_trait_overlap.csv'
            overlap_df.to_csv(overlap_path, index=False)
            print(f"  Wrote {len(overlap_df)} multi-trait independent clusters -> {overlap_path}")

    # ---- Cluster annotations (master file) ----
    print("\n[bonus] Master cluster annotation table")
    # Use any trait's snapshot as the annotation source (they're identical for the 461 clusters)
    annot_source_key = next(iter(per_trait_dfs))
    annot = per_trait_dfs[annot_source_key][[
        'VARIABLE', 'cluster_id_numeric', 'cluster_name', 'supercluster',
        'class_auto', 'neurotransmitter', 'top_regions', 'top_dissections', 'top_genes'
    ]].drop_duplicates('VARIABLE').sort_values('cluster_id_numeric')
    annot_path = out / 'metadata' / 'siletti_cluster_annotations.csv'
    annot.to_csv(annot_path, index=False)
    print(f"  Wrote {len(annot)} cluster annotations -> {annot_path}")

    # ---- Analysis constants ----
    consts_path = out / 'metadata' / 'analysis_constants.txt'
    consts_path.write_text(f"""# Analysis constants for cell-type MAGMA pipeline
# Generated: {date.today().isoformat()}

N_CLUSTERS = {N_CLUSTERS}
BONFERRONI_THRESHOLD = {BONF}
BONFERRONI_THRESHOLD_FORMULA = 0.05 / {N_CLUSTERS}
MAGMA_VERSION = 1.10
MAGMA_MODEL_MARGINAL = direction=greater
MAGMA_MODEL_CONDITIONAL = joint-pairs direction=greater
SPECIFICITY_MATRIX = Siletti_l2_conti_specificity_matrix.txt (continuous)
ATLAS = Siletti 2023 adult human brain, level 2 (461 clusters)
LD_REFERENCE = 1000 Genomes Phase 3 EUR
GENE_ANNOTATION = NCBI 37.3
GENE_WINDOW = 35kb upstream, 10kb downstream
CONDITIONAL_ALGORITHM = Duncan 2025 forward selection (Nat Neurosci)
""")
    print(f"  Wrote analysis constants -> {consts_path}")

    # ---- README ----
    write_readme(out, per_trait_dfs, indep_sets if indep_sets else {})

    print(f"\nDone. Handoff folder: {out}")

def write_readme(out, per_trait_dfs, indep_sets):
    """Write the team-facing README."""
    n_traits = len(per_trait_dfs)
    n_conditional = len(indep_sets)
    n_clusters_total = N_CLUSTERS

    # Build per-trait table
    rows = []
    for t in TRAITS:
        if t['key'] not in per_trait_dfs:
            continue
        df = per_trait_dfs[t['key']]
        n_sig = int(df['bonferroni_significant'].sum())
        n_indep = (df['independent'] == True).sum() if df['independent'].dtype != object else 'n/a'
        rows.append((t['display'], t['category'], t['gwas'], n_sig, n_indep))

    trait_table = "| Trait | Category | GWAS source | Bonferroni-sig | Independent |\n"
    trait_table += "|---|---|---|---|---|\n"
    for r in rows:
        trait_table += f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} |\n"

    readme = f"""# Cross-trait cell-type analysis: team handoff

Generated: {date.today().isoformat()}
Pipeline: MAGMA cell-type analysis + Duncan 2025 conditional forward selection
Atlas: Siletti 2023 adult human brain, level 2 ({n_clusters_total} clusters)

## Contents

```
per_trait/                      One CSV per trait, uniform schema.
  01_scz_pgc3_2022.csv          Schizophrenia (PGC3)
  02_chronotype_loh2018_ukb.csv Chronotype (UKB)
  03_insomnia_ukb.csv           Insomnia (UKB)
  04_t2d_xue2018.csv            T2D (Xue 2018)
  05_ibs_eijsbouts2021.csv      IBS (Eijsbouts 2021)
  06_bmi_yengo2018.csv          BMI (Yengo 2018) — Akshita's gene-level + our cell-type
  07-12_*_finngen.csv           Armaan's 6 FinnGen traits

cross_trait/
  all_traits_long.csv                  Long format: one row per (trait, cluster)
  all_traits_wide_neglog10P.csv        Wide format: rows=clusters, cols=traits
  independent_clusters_per_trait.csv   Forward-selection independent clusters with biology
  independent_cross_trait_overlap.csv  Clusters independent in 2+ traits (the key finding)

metadata/
  siletti_cluster_annotations.csv  All 461 clusters + supercluster/region/dissection
  analysis_constants.txt           Bonferroni threshold, MAGMA params, etc.
```

## Per-trait summary

{trait_table}

## Headline finding

Cluster132 (retrosplenial cortex, A29-A30 excitatory neuron) is forward-selection
independent across SCZ, Chronotype, and BMI. Cluster228 (amygdala eccentric MSN)
and Cluster428 (basal forebrain septum splatter) are independent in both SCZ and
BMI. See `cross_trait/independent_cross_trait_overlap.csv`.

## Schema for per-trait CSVs

| Column | Description |
|---|---|
| trait_key, trait_display, trait_category | Trait identifiers |
| VARIABLE | Cluster name as string, e.g. "Cluster239" |
| cluster_id_numeric | Same as integer, e.g. 239 (for plotting x-axis) |
| cluster_name, supercluster, class_auto, neurotransmitter | Biological annotation from Siletti Table S3 |
| top_regions, top_dissections, top_genes | Annotation strings (% breakdowns) |
| NGENES, BETA, BETA_STD, SE, P | MAGMA cell-type test outputs |
| log10_neg_P | -log10(P), pre-computed for plotting |
| bonferroni_significant | bool, P < {BONF:.3e} (= 0.05/461) |
| rank_within_trait | Integer rank by P ascending (1 = top hit) |
| independent | bool / NaN if conditional analysis not run for this trait |
| independent_rank | 1..N for clusters in the forward-selection list |

## Quick plotting starter code

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('per_trait/01_scz_pgc3_2022.csv')

# Cluster-id scatter coloured by supercluster, Duncan Fig 3 style
fig, ax = plt.subplots(figsize=(12, 4))
for sup, sdf in df.groupby('supercluster'):
    ax.scatter(sdf['cluster_id_numeric'], sdf['log10_neg_P'],
               s=20, label=sup, alpha=0.7)
ax.axhline(-np.log10({BONF}), ls='--', c='grey', label='Bonferroni')
indep = df[df['independent'] == True]
ax.scatter(indep['cluster_id_numeric'], indep['log10_neg_P'],
           s=100, facecolors='none', edgecolors='black', linewidths=1.5,
           label='Independent', zorder=10)
ax.set_xlabel('Cluster ID')
ax.set_ylabel('-log10(P)')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig('scz_celltype.png', dpi=150)
```

## Methods

1. **Gene-level analysis** (MAGMA Steps 1-2)
   - Per-SNP P values mapped to genes within 35kb upstream / 10kb downstream window
   - LD reference: 1000 Genomes Phase 3 EUR
   - Output: per-gene Z-statistics

2. **Marginal cell-type analysis** (MAGMA gene-property)
   - 461 Siletti l2 clusters tested individually
   - Continuous specificity matrix (Siletti_l2_conti_specificity_matrix.txt)
   - Flags: `--model direction=greater` (one-sided, positive set + covar)
   - Bonferroni threshold: 0.05 / 461 = {BONF:.3e}

3. **Conditional analysis** (Duncan 2025 Stages 4-5)
   - Spec matrix filtered to Bonferroni-significant clusters (per trait)
   - MAGMA `--model joint-pairs direction=greater` runs all pairwise conditional tests
   - Forward selection (Duncan/Tayden Li R algorithm) identifies independent clusters
     using PS = log10(P_cond) / log10(P_marg) with thresholds 0.5 and 0.8
   - Only run for traits with ≥2 Bonferroni-significant marginal clusters

## Caveats

1. **Gene-set drift from Duncan's reference**: our marginal analysis uses 16,333
   genes; Duncan's reference SCZ analysis used 16,641 (~1.8% gene loss).
   Recovers same top hits and biology, but a few middle-rank clusters flip near
   the Bonferroni boundary.

2. **Sample size for conditional analysis**: 3 of 12 traits have enough sig
   clusters for forward selection (SCZ, Chronotype, BMI). T2D (0 sig), IBS (0),
   Insomnia (1) are marginal-only.

3. **BMI provenance**: gene-level analysis by Akshita (Yengo 2018 sumstats);
   cell-type analysis re-run by Nick to use the same specificity matrix as
   the other 5 traits for cross-trait consistency.

4. **FinnGen traits**: marginal-only from Armaan; no conditional analysis yet.

5. **Cluster numbers are Siletti l2 IDs**, not arbitrary. Cluster239 in our
   files == Cluster239 in Duncan's published analysis (same MGE interneuron).

## Contact

Nick Monaco, UCL CCMN (Brierley lab)
"""
    (out / 'README.md').write_text(readme)
    print(f"\nWrote README -> {out / 'README.md'}")

if __name__ == '__main__':
    main()
