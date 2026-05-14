# Cross-trait cell-type analysis: team handoff

Generated: 2026-05-13
Pipeline: MAGMA cell-type analysis + Duncan 2025 conditional forward selection
Atlas: Siletti 2023 adult human brain, level 2 (461 clusters)

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

| Trait | Category | GWAS source | Bonferroni-sig | Independent |
|---|---|---|---|---|
| SCZ | psychiatric / brain-direct | PGC3 wave3 EUR public v3 (Trubetskoy 2022) | 105 | 9 |
| Chronotype | sleep / brain-direct | UKB chronotype (Loh 2018) | 48 | 11 |
| Insomnia | sleep / brain-direct | UKB insomnia (UKB-only release) | 1 | n/a |
| T2D | metabolic / body-mediated | T2D (Xue 2018) | 0 | n/a |
| IBS | gastrointestinal / body-mediated | IBS (Eijsbouts 2021) | 0 | n/a |
| BMI | metabolic / brain-direct | BMI (Yengo 2018) | 81 | 12 |
| Menorrhagia | reproductive / body-mediated | FinnGen menorrhagia | 0 | n/a |
| Excessive menstruation | reproductive / body-mediated | FinnGen excessive menstruation | 0 | n/a |
| Irregular menstruation | reproductive / body-mediated | FinnGen irregular menstruation | 0 | n/a |
| Migraine | neurological / brain-direct | FinnGen migraine | 0 | n/a |
| Sleep apnea | sleep / mixed | FinnGen sleep apnea | 0 | n/a |
| Sleep combined | sleep / mixed | FinnGen sleep combined | 0 | n/a |


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
| bonferroni_significant | bool, P < 1.085e-04 (= 0.05/461) |
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
ax.axhline(-np.log10(0.00010845986984815619), ls='--', c='grey', label='Bonferroni')
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
   - Bonferroni threshold: 0.05 / 461 = 1.085e-04

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
