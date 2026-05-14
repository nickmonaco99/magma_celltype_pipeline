"""
Driver-gene profiling for BMI independent clusters.

For each of the 12 BMI independent clusters, computes:
  - Top 5 markers by specificity (defines what the cluster IS)
  - Top 10 BMI drivers by ZSTAT * specificity (drives the BMI signal)
  - Top canonical feeding/obesity gene present
  - Counts of developmental TFs vs feeding genes in top 10 drivers
  - Heuristic confidence tier (1/2/3)

Outputs:
  - Compact summary table to console
  - results/BMI_Yengo2018/conditional/cluster_driver_profiles_summary.csv
  - results/BMI_Yengo2018/conditional/cluster_driver_profiles_detail/<cluster>.csv
"""
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---- Load reference files ----
print("Loading reference files...")
spec = pd.read_csv(PROJECT_ROOT / 'gene-level/Siletti_l2_conti_specificity_matrix.txt', sep=r'\s+')
bmi = pd.read_csv(PROJECT_ROOT / 'Akshita BMI data/bmi_yengo2018.no_heading.step2.genes.out',
                  sep=r'\s+', comment='#')
loc = pd.read_csv(PROJECT_ROOT / 'auxfiles/NCBI37.3.gene.loc', sep=r'\s+', header=None,
                  names=['entrez','chr','start','stop','strand','symbol'])
ak = pd.read_csv(PROJECT_ROOT / 'Akshita BMI data/celltype_results_BMI_2018.csv')
s3 = pd.read_excel(PROJECT_ROOT / 'auxfiles/science.add7046_table_s3.xlsx', sheet_name=0)

# ---- BMI independent clusters in algorithm order (from independent_clusters.txt) ----
bmi_indep = ['Cluster251','Cluster228','Cluster286','Cluster132','Cluster389',
             'Cluster428','Cluster386','Cluster305','Cluster410','Cluster439',
             'Cluster332','Cluster311']

# ---- Canonical feeding / obesity genes ----
# Hypothalamic neuropeptides, feeding receptors, monogenic obesity loci, big BMI GWAS hits
FEEDING_GENES = {
    # Neuropeptides
    'POMC','AGRP','NPY','CARTPT','PMCH','TRH','CRH','OXT','AVP','HCRT',
    'GHRH','QRFP','GAL','VGF','ADCYAP1','PYY','NTS','GIP','GLP1','GHRL',
    # Receptors
    'MC4R','MC3R','MC2R','LEPR','GHSR','HCRTR1','HCRTR2','GLP1R','CCKAR',
    'CCKBR','NPY1R','NPY5R','GALR1','GALR2','MRAP','GPR75',
    # Processing / signaling
    'PCSK1','LEP','BDNF','INS','UCP1','SH2B1',
    # Hypothalamic TFs (definitionally feeding)
    'SIM1','FEZF1','NKX2-1','OTP','NR5A1',
    # Major BMI GWAS hits
    'FTO','TFAP2B','SKOR1','KLF14','FAIM2','TMEM18',
}

# ---- Developmental positional TFs (low functional interpretation) ----
def is_developmental(symbol):
    if not isinstance(symbol, str):
        return False
    return (symbol.startswith('HOX') or
            symbol.startswith('POU4F') or
            symbol.startswith('HMX') or
            symbol.startswith('BARX') or
            (symbol.startswith('EN') and len(symbol) == 3 and symbol[2:].isdigit()))

# ---- Output paths ----
out_dir = PROJECT_ROOT / 'results/BMI_Yengo2018/conditional'
detail_dir = out_dir / 'cluster_driver_profiles_detail'
detail_dir.mkdir(parents=True, exist_ok=True)

# ---- Loop ----
print(f"Profiling {len(bmi_indep)} BMI independent clusters...\n")
rows = []

for cluster in bmi_indep:
    cid = int(cluster.replace('Cluster',''))
    s3_row = s3[s3['Cluster ID'] == cid].iloc[0]
    ak_row = ak[ak['VARIABLE'] == cluster].iloc[0]

    # Combine specificity, gene symbols, BMI Z-scores
    m = spec[['GENE', cluster]].rename(columns={cluster: 'specificity'})
    m = m.merge(loc[['entrez','symbol']], left_on='GENE', right_on='entrez', how='left')
    m = m.merge(bmi[['GENE','ZSTAT','P']], on='GENE', how='left')
    m['score'] = m['ZSTAT'] * m['specificity']

    # Top markers and drivers
    top_markers = m.nlargest(25, 'specificity')[['symbol','specificity','ZSTAT','P']]
    top_drivers = m.nlargest(25, 'score')[['symbol','specificity','ZSTAT','P','score']]

    # Save detailed per-cluster CSVs
    top_markers.to_csv(detail_dir / f'{cluster}_markers.csv', index=False)
    top_drivers.to_csv(detail_dir / f'{cluster}_bmi_drivers.csv', index=False)

    # Heuristic tier classification (looking at top 10 drivers)
    top10 = top_drivers.head(10)
    n_dev = sum(top10['symbol'].apply(is_developmental))
    n_feeding = sum(top10['symbol'].isin(FEEDING_GENES))

    if n_feeding >= 3:
        tier, tier_label = 1, 'Functional'
    elif n_feeding >= 1 and n_dev <= 1:
        tier, tier_label = 2, 'Mixed'
    elif n_dev >= 3:
        tier, tier_label = 3, 'Developmental'
    else:
        tier, tier_label = 2, 'Mixed'

    # Best feeding gene present
    feeding_hits = m[m['symbol'].isin(FEEDING_GENES) & m['ZSTAT'].notna()].sort_values('score', ascending=False)
    top_feeding = feeding_hits.iloc[0] if len(feeding_hits) > 0 else None

    rows.append({
        'cluster': cluster,
        'cluster_id': cid,
        'supercluster': s3_row['Supercluster'],
        'top_region': str(s3_row['Top three regions']).split(',')[0].strip(),
        'top_dissection': str(s3_row['Top three dissections']).split(',')[0].replace('Human ','').strip(),
        'P_bmi': float(ak_row['P']),
        'top_3_markers': ', '.join(top_markers.head(3)['symbol'].astype(str).tolist()),
        'top_5_drivers': ', '.join(top_drivers.head(5)['symbol'].astype(str).tolist()),
        'top_feeding_gene': top_feeding['symbol'] if top_feeding is not None else None,
        'top_feeding_spec': round(top_feeding['specificity'], 4) if top_feeding is not None else None,
        'top_feeding_Z': round(top_feeding['ZSTAT'], 2) if top_feeding is not None else None,
        'n_dev_TFs_in_top10': n_dev,
        'n_feeding_genes_in_top10': n_feeding,
        'tier': tier,
        'tier_label': tier_label,
    })

summary = pd.DataFrame(rows)

# ---- Save summary CSV ----
summary_path = out_dir / 'cluster_driver_profiles_summary.csv'
summary.to_csv(summary_path, index=False)
print(f"Saved summary -> {summary_path}")
print(f"Saved per-cluster details -> {detail_dir}/")
print()

# ---- Pretty-print to console ----
print("="*120)
print("BMI INDEPENDENT CLUSTER DRIVER PROFILES")
print("="*120)

# Group by tier
for tier_num, tier_name in [(1,'Functional (feeding genes dominate top drivers)'),
                            (2,'Mixed (some feeding genes, some non-functional drivers)'),
                            (3,'Developmental / positional (HOX/POU TFs dominate)')]:
    sub = summary[summary['tier'] == tier_num].sort_values('P_bmi')
    if len(sub) == 0:
        continue
    print(f"\n--- TIER {tier_num}: {tier_name} ({len(sub)} clusters) ---\n")
    for _, r in sub.iterrows():
        print(f"  {r['cluster']:>12s}  P={r['P_bmi']:.2e}  {r['supercluster']:<30s}  "
              f"top region: {r['top_region']}")
        print(f"               top markers:  {r['top_3_markers']}")
        print(f"               top BMI drivers: {r['top_5_drivers']}")
        if r['top_feeding_gene']:
            print(f"               top feeding gene: {r['top_feeding_gene']} "
                  f"(spec={r['top_feeding_spec']}, Z={r['top_feeding_Z']})")
        print(f"               dev TFs in top 10: {r['n_dev_TFs_in_top10']},  "
              f"feeding genes in top 10: {r['n_feeding_genes_in_top10']}")
        print()
