# %% [markdown]
# # 02_qc_gene_results — Validate MAGMA gene-level output
#
# **Runs AFTER**: `02_run_magma_step1and2.sh` produces `*.step2.genes.out`
# **Runs BEFORE**: `03_run_magma_celltype.sh` (cell-type analysis)
#
# **Purpose**: Formal Bonferroni assertion on sentinel genes. If insomnia signal
# isn't where it should be at the gene level, we don't waste 5 min on cell-type
# analysis (and more importantly, the cell-type result wouldn't be trustworthy).
#
# **Pass criteria** (anchored to Jansen 2019 UKB-only, N~386k, comparable to ours):
#   1. ALL 5 sentinel genes have gene-level P < 2.5e-6 (Bonferroni: 0.05/20000)
#   2. At least 3 of 5 have gene-level P < 1e-8
#
# **Also prints**: top 30 gene-level hits for eyeball comparison against
# Jansen 2019 Supplementary Table 8.

# %%
# =============================================================================
# Imports & config
# =============================================================================
import sys
import yaml
import pathlib
from datetime import datetime, timezone

import numpy as np
import pandas as pd

pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 200)
pd.set_option('display.max_rows', 50)

with open('00_config.yaml') as f:
    config = yaml.safe_load(f)

project_root = pathlib.Path(config['project_root'])
gwas_key     = config['active_gwas']
gwas         = config['gwas'][gwas_key]
results_dir  = project_root / config['paths']['results_dir'] / gwas['output_folder']
provenance_log = results_dir / 'run.log'

print(f"GWAS:        {gwas['name']}")
print(f"Results dir: {results_dir}")


# %%
# =============================================================================
# 1. Locate the .step2.genes.out file
# =============================================================================
# Duncan's naming: <sumstats_filename>.step2.genes.out
# Our sumstats is named insomnia_ukb.no_heading → output is
# insomnia_ukb.no_heading.step2.genes.out

genes_out_path = results_dir / "insomnia_ukb.no_heading.step2.genes.out"
genes_log_path = results_dir / "insomnia_ukb.no_heading.step2.log"

assert genes_out_path.exists(), \
    f"\nGene-level output not found: {genes_out_path}\n" \
    f"Did you run 02_run_magma_step1and2.sh successfully?"

print(f"Reading: {genes_out_path}")
print(f"  Size: {genes_out_path.stat().st_size / 1e6:.1f} MB")


# %%
# =============================================================================
# 2. Load gene-level results
# =============================================================================
# Format (from Duncan's Example_results):
#   GENE     CHR    START    STOP   NSNPS  NPARAM   N      ZSTAT      P
#   148398   1      824993   889961 188    25       51186  1.4604     0.072086
# GENE is Entrez ID. Whitespace-separated, single header row.

genes = pd.read_csv(genes_out_path, sep=r'\s+', dtype={'GENE': str, 'CHR': str})
print(f"Loaded {len(genes):,} genes")
print(f"\nColumns: {list(genes.columns)}")
print(f"\n--- First 5 rows ---")
print(genes.head())
print(f"\n--- p-value distribution ---")
print(f"  min:    {genes['P'].min():.2e}")
print(f"  median: {genes['P'].median():.4f}")
print(f"  max:    {genes['P'].max():.4f}")


# %%
# =============================================================================
# 3. Attach gene symbols (for human-readable output)
# =============================================================================
gene_loc = pd.read_csv(
    project_root / config['paths']['gene_loc'],
    sep='\t', header=None,
    names=['ENTREZ', 'CHR_loc', 'START_loc', 'STOP_loc', 'STRAND', 'SYMBOL'],
    dtype={'ENTREZ': str, 'CHR_loc': str}
)
print(f"Loaded {len(gene_loc):,} gene locations from NCBI37.3.gene.loc")

genes = genes.merge(
    gene_loc[['ENTREZ', 'SYMBOL']].rename(columns={'ENTREZ': 'GENE'}),
    on='GENE', how='left'
)
n_with_symbol = genes['SYMBOL'].notna().sum()
print(f"  {n_with_symbol:,} / {len(genes):,} genes have a symbol attached")


# %%
# =============================================================================
# 4. Top 30 gene-level hits — for eyeball comparison vs Jansen 2019 Supp Table 8
# =============================================================================
# Jansen 2019 Supplementary Table 8 lists the top gene-level MAGMA hits from
# their UKB-only analysis. We expect overlap on most of the top genes — exact
# rank order will differ slightly because of MAGMA version (v1.07 vs v1.10),
# gene window (2,1 vs 35,10), and gene.loc file updates.

top30 = genes.nsmallest(30, 'P')[['GENE', 'SYMBOL', 'CHR', 'START', 'NSNPS', 'ZSTAT', 'P']]
top30 = top30.reset_index(drop=True)
top30.index = top30.index + 1
top30.index.name = 'rank'

print("\n" + "=" * 80)
print(f"TOP 30 GENE-LEVEL HITS — {gwas['name']}")
print("=" * 80)
print(top30.to_string(formatters={'P': lambda x: f"{x:.2e}", 'ZSTAT': lambda x: f"{x:.2f}"}))

print("\nCompare against Jansen 2019 Supp Table 8 (UKB-only top genes).")
print("Expect MEIS1, BTBD9, PTPRD, AUTS2, LSAMP, MEF2C, NEGR1 to appear in the top tier.")


# %%
# =============================================================================
# 5. FORMAL SENTINEL ASSERTION — the gate to cell-type analysis
# =============================================================================
# Pass rule (from config, anchored to Jansen 2019 UKB-only):
#   - ALL 5 sentinels at gene-level P < 2.5e-6
#   - At least 3 of 5 at gene-level P < 1e-8

print("\n" + "=" * 80)
print("FORMAL SENTINEL ASSERTION")
print("=" * 80)
print(f"{'Symbol':<10}{'Entrez':<10}{'Rank':<8}{'P':<15}{'≤2.5e-6?':<12}{'≤1e-8?':<10}{'Note'}")
print("-" * 80)

# Build the rank lookup once
genes_sorted = genes.sort_values('P').reset_index(drop=True)
genes_sorted['rank'] = range(1, len(genes_sorted) + 1)
rank_lookup = dict(zip(genes_sorted['GENE'], genes_sorted['rank']))
p_lookup    = dict(zip(genes_sorted['GENE'], genes_sorted['P']))

results = []
for sentinel in gwas['sentinel_genes']:
    entrez = sentinel['entrez']
    symbol = sentinel['symbol']
    if entrez not in p_lookup:
        print(f"{symbol:<10}{entrez:<10}{'-':<8}{'MISSING':<15}{'FAIL':<12}{'FAIL':<10}{sentinel['note']}")
        results.append({'symbol': symbol, 'p': np.nan,
                        'bonf_pass': False, 'strict_pass': False})
        continue

    p = p_lookup[entrez]
    rank = rank_lookup[entrez]
    bonf_pass   = p < sentinel['gene_level_p_max_bonf']     # 2.5e-6
    strict_pass = p < sentinel['gene_level_p_max_strict']    # 1e-8

    bonf_str   = '✓ PASS' if bonf_pass else '✗ FAIL'
    strict_str = '✓ PASS' if strict_pass else '✗ FAIL'

    print(f"{symbol:<10}{entrez:<10}{rank:<8}{p:<15.2e}{bonf_str:<12}{strict_str:<10}{sentinel['note']}")
    results.append({'symbol': symbol, 'p': p,
                    'bonf_pass': bonf_pass, 'strict_pass': strict_pass})

print("-" * 80)

n_bonf   = sum(r['bonf_pass'] for r in results)
n_strict = sum(r['strict_pass'] for r in results)
print(f"\nSummary:")
print(f"  Bonferroni (P < 2.5e-6):    {n_bonf}/5  (need: 5)")
print(f"  Strict     (P < 1e-8):      {n_strict}/5  (need: ≥ 3)")


# %%
# =============================================================================
# 6. THE ASSERTION
# =============================================================================
# This is the formal gate. AssertionError here means do NOT proceed.

bonf_pass_all   = (n_bonf == 5)
strict_pass_min = (n_strict >= 3)

# Provenance logging
ts = datetime.now(timezone.utc).isoformat()
with open(provenance_log, 'a') as f:
    f.write(f"[{ts}] 02_qc_gene_results: bonf={n_bonf}/5, strict={n_strict}/5\n")
    for r in results:
        f.write(f"[{ts}]   {r['symbol']}: P={r['p']:.2e} "
                f"bonf={r['bonf_pass']} strict={r['strict_pass']}\n")

if not (bonf_pass_all and strict_pass_min):
    failure_msg = f"""

*** GENE-LEVEL SENTINEL ASSERTION FAILED ***

Pass criteria:
  - ALL 5 sentinels at P < 2.5e-6: {n_bonf}/5
  - ≥3 of 5 at P < 1e-8:           {n_strict}/5

Do NOT proceed to cell-type analysis. Diagnose first.

Likely causes (ranked by probability):
  1. Wrong gene.loc file (build mismatch). Confirm NCBI37.3.gene.loc (hg19).
  2. Wrong reference panel ancestry. Confirm g1000_eur.bim. EUR for European GWAS.
  3. Wrong column in --pval (e.g. picked an effect-size column instead of P).
     Check log: results/INS_UKB/insomnia_ukb.no_heading.step2.log
  4. NEFF column got corrupted in prep. Check first 5 rows of insomnia_ukb.no_heading.
  5. Insomnia signal at these specific genes is just genuinely weaker in
     Watanabe UKB-only than in Jansen 2019 (different cohort definitions).
     If 4/5 pass Bonferroni and the single failure is borderline (P ≈ 5e-6),
     accept and document. If multiple sentinels fail by orders of magnitude,
     there's a real bug.
"""
    print(failure_msg)
    raise AssertionError(failure_msg.strip())

print("\n" + "=" * 80)
print(f"✓ SENTINEL ASSERTION PASSED")
print(f"  Bonferroni: 5/5 below 2.5e-6")
print(f"  Strict:     {n_strict}/5 below 1e-8 (need ≥3)")
print(f"  Safe to proceed to cell-type analysis (03_run_magma_celltype.sh)")
print("=" * 80)


# %%
# =============================================================================
# 7. Save the top-30 table for the methods write-up
# =============================================================================
top30_out = results_dir / "top30_gene_level_hits.tsv"
top30.to_csv(top30_out, sep='\t')
print(f"\nSaved top-30 table: {top30_out}")
print(f"Compare manually against Jansen 2019 Supp Table 8 (UKB-only column).")
