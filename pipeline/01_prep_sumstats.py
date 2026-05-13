# %% [markdown]
# # 01: Prepare Watanabe Insomnia UKB Sumstats for MAGMA
#
# **Goal**: Convert raw Watanabe et al. 2022 insomnia UKB-only sumstats into
# the format Duncan's MAGMA scripts expect, with QC at every step.
#
# **What this notebook does**:
# 1. Load config, verify MAGMA binary version matches Duncan's v1.10
# 2. Inspect raw sumstats file structure (run-then-update-config-then-rerun)
# 3. Apply column mapping, compute NEFF, filter, write Duncan-format outputs
# 4. QC: λ_GC, SNP-level proxy check on Jansen-2019 sentinel genes
# 5. **Pre-flight overlap gate**: assert ≥70% of our SNPs are in g1000_eur.bim
#    (catches build mismatches before 45 min of MAGMA compute)
#
# **Two-phase first run**:
# - Phase 1: cells 1-5 (inspect file). Then update `columns:` in 00_config.yaml.
# - Phase 2: re-run from cell 1 through end.

# %%
# =============================================================================
# Imports & setup
# =============================================================================
import os
import sys
import gzip
import yaml
import json
import shlex
import pathlib
import platform
import subprocess
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2

pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 200)

print(f"Python:    {sys.version.split()[0]}")
print(f"pandas:    {pd.__version__}")
print(f"numpy:     {np.__version__}")
print(f"Platform:  {platform.system()} {platform.machine()}")
print(f"Timestamp: {datetime.now(timezone.utc).isoformat()}")


# %%
# =============================================================================
# 1. Load configuration
# =============================================================================
with open('00_config.yaml') as f:
    config = yaml.safe_load(f)

project_root = pathlib.Path(config['project_root'])
assert project_root.exists(), f"project_root does not exist: {project_root}"

gwas_key = config['active_gwas']
gwas = config['gwas'][gwas_key]

print(f"Active GWAS: {gwas['name']}")
print(f"  N_cases:    {gwas['n_cases']:,}")
print(f"  N_controls: {gwas['n_controls']:,}")
print(f"  N_eff:      {gwas['n_effective']:,}")
print(f"  Build:      {gwas['build']}")

# Resolve paths
gwas_raw_path        = project_root / config['paths']['gwas_dir'] / gwas['raw_subfolder'] / gwas['raw_filename']
results_dir          = project_root / config['paths']['results_dir'] / gwas['output_folder']
results_dir.mkdir(parents=True, exist_ok=True)
file_prefix = config['gwas'][gwas_key].get('file_prefix', gwas_key)
output_sumstats_path = results_dir / f"{file_prefix}.no_heading"
output_snploc_path   = results_dir / f"snploc_{file_prefix}"
provenance_log       = results_dir / "run.log"


# %%
# =============================================================================
# 2. MAGMA version pin — gracefully verify v1.10
# =============================================================================
# Strategy: try to execute the binary. If it runs (Linux host or inside Docker),
# verify the version. If it FAILS to execute (macOS host with Linux binary),
# fall back to trusting the Dockerfile build-time version check.
#
# Why this works: when `docker build -t magma_celltype:latest .` succeeds, the
# Dockerfile's own RUN step (`magma --version 2>&1 | grep -q "v1.10"`) has
# already verified the binary is v1.10. So if you built the image, the version
# is guaranteed correct.
#
# The check only needs to happen at MAGMA invocation time, which is always
# inside Docker (via run_in_docker.sh). The host doesn't need a working MAGMA.

magma_binary_path = project_root / config['paths']['magma_binary']
required_version  = config['magma']['required_version']

# Pick path: prefer in-container path if it exists, otherwise host path
if pathlib.Path('/usr/local/bin/magma').exists():
    magma_cmd = '/usr/local/bin/magma'
    runtime_env = 'docker'
elif magma_binary_path.exists():
    magma_cmd = str(magma_binary_path)
    runtime_env = 'host'
else:
    raise FileNotFoundError(
        f"\nMAGMA binary not found at either location:\n"
        f"  Container: /usr/local/bin/magma\n"
        f"  Host:      {magma_binary_path}\n"
        f"Download MAGMA v1.10 Linux static binary to magma_binary/magma."
    )

print(f"Runtime environment: {runtime_env}")
print(f"MAGMA binary path:   {magma_cmd}")

# Try to run --version. If the binary can't execute on this host (Linux binary
# on macOS), gracefully degrade.
magma_version_output = None
try:
    result = subprocess.run(
        [magma_cmd, '--version'],
        capture_output=True, text=True, timeout=5
    )
    # MAGMA prints version to stderr, not stdout
    magma_version_output = (result.stdout + result.stderr).strip()

    if required_version in magma_version_output:
        print(f"\n✓ MAGMA version {required_version} confirmed (binary executed)")
        print(f"  Full output: {magma_version_output.splitlines()[0]}")
    else:
        raise AssertionError(
            f"MAGMA version mismatch!\n"
            f"  Required: {required_version}\n"
            f"  Found:    {magma_version_output}"
        )

except (OSError, PermissionError) as e:
    # Binary won't execute on this host — almost certainly Linux binary on macOS
    binary_size = magma_binary_path.stat().st_size if magma_binary_path.exists() else 0
    print(f"\n⚠️  Cannot execute MAGMA binary on this host ({platform.system()} {platform.machine()})")
    print(f"    Binary exists ({binary_size:,} bytes) but is the wrong architecture.")
    print(f"    This is expected if you're on macOS with the Linux static binary.")
    print(f"")
    print(f"    Version pin is still enforced via:")
    print(f"      1. Dockerfile build-time check (passed at docker build)")
    print(f"      2. Will re-verify when bash scripts invoke MAGMA via Docker")
    print(f"")
    print(f"    Proceeding with prep notebook (host-side MAGMA is never invoked).")
    magma_version_output = f"deferred-to-docker ({platform.system()} host)"


# %%
# =============================================================================
# 3. Initialize provenance log
# =============================================================================
def log_provenance(message: str, log_path: pathlib.Path = provenance_log) -> None:
    """Append a timestamped message to the run log."""
    ts = datetime.now(timezone.utc).isoformat()
    with open(log_path, 'a') as f:
        f.write(f"[{ts}] {message}\n")

with open(provenance_log, 'w') as f:
    f.write(f"# MAGMA Cell-Type Pipeline — Run Log\n")
    f.write(f"# GWAS: {gwas['name']}\n")
    f.write(f"# Started: {datetime.now(timezone.utc).isoformat()}\n")
    f.write(f"# MAGMA: {magma_version_output}\n")
    f.write(f"# Runtime: {runtime_env}\n")
    f.write(f"# Host: {platform.system()} {platform.machine()}\n\n")

log_provenance(f"01_prep_sumstats started; raw file = {gwas_raw_path.name}")
print(f"Provenance log: {provenance_log}")


# %%
# =============================================================================
# 4. Inspect raw file structure
# =============================================================================
assert gwas_raw_path.exists(), \
    f"\nRAW FILE NOT FOUND: {gwas_raw_path}\n" \
    f"Download from {gwas['download_url']} into {gwas_raw_path.parent}/"

print(f"Inspecting: {gwas_raw_path}")
print(f"File size:  {gwas_raw_path.stat().st_size / 1e6:.1f} MB")

opener = gzip.open if str(gwas_raw_path).endswith('.gz') else open
with opener(gwas_raw_path, 'rt') as f:
    first_lines = [next(f) for _ in range(5)]

print("\n--- First 5 raw lines ---")
for i, line in enumerate(first_lines):
    print(f"[{i}] {line.rstrip()[:200]}")

df_peek = pd.read_csv(gwas_raw_path, sep=None, engine='python', nrows=1000)
print(f"\n--- Sample shape: {df_peek.shape} ---")
print(f"\nColumns: {list(df_peek.columns)}")
print(f"\nDtypes:\n{df_peek.dtypes}")
print(f"\n--- First 5 rows ---")
print(df_peek.head())


# %% [markdown]
# ### STOP HERE on first run
#
# Look at column names above. If they don't match the `columns:` block in
# `00_config.yaml`, update the config and re-run from cell 1.

# %%
# =============================================================================
# 5. Validate column mapping
# =============================================================================
col_map = gwas['columns']
required = ['rsid_col', 'pvalue_col', 'chr_col', 'bp_col']
missing  = [c for c in required if not col_map.get(c)]
if missing:
    raise ValueError(
        f"Required column mappings not set: {missing}\n"
        f"Update `columns:` in 00_config.yaml to match observed columns: "
        f"{list(df_peek.columns)}"
    )

print("Column mapping:")
for k, v in col_map.items():
    print(f"  {k:12s} → {v}")

specified = {v for v in col_map.values() if v is not None}
missing_in_file = specified - set(df_peek.columns)
if missing_in_file:
    raise ValueError(f"Columns in config but not in file: {missing_in_file}")
print("\n✓ All specified columns are present in the file")


# %%
# =============================================================================
# 6. Load full file
# =============================================================================
print("Loading full sumstats (may take 1-3 min for ~10M SNPs)...")
df = pd.read_csv(gwas_raw_path, sep=None, engine='python')
print(f"Loaded {len(df):,} rows × {df.shape[1]} columns "
      f"({df.memory_usage(deep=True).sum() / 1e9:.2f} GB)")
log_provenance(f"loaded raw sumstats: {len(df):,} rows")


# %%
# =============================================================================
# 7. Rename columns to Duncan's expected schema
# =============================================================================
rename_map = {
    col_map['rsid_col']:   'ID',
    col_map['pvalue_col']: 'PVAL',
    col_map['chr_col']:    'CHR',
    col_map['bp_col']:     'BP',
}
if col_map.get('a1_col'):   rename_map[col_map['a1_col']]   = 'A1'
if col_map.get('a2_col'):   rename_map[col_map['a2_col']]   = 'A2'
if col_map.get('maf_col'):  rename_map[col_map['maf_col']]  = 'MAF'
if col_map.get('info_col'): rename_map[col_map['info_col']] = 'INFO'

df = df.rename(columns=rename_map)
print("After renaming:", list(df.columns))


# %%
# =============================================================================
# 8. Compute NEFF per SNP
# =============================================================================
if col_map.get('n_col'):
    df = df.rename(columns={col_map['n_col']: 'NEFF'})
    print(f"Using per-SNP N column ('NEFF'); "
          f"min={df['NEFF'].min():.0f}, median={df['NEFF'].median():.0f}, max={df['NEFF'].max():.0f}")
else:
    n_eff = 4.0 / (1.0 / gwas['n_cases'] + 1.0 / gwas['n_controls'])
    df['NEFF'] = n_eff
    print(f"No per-SNP N column; using constant NEFF = {n_eff:.0f}")


# %%
# =============================================================================
# 9. Filter
# =============================================================================
print(f"Starting rows:                       {len(df):>12,}")

n_before = len(df)
df = df.dropna(subset=['ID', 'PVAL', 'CHR', 'BP'])
print(f"After dropping NA in ID/PVAL/CHR/BP: {len(df):>12,}  (-{n_before - len(df):,})")

df['CHR'] = df['CHR'].astype(str).str.replace('chr', '', case=False)
keep_chrs = [str(c) for c in config['qc']['keep_chromosomes']]
n_before = len(df)
df = df[df['CHR'].isin(keep_chrs)].copy()
print(f"After chromosome filter:             {len(df):>12,}  (-{n_before - len(df):,})")

if 'A1' in df.columns and 'A2' in df.columns and config['qc']['drop_indels']:
    n_before = len(df)
    is_snv = (df['A1'].str.len() == 1) & (df['A2'].str.len() == 1) \
           & df['A1'].isin(['A', 'C', 'G', 'T']) & df['A2'].isin(['A', 'C', 'G', 'T'])
    df = df[is_snv].copy()
    print(f"After indel filter (SNV only):       {len(df):>12,}  (-{n_before - len(df):,})")
else:
    print(f"Indel filter SKIPPED (no A1/A2 or drop_indels=false)")

if 'MAF' in df.columns and config['qc']['maf_min'] is not None:
    n_before = len(df)
    df = df[df['MAF'] >= config['qc']['maf_min']].copy()
    print(f"After MAF >= {config['qc']['maf_min']}:                {len(df):>12,}  (-{n_before - len(df):,})")

if 'INFO' in df.columns and config['qc']['info_min'] is not None:
    n_before = len(df)
    df = df[df['INFO'] >= config['qc']['info_min']].copy()
    print(f"After INFO >= {config['qc']['info_min']}:               {len(df):>12,}  (-{n_before - len(df):,})")

n_before = len(df)
df = df.sort_values('PVAL').drop_duplicates(subset='ID', keep='first')
print(f"After dropping duplicate rsIDs:      {len(df):>12,}  (-{n_before - len(df):,})")

n_before = len(df)
df = df[(df['PVAL'] > 0) & (df['PVAL'] <= 1)].copy()
print(f"After valid p-value (0 < P <= 1):    {len(df):>12,}  (-{n_before - len(df):,})")

df['CHR_sort'] = df['CHR'].replace({'X': 23, 'Y': 24}).astype(int)
df['BP'] = df['BP'].astype(int)
df = df.sort_values(['CHR_sort', 'BP']).drop(columns='CHR_sort').reset_index(drop=True)

print(f"\nFinal row count:                     {len(df):>12,}")
log_provenance(f"after filtering: {len(df):,} rows")


# %%
# =============================================================================
# 10. Genomic inflation factor (λ_GC)
# =============================================================================
chi2_stats = chi2.isf(df['PVAL'].values, df=1)
lambda_gc = np.nanmedian(chi2_stats) / 0.4549
print(f"λ_GC = {lambda_gc:.4f}")
print("  Watanabe 2022 reports λ ≈ 1.5 (full meta); UKB-only expected ~1.1-1.3.")
print("  >1.5 generally fine for polygenic traits; LDSC intercept >1.1 would suggest confounding.")
log_provenance(f"lambda_GC = {lambda_gc:.4f}")


# %%
# =============================================================================
# 11. SNP-level proxy check on Jansen 2019 sentinel genes
# =============================================================================
gene_loc = pd.read_csv(
    project_root / config['paths']['gene_loc'],
    sep='\t', header=None,
    names=['ENTREZ', 'CHR', 'START', 'STOP', 'STRAND', 'SYMBOL'],
    dtype={'ENTREZ': str, 'CHR': str}
)
print(f"Loaded {len(gene_loc):,} gene locations")


def best_snp_near_gene(entrez: str, window_kb: int = 50) -> Optional[dict]:
    row = gene_loc[gene_loc['ENTREZ'] == entrez]
    if row.empty:
        return None
    chrom = row.iloc[0]['CHR']
    start = row.iloc[0]['START'] - window_kb * 1000
    stop  = row.iloc[0]['STOP']  + window_kb * 1000
    nearby = df[(df['CHR'] == str(chrom)) & (df['BP'] >= start) & (df['BP'] <= stop)]
    if nearby.empty:
        return None
    best = nearby.nsmallest(1, 'PVAL').iloc[0]
    return {'chr': chrom, 'best_snp': best['ID'], 'best_p': best['PVAL'],
            'n_in_window': len(nearby)}


print("\nSNP-level proxy check (best SNP within 50kb of gene):")
print("-" * 80)
snp_check_results = []
for sentinel in gwas['sentinel_genes']:
    result = best_snp_near_gene(sentinel['entrez'])
    if result is None:
        print(f"  {sentinel['symbol']:8s}: NO SNPS IN WINDOW (Entrez {sentinel['entrez']})")
        snp_check_results.append((sentinel['symbol'], False))
        continue
    expected = sentinel['snp_proxy_p_max']
    passed = result['best_p'] <= expected
    status = "✓" if passed else "⚠"
    print(f"  {status} {sentinel['symbol']:8s} chr{result['chr']:>2s}: "
          f"best={result['best_snp']:>15s} P={result['best_p']:.2e} "
          f"(threshold ≤ {expected:.0e})")
    snp_check_results.append((sentinel['symbol'], passed))

n_passed = sum(1 for _, p in snp_check_results if p)
print(f"\nSNP-level proxy: {n_passed}/{len(snp_check_results)} sentinels show signal")
print("Note: this is a weak check. Formal Bonferroni assertion runs in 02_qc_gene_results")
print("      after MAGMA gene analysis.")
log_provenance(f"SNP-level proxy: {n_passed}/{len(snp_check_results)} sentinels passed")


# %%
# =============================================================================
# 12. Manhattan plot
# =============================================================================
fig, ax = plt.subplots(figsize=(14, 4))
df_plot = df.copy()
df_plot['CHR_num'] = df_plot['CHR'].replace({'X': 23}).astype(int)
df_plot = df_plot.sort_values(['CHR_num', 'BP'])
df_plot['neglog10P'] = -np.log10(df_plot['PVAL'].clip(lower=1e-300))

sig = df_plot[df_plot['PVAL'] < 0.01]
nom = df_plot[df_plot['PVAL'] >= 0.01].sample(
    min(200000, (df_plot['PVAL'] >= 0.01).sum()), random_state=0)
sub = pd.concat([sig, nom]).sort_values(['CHR_num', 'BP'])

chr_offsets = {}
running = 0
for c in range(1, 24):
    chr_offsets[c] = running
    in_chr = sub[sub['CHR_num'] == c]
    if not in_chr.empty:
        running += in_chr['BP'].max()
sub['x_cum'] = sub.apply(lambda r: chr_offsets[r['CHR_num']] + r['BP'], axis=1)

colors = ['#1f77b4', '#ff7f0e']
for c in range(1, 24):
    cd = sub[sub['CHR_num'] == c]
    if not cd.empty:
        ax.scatter(cd['x_cum'], cd['neglog10P'], c=colors[c % 2], s=2, alpha=0.6)
ax.axhline(-np.log10(5e-8), color='red', ls='--', lw=0.8, label='GWS (5e-8)')
ax.set_xlabel('Chromosome')
ax.set_ylabel(r'$-\log_{10}(P)$')
ax.set_title(f"Manhattan — {gwas['name']}")
ax.legend()
midpoints = [(c, sub[sub['CHR_num']==c]['x_cum'].median())
             for c in range(1, 24) if not sub[sub['CHR_num']==c].empty]
ax.set_xticks([m for _, m in midpoints])
ax.set_xticklabels(['X' if c==23 else str(c) for c, _ in midpoints], fontsize=8)
plt.tight_layout()

out_fig = project_root / 'figures' / f"01_manhattan_{gwas['output_folder']}.png"
out_fig.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_fig, dpi=150, bbox_inches='tight')
print(f"Saved: {out_fig}")
plt.close(fig)


# %%
# =============================================================================
# 13. QQ plot
# =============================================================================
fig, ax = plt.subplots(figsize=(5, 5))
n_sample = min(500000, len(df))
sample = df.sample(n_sample, random_state=0)
observed = -np.log10(np.sort(sample['PVAL'].values))
expected_pvals = -np.log10(np.linspace(1/n_sample, 1, n_sample))
ax.scatter(expected_pvals, observed, s=2, alpha=0.5)
lim = max(observed.max(), expected_pvals.max())
ax.plot([0, lim], [0, lim], 'r--', lw=0.8, label='y=x')
ax.set_xlabel(r'Expected $-\log_{10}(P)$')
ax.set_ylabel(r'Observed $-\log_{10}(P)$')
ax.set_title(f"QQ — {gwas['name']}\nλ_GC = {lambda_gc:.3f}")
ax.legend()
plt.tight_layout()
out_qq = project_root / 'figures' / f"01_qqplot_{gwas['output_folder']}.png"
plt.savefig(out_qq, dpi=150, bbox_inches='tight')
print(f"Saved: {out_qq}")
plt.close(fig)


# %%
# =============================================================================
# 14. Write output files in Duncan format
# =============================================================================
output_cols = ['ID', 'CHR', 'BP', 'PVAL', 'NEFF']
if 'A1' in df.columns: output_cols.append('A1')
if 'A2' in df.columns: output_cols.append('A2')

df['NEFF'] = df['NEFF'].round().astype(int)

df[output_cols].to_csv(output_sumstats_path, sep='\t', index=False)
print(f"✓ Wrote sumstats: {output_sumstats_path}")
print(f"  ({len(df):,} rows × {len(output_cols)} columns, "
      f"{output_sumstats_path.stat().st_size / 1e6:.1f} MB)")

df[['ID', 'CHR', 'BP']].to_csv(output_snploc_path, sep='\t', index=False, header=False)
print(f"✓ Wrote snploc:   {output_snploc_path}")
log_provenance(f"wrote sumstats ({len(df):,} rows) and snploc")


# %%
# =============================================================================
# 15. PRE-FLIGHT OVERLAP GATE — assert ≥70% overlap with g1000_eur.bim
# =============================================================================
g1000_bim_path = project_root / (config['paths']['g1000_ref_prefix'] + '.bim')
assert g1000_bim_path.exists(), \
    f"\nReference panel not found: {g1000_bim_path}\n" \
    f"Download g1000_eur.{{bed,bim,fam,synonyms}} into auxfiles/"

print(f"Loading SNP IDs from reference panel: {g1000_bim_path}")
print(f"  ({g1000_bim_path.stat().st_size / 1e6:.0f} MB; may take 30-60s)")

ref_rsids = pd.read_csv(g1000_bim_path, sep='\t', header=None, usecols=[1],
                       names=['ID'], dtype={'ID': str})
print(f"Reference panel: {len(ref_rsids):,} SNPs")

our_rsids = set(df['ID'].values)
ref_set   = set(ref_rsids['ID'].values)
overlap   = our_rsids & ref_set
overlap_frac = len(overlap) / len(our_rsids)

print(f"\nOverlap with reference panel:")
print(f"  Our SNPs:       {len(our_rsids):>12,}")
print(f"  Ref SNPs:       {len(ref_set):>12,}")
print(f"  Overlap:        {len(overlap):>12,}")
print(f"  Overlap fraction: {overlap_frac:.2%}")

threshold = config['qc']['min_g1000_overlap_frac']
log_provenance(f"g1000 overlap: {overlap_frac:.2%} ({len(overlap):,}/{len(our_rsids):,})")

if overlap_frac < threshold:
    log_provenance(f"FAIL: overlap {overlap_frac:.2%} below threshold {threshold:.0%}")
    raise AssertionError(
        f"\n\n*** PRE-FLIGHT OVERLAP GATE FAILED ***\n"
        f"Only {overlap_frac:.2%} of our SNPs are in g1000_eur reference.\n"
        f"Required threshold: ≥ {threshold:.0%}.\n\n"
        f"Probable causes:\n"
        f"  1. Build mismatch: sumstats are hg38 but reference is hg19. Liftover needed.\n"
        f"  2. rsID staleness: sumstats use old/deprecated rsIDs. Update via dbSNP merge table.\n"
        f"  3. Wrong reference panel: e.g. EAS or AFR rsIDs vs. EUR file.\n"
        f"  4. ID format mismatch: 'chr:pos:ref:alt' vs 'rsNNN'. Standardize before continuing.\n\n"
        f"Do NOT proceed to MAGMA until this passes."
    )
else:
    print(f"\n✓ Overlap gate PASSED ({overlap_frac:.2%} >= {threshold:.0%})")
    log_provenance(f"PASS: overlap {overlap_frac:.2%} meets threshold")


# %%
# =============================================================================
# 16. Summary
# =============================================================================
print("\n" + "=" * 70)
print("01_prep_sumstats — DONE")
print("=" * 70)
print(f"GWAS:              {gwas['name']}")
print(f"MAGMA version:     {required_version} ({magma_version_output})")
print(f"Final SNPs:        {len(df):,}")
print(f"λ_GC:              {lambda_gc:.3f}")
print(f"SNP-proxy passes:  {n_passed}/{len(snp_check_results)}")
print(f"g1000 overlap:     {overlap_frac:.2%}")
print(f"\nOutputs:")
print(f"  • {output_sumstats_path}")
print(f"  • {output_snploc_path}")
print(f"  • {provenance_log}")
print(f"\nNext: run 02_run_magma_step1and2.sh (in Docker)")

log_provenance("01_prep_sumstats completed successfully")
