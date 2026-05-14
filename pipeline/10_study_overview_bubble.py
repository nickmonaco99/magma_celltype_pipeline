"""
Study overview figure for cross-trait MAGMA gene-level analysis.

Scientific paper aesthetic via SciencePlots:
  - white background, sans-serif (Nature-style)
  - NPG (Nature Publishing Group) palette for trait categories
  - thin black spines, tick marks inward, minor ticks visible
  - figure caption written for placement BELOW the figure in the manuscript

X-axis: gene-level Bonferroni-significant gene hits
Y-axis: studies, ordered from largest N_eff (top) to smallest (bottom)
Bubble area: effective sample size

Outputs PNG + SVG (vector) to figures/.

Install once: pip install SciencePlots
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots  # registers 'science', 'nature', 'ieee' etc. styles
from matplotlib.lines import Line2D
from pathlib import Path

# Set style: 'science' for the scientific baseline, 'nature' for sans-serif,
# 'no-latex' so the script runs without a LaTeX install on disk.
plt.style.use(['science', 'nature', 'no-latex'])

# Force off the top/right spines AND ticks globally.
# SciencePlots' 'science' style turns these on; we override here.
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.top'] = False
plt.rcParams['ytick.right'] = False

PROJECT_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = PROJECT_ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)


# ----------------------------------------------------------------------------
# NPG palette (Nature Publishing Group, from R's ggsci::pal_npg)
# ----------------------------------------------------------------------------

NPG = {
    'red':    '#E64B35',
    'cyan':   '#4DBBD5',
    'teal':   '#00A087',
    'navy':   '#3C5488',
    'peach':  '#F39B7F',
    'slate':  '#8491B4',
    'mint':   '#91D1C2',
    'brown':  '#7E6148',
    'khaki':  '#B09C85',
}

CAT_COLORS = {
    'Brain-direct':         NPG['navy'],
    'Sleep (mixed)':        NPG['cyan'],
    'Body metabolic / GI':  NPG['red'],
    'Reproductive':         NPG['peach'],
}


# ----------------------------------------------------------------------------
# Trait registry
# ----------------------------------------------------------------------------

TRAITS = [
    {'key': 'scz_pgc3_2022',           'display': 'Schizophrenia',
     'source': 'PGC3 2022',  'category': 'Brain-direct',
     'genes_out': 'results/SCZ_PGC3/scz_pgc3_2022.no_heading.step2.genes.out'},
    {'key': 'chronotype_loh2018_ukb',  'display': 'Chronotype',
     'source': 'Loh 2018, UKB',  'category': 'Brain-direct',
     'genes_out': 'results/Chronotype_Loh2018/chronotype_loh2018_ukb.no_heading.step2.genes.out'},
    {'key': 'insomnia_ukb',            'display': 'Insomnia',
     'source': 'UKB',  'category': 'Brain-direct',
     'genes_out': 'results/INS_UKB/insomnia_ukb.no_heading.step2.genes.out'},
    {'key': 'bmi_yengo2018',           'display': 'BMI',
     'source': 'Yengo 2018',  'category': 'Brain-direct',
     'genes_out': 'Akshita BMI data/bmi_yengo2018.no_heading.step2.genes.out'},
    {'key': 'migraine_finngen',        'display': 'Migraine',
     'source': 'FinnGen',  'category': 'Brain-direct',
     'genes_out': 'Armaan Snapshots/migraine_finngen_v1_2026-05-13/gene_level/genes.out'},
    {'key': 't2d_xue2018',             'display': 'Type 2 diabetes',
     'source': 'Xue 2018',  'category': 'Body metabolic / GI',
     'genes_out': 'results/T2D_Xue2018/t2d_xue2018.no_heading.step2.genes.out'},
    {'key': 'ibs_eijsbouts2021',       'display': 'IBS',
     'source': 'Eijsbouts 2021',  'category': 'Body metabolic / GI',
     'genes_out': 'results/IBS_Eijsbouts2021/ibs_eijsbouts2021.no_heading.step2.genes.out'},
    {'key': 'sleep_apnea_finngen',     'display': 'Sleep apnea',
     'source': 'FinnGen',  'category': 'Sleep (mixed)',
     'genes_out': 'Armaan Snapshots/sleep_apnea_finngen_v1_2026-05-13/gene_level/genes.out'},
    {'key': 'sleep_combined_finngen',  'display': 'Sleep combined',
     'source': 'FinnGen',  'category': 'Sleep (mixed)',
     'genes_out': 'Armaan Snapshots/sleep_combined_finngen_v1_2026-05-13/gene_level/genes.out'},
    {'key': 'menorrhagia_finngen',     'display': 'Menorrhagia',
     'source': 'FinnGen',  'category': 'Reproductive',
     'genes_out': 'Armaan Snapshots/menorrhagia_finngen_v1_2026-05-13/gene_level/genes.out'},
    {'key': 'menst_excessive_finngen', 'display': 'Heavy menses',
     'source': 'FinnGen',  'category': 'Reproductive',
     'genes_out': 'Armaan Snapshots/menst_excessive_finngen_v1_2026-05-13/gene_level/genes.out'},
    {'key': 'menst_irregular_finngen', 'display': 'Irregular menses',
     'source': 'FinnGen',  'category': 'Reproductive',
     'genes_out': 'Armaan Snapshots/menst_irregular_finngen_v1_2026-05-13/gene_level/genes.out'},
]


# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------

print('Extracting gene-level statistics per trait...\n')
rows = []
for t in TRAITS:
    p = PROJECT_ROOT / t['genes_out']
    if not p.exists():
        print(f"  [skip] {t['key']}: not found at {t['genes_out']}")
        continue
    df = pd.read_csv(p, sep=r'\s+', comment='#')
    if 'P' not in df.columns or 'N' not in df.columns:
        print(f"  [skip] {t['key']}: P or N missing. Cols: {list(df.columns)[:10]}")
        continue
    n_genes_total = len(df)
    gene_bonf = 0.05 / n_genes_total
    n_genes_sig = int((df['P'] < gene_bonf).sum())
    n_eff = float(df['N'].median())
    rows.append({
        'key': t['key'], 'display': t['display'], 'source': t['source'],
        'category': t['category'], 'color': CAT_COLORS[t['category']],
        'n_genes_total': n_genes_total, 'n_genes_sig': n_genes_sig,
        'n_eff': n_eff,
    })
    print(f"  {t['display']:<24s}  total: {n_genes_total:>5d}  "
          f"sig: {n_genes_sig:>5d}  N_eff: {n_eff:>9,.0f}")

if not rows:
    raise SystemExit('No data found. Check genes.out file paths.')

# Sort descending by gene hits and assign y_pos so highest hits sits at TOP
data = pd.DataFrame(rows).sort_values('n_genes_sig', ascending=False).reset_index(drop=True)
data['y_pos'] = len(data) - 1 - np.arange(len(data))
data.to_csv(FIG_DIR / 'study_overview_data.csv', index=False)


# ----------------------------------------------------------------------------
# Plot — single font scale, generous left margin, no clipping
# ----------------------------------------------------------------------------

BASE = 16
TRAIT_FONT       = BASE + 2  # 18 bold (trait name)
SOURCE_FONT      = BASE - 3  # 13 italic (Yengo 2018 etc.)
XTICK_FONT       = BASE      # 16 (axis tick numbers)
XLABEL_FONT      = BASE + 2  # 18 (axis title)
BUBBLE_IN_FONT   = BASE - 1  # 15 white bold (inside large bubble)
BUBBLE_OUT_FONT  = BASE - 2  # 14 (beside small bubble)
LEGEND_FONT      = BASE - 1  # 15
LEGEND_MARKER    = 18

fig, ax = plt.subplots(figsize=(14, 10))

# Bubble area: linear in N_eff so area faithfully encodes N
n_max, n_min = data['n_eff'].max(), data['n_eff'].min()
AREA_MAX, AREA_MIN = 1400, 160
data['bubble_area'] = AREA_MIN + (AREA_MAX - AREA_MIN) * \
                      (data['n_eff'] - n_min) / (n_max - n_min)

# Lollipop stems
for _, r in data.iterrows():
    ax.hlines(y=r['y_pos'], xmin=0, xmax=r['n_genes_sig'],
              color='0.75', linewidth=0.7, alpha=1.0, zorder=1)

# Bubbles
for _, r in data.iterrows():
    ax.scatter(r['n_genes_sig'], r['y_pos'],
               s=r['bubble_area'], c=r['color'],
               alpha=0.92, edgecolors='white', linewidths=0.8, zorder=3)

# In-bubble N_eff labels for large bubbles; beside for small
for _, r in data.iterrows():
    label = f"{r['n_eff']/1000:.0f}k"
    if r['bubble_area'] >= 600:
        ax.annotate(label, xy=(r['n_genes_sig'], r['y_pos']),
                    ha='center', va='center', fontsize=BUBBLE_IN_FONT,
                    color='white', fontweight='bold', zorder=4)
    else:
        r_pts = np.sqrt(r['bubble_area'] / np.pi)
        ax.annotate(label, xy=(r['n_genes_sig'], r['y_pos']),
                    xytext=(r_pts + 3, 0), textcoords='offset points',
                    ha='left', va='center', fontsize=BUBBLE_OUT_FONT,
                    color='black', zorder=4)

# Y-axis: trait name as the tick label; source line as a separate annotation
ax.set_yticks(data['y_pos'])
ax.set_yticklabels(data['display'], fontsize=TRAIT_FONT, fontweight='bold')

# Source labels under each trait. Anchor in axes coords (x) + data coords (y)
# so they always sit at the same horizontal x as the y-tick labels.
for _, r in data.iterrows():
    ax.annotate(r['source'], xy=(0, r['y_pos']),
                xytext=(-8, -14), textcoords='offset points',
                xycoords=('axes fraction', 'data'),
                ha='right', va='top', fontsize=SOURCE_FONT,
                color='0.4', fontstyle='italic')

# X-axis
xmax = data['n_genes_sig'].max()
ax.set_xlim(-xmax * 0.04, xmax * 1.10)
ax.set_xlabel('Bonferroni-significant gene hits (MAGMA)',
              fontsize=XLABEL_FONT, labelpad=12)
ax.tick_params(axis='x', labelsize=XTICK_FONT)
ax.tick_params(axis='y', length=0)
ax.set_ylim(-0.7, len(data) - 0.3)

# Remove ALL spines and ticks for a clean borderless look
for side in ('top', 'right', 'bottom', 'left'):
    ax.spines[side].set_visible(False)
ax.tick_params(top=False, right=False, bottom=False, left=False, which='both')

# Category legend (bottom-right, frameless)
cat_handles = [Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=c, markeredgecolor='white',
                       markersize=LEGEND_MARKER, label=cat)
               for cat, c in CAT_COLORS.items()
               if cat in data['category'].unique()]
ax.legend(handles=cat_handles, loc='lower right',
          frameon=False, fontsize=LEGEND_FONT,
          labelspacing=0.9, handletextpad=0.6,
          bbox_to_anchor=(1.0, 0.0))

# Caption to stdout for the manuscript.
print('\nSuggested caption (paste below the figure in your manuscript):\n'
      '  Figure 1. Bonferroni-significant gene hits per GWAS (MAGMA gene-level analysis),\n'
      '  plotted against effective sample size (bubble area). Studies ordered top-to-bottom\n'
      '  by descending gene-hit count. Effective N = 4 / (1/N_case + 1/N_ctrl) for case-control\n'
      '  traits, raw N for continuous traits. For polygenic case-control traits, gene\n'
      '  discovery is also modulated by trait prevalence and effect-size distribution\n'
      '  (Gratten et al. 2014).\n')

# Generous left margin so even "Irregular menses" at 18pt bold has room.
plt.subplots_adjust(left=0.28, right=0.95, top=0.95, bottom=0.10)

png_out = FIG_DIR / 'study_overview_bubble.png'
svg_out = FIG_DIR / 'study_overview_bubble.svg'
plt.savefig(png_out, dpi=300, bbox_inches='tight')
plt.savefig(svg_out, bbox_inches='tight')  # vector for journal submissions
print(f'Saved -> {png_out}')
print(f'Saved -> {svg_out}')
