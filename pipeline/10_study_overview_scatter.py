"""
Study overview figure for cross-trait MAGMA gene-level analysis.

A proper correlation plot:
  x = effective sample size (log scale)
  y = number of Bonferroni-significant gene hits
  dot colour = collapsed trait category (4 groups)
  dot label = trait name + source (placed beside each dot)

Annotations:
  - Spearman rho + p in the top-left
  - reference line: Visscher-style linear scaling, anchored at trait-median N
  - callout on the most prominent outlier (typically SCZ in our data)

Outputs PNG + SVG to figures/.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from scipy.stats import spearmanr
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = PROJECT_ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)


# ----------------------------------------------------------------------------
# Palette
# ----------------------------------------------------------------------------

CAT_COLORS = {
    'Brain-direct':         '#2A4858',   # dark teal-blue
    'Sleep (mixed)':        '#6F7C8E',   # muted slate
    'Body metabolic / GI':  '#C77B4E',   # warm terracotta
    'Reproductive':         '#A23B5D',   # muted burgundy
}
BG_COLOR    = '#FAF7F2'
GRID_COLOR  = '#D9D4CD'
TEXT_COLOR  = '#1F2937'
MUTED_COLOR = '#6B7280'


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
     'genes_out': 'results/Insomnia_UKB/insomnia_ukb.no_heading.step2.genes.out'},
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
    {'key': 'menst_excessive_finngen', 'display': 'Excessive menstruation',
     'source': 'FinnGen',  'category': 'Reproductive',
     'genes_out': 'Armaan Snapshots/menst_excessive_finngen_v1_2026-05-13/gene_level/genes.out'},
    {'key': 'menst_irregular_finngen', 'display': 'Irregular menstruation',
     'source': 'FinnGen',  'category': 'Reproductive',
     'genes_out': 'Armaan Snapshots/menst_irregular_finngen_v1_2026-05-13/gene_level/genes.out'},
]


# Per-trait label placement nudges (in points relative to the dot).
# Hand-tuned because adjustText is finicky and 12 points is small enough to dial in.
# (dx, dy, halign) — flip these to taste once you see the render.
LABEL_OFFSETS = {
    'Schizophrenia':           (10, 6, 'left'),
    'Chronotype':              (10, 4, 'left'),
    'Insomnia':                (10, -10, 'left'),
    'BMI':                     (-10, 8, 'right'),
    'Migraine':                (10, -6, 'left'),
    'Type 2 diabetes':         (10, 6, 'left'),
    'IBS':                     (10, 4, 'left'),
    'Sleep apnea':             (10, 6, 'left'),
    'Sleep combined':          (10, -10, 'left'),
    'Menorrhagia':             (10, 8, 'left'),
    'Excessive menstruation':  (10, -10, 'left'),
    'Irregular menstruation':  (-10, -10, 'right'),
}


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

data = pd.DataFrame(rows)
data.to_csv(FIG_DIR / 'study_overview_data.csv', index=False)


# ----------------------------------------------------------------------------
# Correlation stats
# ----------------------------------------------------------------------------

rho, p_rho = spearmanr(data['n_eff'], data['n_genes_sig'])
print(f"\nSpearman rho (N_eff vs n_genes_sig) = {rho:.2f}  (p = {p_rho:.3g})")


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(11, 7.5), dpi=150)
fig.patch.set_facecolor(BG_COLOR)
ax.set_facecolor(BG_COLOR)

# Scatter
for _, r in data.iterrows():
    ax.scatter(r['n_eff'], r['n_genes_sig'],
               s=110, c=r['color'], alpha=0.92,
               edgecolors='white', linewidths=1.4, zorder=3)

# Trait labels, hand-placed
for _, r in data.iterrows():
    dx, dy, halign = LABEL_OFFSETS.get(r['display'], (10, 6, 'left'))
    ax.annotate(
        f"{r['display']}\n{r['source']}",
        xy=(r['n_eff'], r['n_genes_sig']),
        xytext=(dx, dy), textcoords='offset points',
        ha=halign, va='center',
        fontsize=8.5, color=TEXT_COLOR,
        linespacing=1.25,
        zorder=4,
    )
    # Second line (source) muted
    # (matplotlib doesn't style lines individually in annotate; simplest fix is
    # to add the source as a tiny separate annotation underneath.)

# Optional: faint reference line showing linear scaling anchored at the
# overall median. Disable by setting to False.
SHOW_REF_LINE = True
if SHOW_REF_LINE:
    n_anchor = data['n_eff'].median()
    hits_anchor = data['n_genes_sig'].median()
    # Linear scaling: hits = k * N -> k = hits_anchor / n_anchor
    k = hits_anchor / n_anchor
    x_line = np.logspace(np.log10(data['n_eff'].min() * 0.7),
                          np.log10(data['n_eff'].max() * 1.4), 100)
    y_line = k * x_line
    ax.plot(x_line, y_line, color=MUTED_COLOR, linestyle=':', linewidth=1.2,
            alpha=0.7, zorder=1, label='_nolegend_')
    # Label the line
    x_lab = x_line[int(len(x_line) * 0.85)]
    y_lab = k * x_lab
    ax.annotate('linear scaling\n(anchored at median)',
                xy=(x_lab, y_lab),
                xytext=(0, -22), textcoords='offset points',
                ha='center', va='top',
                fontsize=7.5, color=MUTED_COLOR, fontstyle='italic')

# Outlier callout: identify the largest positive residual from the linear ref
if SHOW_REF_LINE:
    expected = k * data['n_eff']
    residual = data['n_genes_sig'] - expected
    outlier_idx = residual.idxmax()
    outlier = data.loc[outlier_idx]
    # Arrow + text
    ax.annotate(
        f"{outlier['display']} sits far above the trend:\n"
        f"~{int(outlier['n_genes_sig'])} hits at only {outlier['n_eff']/1000:.0f}k N.\n"
        "Dense polygenic architecture, not\nsample size, drives discovery here.",
        xy=(outlier['n_eff'], outlier['n_genes_sig']),
        xytext=(outlier['n_eff'] * 0.55, outlier['n_genes_sig'] * 1.05),
        fontsize=9, color=outlier['color'], fontweight='normal',
        ha='center', va='bottom',
        linespacing=1.4,
        arrowprops=dict(arrowstyle='-', color=outlier['color'],
                        lw=0.8, alpha=0.6,
                        connectionstyle='arc3,rad=0.15'),
        zorder=5,
    )

# Spearman annotation, top-left
ax.text(0.015, 0.97,
        f"Spearman $\\rho$ = {rho:.2f}   (p = {p_rho:.2g})",
        transform=ax.transAxes,
        fontsize=10, color=TEXT_COLOR, fontweight='bold',
        ha='left', va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor=BG_COLOR,
                  edgecolor=GRID_COLOR, linewidth=0.6))

# Axes
ax.set_xscale('log')
ax.set_xlim(data['n_eff'].min() * 0.6, data['n_eff'].max() * 1.6)
ax.set_ylim(-max(data['n_genes_sig']) * 0.04, max(data['n_genes_sig']) * 1.18)

# Custom x-tick labels in "k" units, at scientifically-meaningful round numbers
import matplotlib.ticker as mticker
def k_fmt(x, pos):
    if x >= 1e6:
        return f'{x/1e6:.1f}M'
    if x >= 1e3:
        return f'{x/1e3:.0f}k'
    return f'{x:.0f}'
ax.xaxis.set_major_formatter(mticker.FuncFormatter(k_fmt))
ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=10))
ax.xaxis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=10))
ax.xaxis.set_minor_formatter(mticker.NullFormatter())

ax.set_xlabel('Effective sample size  (log scale)',
              fontsize=10.5, color=TEXT_COLOR, labelpad=10)
ax.set_ylabel('Bonferroni-significant gene hits  (MAGMA)',
              fontsize=10.5, color=TEXT_COLOR, labelpad=10)

# Spines/grid
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color(GRID_COLOR)
    ax.spines[side].set_linewidth(0.8)
ax.tick_params(axis='both', colors=MUTED_COLOR, length=0)
ax.tick_params(axis='both', labelsize=9)
ax.grid(True, axis='both', which='major', color=GRID_COLOR, linewidth=0.6, alpha=0.8, zorder=0)
ax.set_axisbelow(True)

# Title + subtitle
fig.text(0.06, 0.965,
         'Sample size predicts gene discovery, but architecture matters more',
         fontsize=15, fontweight='bold', color=TEXT_COLOR, ha='left', va='top')
fig.text(0.06, 0.925,
         'Bonferroni-significant gene hits (MAGMA) for 12 GWAS, plotted against effective sample size. '
         'Dotted line: linear scaling anchored at the trait-median.',
         fontsize=9.5, color=MUTED_COLOR, ha='left', va='top')

# Category legend (bottom-right, no frame)
cat_handles = [Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=c, markeredgecolor='white',
                       markersize=10, label=cat)
               for cat, c in CAT_COLORS.items()
               if cat in data['category'].unique()]
legend = ax.legend(handles=cat_handles, loc='lower right',
                   frameon=False, fontsize=9, labelspacing=0.6,
                   bbox_to_anchor=(1.0, 0.0))
for text in legend.get_texts():
    text.set_color(TEXT_COLOR)

# Footer
fig.text(0.06, 0.02,
         'Source: MAGMA gene-level analysis. Note: for case-control traits, '
         'discovery is also influenced by trait prevalence and effect-size '
         'distribution (Gratten et al.), so the correlation should be read '
         'as suggestive, not causal.',
         fontsize=8, color=MUTED_COLOR, ha='left', va='bottom', fontstyle='italic')

plt.subplots_adjust(left=0.07, right=0.97, top=0.86, bottom=0.10)

png_out = FIG_DIR / 'study_overview_scatter.png'
svg_out = FIG_DIR / 'study_overview_scatter.svg'
plt.savefig(png_out, dpi=300, bbox_inches='tight', facecolor=BG_COLOR)
plt.savefig(svg_out, bbox_inches='tight', facecolor=BG_COLOR)
print(f'\nSaved figure -> {png_out}')
print(f'Saved figure -> {svg_out}')
