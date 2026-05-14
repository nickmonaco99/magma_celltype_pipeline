"""
Slide figure: Bonferroni -> independent reduction across 3 traits with conditional analysis.

Columns: Trait | Passed Bonferroni | Independent | Reduction
"""
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import scienceplots
    plt.style.use(['science', 'no-latex'])
except ImportError:
    print("Tip: pip install SciencePlots for publication style. Using default.")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = PROJECT_ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

# ---- Data ----
TRAITS = [
    {'name': 'SCZ',         'bonferroni': 105, 'independent': 9,  'color': '#c2185b'},
    {'name': 'Chronotype',  'bonferroni': 48,  'independent': 11, 'color': '#7b1fa2'},
    {'name': 'BMI',         'bonferroni': 81,  'independent': 12, 'color': '#388e3c'},
]

# ---- Layout constants ----
# Vertical layout: title strip (top) | header row | 3 data rows | total row | footer pad
HEADER_Y = 0.0
ROW0_Y   = 1.4
ROW_DY   = 1.0
TOTAL_Y  = ROW0_Y + ROW_DY * len(TRAITS) + 0.3

# Column x-positions (figure-relative, 0-100)
x_trait      = 12
x_bonf = 38
x_indep      = 62
x_reduce     = 86

# ---- Figure ----
fig, ax = plt.subplots(figsize=(10, 4.8), dpi=200)
ax.set_xlim(0, 100)
ax.set_ylim(TOTAL_Y + 1.0, HEADER_Y - 1.4)  # inverted, with padding
ax.axis('off')

# Header row (single-line, won't wrap into row below)
ax.text(x_trait,      HEADER_Y, 'Trait',                    fontsize=12, weight='bold', ha='center', va='center')
ax.text(x_bonf, HEADER_Y, 'Passed Bonferroni',  fontsize=12, weight='bold', ha='center', va='center')
ax.text(x_indep,      HEADER_Y, 'Independent',              fontsize=12, weight='bold', ha='center', va='center')
ax.text(x_reduce,     HEADER_Y, 'Reduction',                fontsize=12, weight='bold', ha='center', va='center')

# Underline below header
underline_y = HEADER_Y + 0.55
ax.plot([2, 98], [underline_y, underline_y], color='black', linewidth=1.5)

# Data rows
for i, t in enumerate(TRAITS):
    y = ROW0_Y + i * ROW_DY
    reduction = t['bonferroni'] / t['independent']

    ax.text(x_trait, y, t['name'], fontsize=13, weight='bold',
            ha='center', va='center', color=t['color'])
    ax.text(x_bonf, y, f"{t['bonferroni']}", fontsize=12, ha='center', va='center')
    ax.text(x_indep,      y, f"{t['independent']}", fontsize=12,
            weight='bold', color=t['color'], ha='center', va='center')
    ax.text(x_reduce,     y, f"{reduction:.0f}×", fontsize=13, weight='bold',
            color='#444444', ha='center', va='center')

    if i < len(TRAITS) - 1:
        sep_y = y + ROW_DY / 2
        ax.plot([2, 98], [sep_y, sep_y], color='lightgrey', linewidth=0.5)

# Totals row, separated by black rule
ax.plot([2, 98], [TOTAL_Y - 0.4, TOTAL_Y - 0.4], color='black', linewidth=1.0)
total_s = sum(t['bonferroni'] for t in TRAITS)
total_i = sum(t['independent'] for t in TRAITS)
ax.text(x_trait,      TOTAL_Y, 'Total', fontsize=12, weight='bold', ha='center', va='center')
ax.text(x_bonf, TOTAL_Y, f"{total_s}", fontsize=12, weight='bold', ha='center', va='center')
ax.text(x_indep,      TOTAL_Y, f"{total_i}", fontsize=12, weight='bold', ha='center', va='center')
ax.text(x_reduce,     TOTAL_Y, f"{total_s/total_i:.0f}×",
        fontsize=13, weight='bold', color='#444444', ha='center', va='center')

# Title and subtitle live ABOVE the axes via fig.suptitle / fig.text so they don't touch the table
fig.suptitle('Conditional analysis: ~100 correlated hits → ~10 non-redundant cell types',
             fontsize=14, weight='bold', y=0.98)
fig.text(0.5, 0.92,
         'Forward selection on pairwise joint-pairs MAGMA (Duncan et al. 2025)',
         fontsize=10, ha='center', style='italic', color='#666666')

# Manual top padding so title doesn't crash into the header row
plt.subplots_adjust(top=0.85, bottom=0.05, left=0.04, right=0.96)

png_out = FIG_DIR / 'conditional_reduction_table.png'
svg_out = FIG_DIR / 'conditional_reduction_table.svg'
plt.savefig(png_out, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(svg_out, bbox_inches='tight', facecolor='white')
print(f"Saved -> {png_out}")
print(f"Saved -> {svg_out}")
