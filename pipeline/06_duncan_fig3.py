"""
Duncan 2025 Fig 3-style multi-panel scatterplots across 12 brain-related GWAS traits.

Reference: Duncan et al., Nat Neurosci 2025 (doi:10.1038/s41593-024-01834-w), Figure 3.

Critical design feature: clusters are reordered along the x-axis by SUPERCLUSTER
(Siletti canonical order, the same 1-31 numbering in Duncan's legend), then by
cluster_id within each supercluster. The x-axis is an INDEX position 0..460, not
the cluster_id itself, and is not labelled. This is what produces the "rainbow
blocks" that make the figure readable.

Aesthetic: black panel background, white axes/text, solid white Bonferroni line,
no grid, no x-ticks, top hit annotated with cluster_name in its supercluster colour.

Outputs:
  - figures/fig3_1x12.{png,pdf}  : single column stack (primary, most Duncan-faithful)
  - figures/fig3_2x6.{png,pdf}   : 2 rows x 6 cols (slide-friendly secondary)

Author: Nick Monaco (rewritten 2026-05-13)
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

N_CLUSTERS = 461
BONF_THRESHOLD = 0.05 / N_CLUSTERS           # P < 1.085e-4
BONF_LOG = -np.log10(BONF_THRESHOLD)         # ~3.964

# Font scale: change BASE and everything scales together
BASE = 14
PANEL_TITLE_FONT  = BASE + 1   # 15 bold
ANNOT_FONT        = BASE - 2   # 12 bold (supercluster labels on points)
YTICK_FONT        = BASE - 4   # 10
YLABEL_FONT       = BASE - 3   # 11
LEGEND_FONT       = BASE - 3   # 11
LEGEND_TITLE_FONT = BASE - 2   # 12
SUPTITLE_FONT     = BASE + 1   # 15

# Siletti canonical supercluster order (1..31). This drives x-axis layout and
# the palette. Order matches the numbered legend in Duncan Fig 3.
SUPERCLUSTER_ORDER = [
    'Miscellaneous',                          # 1
    'Microglia',                              # 2
    'Vascular',                               # 3
    'Fibroblast',                             # 4
    'Oligodendrocyte precursor',              # 5
    'Committed oligodendrocyte precursor',    # 6
    'Oligodendrocyte',                        # 7
    'Bergmann glia',                          # 8
    'Astrocyte',                              # 9
    'Ependymal',                              # 10
    'Choroid plexus',                         # 11
    'Deep-layer near-projecting',             # 12
    'Deep-layer corticothalamic and 6b',      # 13
    'Hippocampal CA1-3',                      # 14
    'Upper-layer intratelencephalic',         # 15
    'Deep-layer intratelencephalic',          # 16
    'Amygdala excitatory',                    # 17
    'Hippocampal CA4',                        # 18
    'Hippocampal dentate gyrus',              # 19
    'Medium spiny neuron',                    # 20
    'Eccentric medium spiny neuron',          # 21
    'Splatter',                               # 22
    'MGE interneuron',                        # 23
    'LAMP5-LHX6 and Chandelier',              # 24
    'CGE interneuron',                        # 25
    'Upper rhombic lip',                      # 26
    'Cerebellar inhibitory',                  # 27
    'Lower rhombic lip',                      # 28
    'Mammillary body',                        # 29
    'Thalamic excitatory',                    # 30
    'Midbrain-derived inhibitory',            # 31
]
SC_RANK = {s: i for i, s in enumerate(SUPERCLUSTER_ORDER)}

CANONICAL_COLS = ['cluster_id', 'cluster_name', 'supercluster', 'P']


# Hand-tuned per-panel label rules for the high-signal traits.
# Each entry is (supercluster_name, n_subclusters_to_show).
# 'mode' = 'stacked'      -> supercluster bold on top, cluster_names smaller below
# 'mode' = 'inline_paren' -> "{supercluster} ({cluster_name})" on one line
# Traits NOT in this dict fall back to the automatic 2-tier rule.
PANEL_LABEL_OVERRIDES = {
    'SCZ': {
        'mode': 'stacked',
        'labels': [
            # (supercluster, n_subclusters, x_offset_pts [optional, default 0])
            ('Eccentric medium spiny neuron', 2),
            ('MGE interneuron',               1),
            ('CGE interneuron',               1, 40),  # nudge right of MGE
        ],
        'show_count': True,
    },
    'Chronotype': {
        'mode': 'inline_paren',
        'labels': [
            ('Splatter',            1),
            ('Medium spiny neuron', 1),
        ],
        'show_count': True,
    },
    'BMI': {
        'mode': 'inline_paren',
        'labels': [
            ('Splatter',          1),
            ('CGE interneuron',   1),
            ('Lower rhombic lip', 1),
        ],
        'show_count': True,
    },
}


# ----------------------------------------------------------------------------
# Trait inventory
# ----------------------------------------------------------------------------

TRAITS = [
    ('SCZ',             'canonical', 'results/snapshots/scz_pgc3_2022_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Chronotype',      'canonical', 'results/snapshots/chronotype_loh2018_ukb_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Insomnia',        'canonical', 'results/snapshots/insomnia_ukb_v2_2026-05-13/annotated/celltype_annotated.csv'),
    ('Sleep apnea',     'canonical', 'Armaan Snapshots/sleep_apnea_finngen_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Sleep combined',  'canonical', 'Armaan Snapshots/sleep_combined_finngen_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Migraine',        'canonical', 'Armaan Snapshots/migraine_finngen_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Menorrhagia',     'canonical', 'Armaan Snapshots/menorrhagia_finngen_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Menst. excess.',  'canonical', 'Armaan Snapshots/menst_excessive_finngen_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('Menst. irreg.',   'canonical', 'Armaan Snapshots/menst_irregular_finngen_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('BMI',             'akshika',   'Akshita BMI data/celltype_results_BMI_2018.csv'),
    ('T2D',             'canonical', 'results/snapshots/t2d_xue2018_v1_2026-05-13/annotated/celltype_annotated.csv'),
    ('IBS',             'canonical', 'results/snapshots/ibs_eijsbouts2021_v1_2026-05-13/annotated/celltype_annotated.csv'),
]

ORDER_1X12 = [t for t, _, _ in TRAITS]
ORDER_2X6 = [
    'SCZ', 'Chronotype', 'Insomnia', 'Sleep apnea', 'Sleep combined', 'Migraine',
    'Menorrhagia', 'Menst. excess.', 'Menst. irreg.', 'BMI', 'T2D', 'IBS',
]
ORDER_4X3 = [
    'SCZ',          'BMI',      'Menorrhagia',
    'Chronotype',   'T2D',      'Menst. excess.',
    'Insomnia',     'IBS',      'Menst. irreg.',
    'Sleep apnea',  'Migraine', 'Sleep combined',
]


# ----------------------------------------------------------------------------
# Schema-aware loaders
# ----------------------------------------------------------------------------

def load_canonical(path: Path, trait: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[CANONICAL_COLS].copy()
    df['trait'] = trait
    return df


def load_akshika(path: Path, trait: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={
        'Cluster ID': 'cluster_id',
        'Cluster name': 'cluster_name',
        'Supercluster': 'supercluster',
    })
    df['cluster_id'] = df['cluster_id'].astype(int)
    df = df[CANONICAL_COLS].copy()
    df['trait'] = trait
    return df


def load_all_traits() -> pd.DataFrame:
    dfs, missing = [], []
    for trait, kind, path in TRAITS:
        p = Path(path)
        if not p.exists():
            missing.append(f'  - {trait}: {p}')
            continue
        loader = load_canonical if kind == 'canonical' else load_akshika
        dfs.append(loader(p, trait))
    if missing:
        print('WARNING: missing files (skipped):')
        print('\n'.join(missing))

    combined = pd.concat(dfs, ignore_index=True)
    combined['logP'] = -np.log10(combined['P'].clip(lower=1e-50))
    combined['sc_rank'] = combined['supercluster'].map(SC_RANK)

    if combined['sc_rank'].isna().any():
        unknown = combined[combined['sc_rank'].isna()]['supercluster'].unique().tolist()
        raise ValueError(f'Unknown supercluster strings (not in canonical order): {unknown}')

    # Reorder per trait by (sc_rank, cluster_id), then assign x_index 0..460
    combined = combined.sort_values(['trait', 'sc_rank', 'cluster_id']).reset_index(drop=True)
    combined['x_index'] = combined.groupby('trait').cumcount()
    return combined


# ----------------------------------------------------------------------------
# Palette
# ----------------------------------------------------------------------------

def build_palette() -> dict:
    """Siletti-style palette:
      - 1 (Miscellaneous): brown/tan
      - 2-11 (glial / non-neuronal): greyscale gradient (visible on black, quiet)
      - 12-31 (neuronal): full rainbow
    """
    palette = {}

    # 1: Miscellaneous - brown/tan
    palette[SUPERCLUSTER_ORDER[0]] = (0.62, 0.38, 0.18)

    # 2-11: glial / non-neuronal - mid-to-light greys
    glial = SUPERCLUSTER_ORDER[1:11]
    greys = np.linspace(0.45, 0.85, len(glial))
    for sc, g in zip(glial, greys):
        palette[sc] = (g, g, g)

    # 12-31: neuronal - rainbow
    neurons = SUPERCLUSTER_ORDER[11:]
    for i, sc in enumerate(neurons):
        palette[sc] = cm.gist_rainbow(i / (len(neurons) - 1))

    return palette


def brighten(rgba, amount: float = 0.45):
    """Blend an RGBA toward white. Used for annotation text so dark reds/oranges
    stay legible against the black panel background."""
    r, g, b = rgba[:3]
    return (r + (1 - r) * amount,
            g + (1 - g) * amount,
            b + (1 - b) * amount)


# ----------------------------------------------------------------------------
# Panel rendering (Duncan-style: black bg, white axes, no x-ticks)
# ----------------------------------------------------------------------------

def render_panel(ax, df_trait: pd.DataFrame, title: str, sc_colors: dict,
                 ymax: float, title_inside: bool = True) -> None:
    sub = df_trait.sort_values('x_index')

    sig_mask = sub['logP'] > BONF_LOG
    non_sig = sub[~sig_mask]
    sig     = sub[ sig_mask]

    # Non-sig points: small and dim (the background)
    if len(non_sig):
        ax.scatter(non_sig['x_index'], non_sig['logP'],
                   c=[sc_colors[s] for s in non_sig['supercluster']],
                   s=11, alpha=0.55, linewidths=0, zorder=2)

    # Sig points: 1.5x bigger, full opacity (the foreground)
    if len(sig):
        ax.scatter(sig['x_index'], sig['logP'],
                   c=[sc_colors[s] for s in sig['supercluster']],
                   s=17, alpha=1.0, linewidths=0, zorder=3)

    # Bonferroni: solid white line
    ax.axhline(BONF_LOG, color='white', lw=0.8, alpha=0.85, zorder=1.5)

    # Labels: per-trait override if present, else automatic 2-tier rule.
    LABEL_CAP = 4
    THRESHOLD_FOR_BROAD = 3

    def _ha_for(x):
        x_frac = x / N_CLUSTERS
        if x_frac < 0.06:
            return 'left'
        if x_frac > 0.94:
            return 'right'
        return 'center'

    override = PANEL_LABEL_OVERRIDES.get(title)

    if override and len(sig):
        # Use hand-tuned per-panel rules
        mode = override['mode']
        for entry in override['labels']:
            # Allow optional 3rd tuple element for x-offset in points
            if len(entry) == 2:
                sc_name, n_subclusters = entry
                x_off = 0
            else:
                sc_name, n_subclusters, x_off = entry

            group = sig[sig['supercluster'] == sc_name]
            if len(group) == 0:
                continue  # no sig clusters in this supercluster for this trait
            top_in_sc = group.loc[group['logP'].idxmax()]
            x, y = top_in_sc['x_index'], top_in_sc['logP']
            ha = _ha_for(x)
            colour = brighten(sc_colors[sc_name])

            if mode == 'inline_paren':
                label = f"{sc_name} ({top_in_sc['cluster_name']})"
                ax.annotate(label, xy=(x, y),
                            xytext=(x_off, 6), textcoords='offset points',
                            ha=ha, fontsize=ANNOT_FONT, fontweight='bold',
                            color=colour, zorder=4)

            elif mode == 'stacked':
                sub_clusters = group.nlargest(n_subclusters, 'logP')
                sub_text = ", ".join(sub_clusters['cluster_name'])
                ax.annotate(sc_name, xy=(x, y),
                            xytext=(x_off, 22), textcoords='offset points',
                            ha=ha, va='bottom',
                            fontsize=ANNOT_FONT, fontweight='bold',
                            color=colour, zorder=4)
                ax.annotate(sub_text, xy=(x, y),
                            xytext=(x_off, 8), textcoords='offset points',
                            ha=ha, va='bottom',
                            fontsize=ANNOT_FONT - 3,
                            color=colour, zorder=4)

    elif len(sig):
        # Automatic 2-tier rule for traits without overrides
        candidates = []
        for sc, group in sig.groupby('supercluster'):
            top_in_sc = group.loc[group['logP'].idxmax()]
            n_in_sc = len(group)
            if n_in_sc >= THRESHOLD_FOR_BROAD:
                label_text = sc
            else:
                label_text = f"{sc} ({top_in_sc['cluster_name']})"
            candidates.append({
                'sc': sc, 'label': label_text,
                'x': top_in_sc['x_index'], 'y': top_in_sc['logP'],
            })
        candidates.sort(key=lambda d: -d['y'])
        candidates = candidates[:LABEL_CAP]
        for c in candidates:
            ax.annotate(
                c['label'],
                xy=(c['x'], c['y']),
                xytext=(0, 6), textcoords='offset points',
                ha=_ha_for(c['x']),
                fontsize=ANNOT_FONT, fontweight='bold',
                color=brighten(sc_colors[c['sc']]),
                zorder=4,
            )

    # Black background, white spines/ticks, no x-ticks, no grid
    ax.set_facecolor('black')
    for spine in ax.spines.values():
        spine.set_color('white')
        spine.set_linewidth(0.6)
    ax.tick_params(axis='y', colors='white', labelsize=YTICK_FONT,
                   length=2.5, width=0.5)
    ax.set_xticks([])
    ax.grid(False)

    # Panel title in the top-left of the panel, white
    if title_inside:
        ax.text(0.012, 0.93, title, transform=ax.transAxes,
                fontsize=PANEL_TITLE_FONT, color='white',
                ha='left', va='top', fontweight='bold')
        # For the high-signal traits, show how many cell types passed Bonferroni
        if override and override.get('show_count'):
            # Position the count to the right of the title; offset depends on
            # title length (each char ~0.012 axes fraction at PANEL_TITLE_FONT)
            count_x = 0.012 + len(title) * 0.013 + 0.01
            ax.text(count_x, 0.93,
                    f"{len(sig)} cell types pass Bonferroni",
                    transform=ax.transAxes,
                    fontsize=ANNOT_FONT - 2, color='white',
                    ha='left', va='top')
    else:
        ax.set_title(title, fontsize=PANEL_TITLE_FONT, color='black',
                     loc='left', pad=2, fontweight='bold')

    ax.set_xlim(-3, N_CLUSTERS + 3)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(r'$-\log_{10}(P)$', fontsize=YLABEL_FONT)


# ----------------------------------------------------------------------------
# Figure assembly
# ----------------------------------------------------------------------------

def build_figure(layout_rows: int, layout_cols: int, traits_order: list,
                 fname: str, combined: pd.DataFrame, sc_colors: dict,
                 figsize: tuple) -> None:
    fig, axes = plt.subplots(layout_rows, layout_cols, figsize=figsize, squeeze=False)
    fig.patch.set_facecolor('white')

    # Per-panel y-max so weak signals (T2D, IBS) stay legible alongside strong ones (SCZ).
    # Stacked labels need more vertical headroom than single-line labels.
    # Floor at 1.5x Bonferroni so an empty panel still shows the threshold line.
    per_panel_ymax = {}
    for t in traits_order:
        sub = combined[combined['trait'] == t]
        override = PANEL_LABEL_OVERRIDES.get(t, {})
        headroom = 1.42 if override.get('mode') == 'stacked' else 1.22
        if len(sub) == 0:
            per_panel_ymax[t] = BONF_LOG * 1.5
        else:
            per_panel_ymax[t] = max(sub['logP'].max() * headroom, BONF_LOG * 1.5)

    for idx, trait in enumerate(traits_order):
        r, c = divmod(idx, layout_cols)
        ax = axes[r, c]
        sub = combined[combined['trait'] == trait]
        if len(sub) == 0:
            ax.set_visible(False)
            continue
        render_panel(ax, sub, trait, sc_colors, per_panel_ymax[trait])

    for idx in range(len(traits_order), layout_rows * layout_cols):
        r, c = divmod(idx, layout_cols)
        axes[r, c].set_visible(False)

    # Numbered bottom legend (1..31) in 4 columns
    handles = [plt.scatter([], [], c=[sc_colors[sc]], s=40, marker='s')
               for sc in SUPERCLUSTER_ORDER]
    labels = [f'{i:>2}  {sc}' for i, sc in enumerate(SUPERCLUSTER_ORDER, start=1)]
    fig.legend(handles=handles, labels=labels,
               loc='lower center', bbox_to_anchor=(0.5, 0.005),
               fontsize=LEGEND_FONT, frameon=False, ncol=4,
               columnspacing=1.4, handletextpad=0.5,
               title=f'Superclusters ({len(SUPERCLUSTER_ORDER)} total) encompassing the {N_CLUSTERS} cell types',
               title_fontsize=LEGEND_TITLE_FONT)

    # Small top-left figure caption
    fig.suptitle('Cell-type associations across 12 brain-related phenotypes',
                 fontsize=SUPTITLE_FONT, fontweight='bold', x=0.01, y=0.995, ha='left')

    # Reserve room for the 8-row legend at bottom and the caption at top.
    # Legend has a fixed pixel height, so shorter figures need a larger fractional margin.
    if layout_rows >= 6:
        bottom = 0.13
    elif layout_rows >= 4:
        bottom = 0.17
    else:
        bottom = 0.22
    plt.tight_layout(rect=[0.01, bottom, 0.99, 0.97])

    Path('figures').mkdir(exist_ok=True)
    for ext in ('png', 'pdf'):
        out = Path('figures') / f'{fname}.{ext}'
        plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {out}')
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    print('Loading trait data...')
    combined = load_all_traits()
    print(f'  Loaded {combined["trait"].nunique()} traits, {len(combined):,} rows')

    print('\nPer-trait summary:')
    for trait, _, _ in TRAITS:
        sub = combined[combined['trait'] == trait]
        if len(sub) == 0:
            print(f'  {trait:18s}  MISSING')
            continue
        n_sig = int((sub['logP'] > BONF_LOG).sum())
        top = sub.loc[sub['logP'].idxmax()]
        print(f'  {trait:18s}  hits={n_sig:3d}  top={top["cluster_name"]:28s}  P={top["P"]:.2e}')

    sc_colors = build_palette()

    print('\nBuilding figures...')

    # Primary: single-column stack, most Duncan-faithful (wide panels, low height)
    build_figure(
        layout_rows=12, layout_cols=1,
        traits_order=ORDER_1X12,
        fname='fig3_1x12',
        combined=combined, sc_colors=sc_colors,
        figsize=(15, 26),
    )

    # Secondary: 4x3 thematic columns
    build_figure(
        layout_rows=4, layout_cols=3,
        traits_order=ORDER_4X3,
        fname='fig3_4x3',
        combined=combined, sc_colors=sc_colors,
        figsize=(22, 16),
    )

    print('\nDone.')


if __name__ == '__main__':
    main()
