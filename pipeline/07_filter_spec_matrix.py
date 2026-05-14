"""Filter Siletti specificity matrix to Bonferroni-significant clusters from marginal MAGMA output.

Equivalent to creating Duncan's `Siletti_l2_specificity_matrix_sig-only` file.
Used as preprocessing before MAGMA --model joint-pairs (Duncan 2025 Stage 4).
"""
import argparse
import pandas as pd

N_CLUSTERS = 461
BONF = 0.05 / N_CLUSTERS  # 1.085e-4

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--marginal-gsa', required=True,
                help='Path to marginal MAGMA .gsa.out file')
ap.add_argument('--spec-matrix', required=True,
                help='Path to full Siletti specificity matrix (gene-level/Siletti_l2_conti_specificity_matrix.txt)')
ap.add_argument('--output', required=True,
                help='Path for filtered specificity matrix output')
args = ap.parse_args()

# Identify Bonferroni-significant clusters from marginal results
marg = pd.read_csv(args.marginal_gsa, sep=r'\s+', comment='#')
sig = marg.loc[marg['P'].astype(float) < BONF, 'VARIABLE'].tolist()
print(f"  {len(sig)} Bonferroni-significant clusters in marginal (P < {BONF:.2e})")

if len(sig) < 2:
    raise SystemExit(f"  Cannot run conditional analysis with only {len(sig)} significant clusters")

# Read full specificity matrix (tab-separated, first column is gene name, rest are clusters)
spec = pd.read_csv(args.spec_matrix, sep=r'\s+')

# Keep gene-name column (first) + the significant cluster columns
gene_col = spec.columns[0]
keep = [gene_col] + [c for c in spec.columns[1:] if c in sig]
missing = [c for c in sig if c not in spec.columns]
if missing:
    print(f"  WARNING: {len(missing)} significant clusters absent from spec matrix: {missing[:5]}...")

spec[keep].to_csv(args.output, sep='\t', index=False)
print(f"  Wrote {len(keep)-1} cluster columns × {len(spec)} genes -> {args.output}")
