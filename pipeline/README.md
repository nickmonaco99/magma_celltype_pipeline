# MAGMA Cell-Type Pipeline — Local + Docker

Reproduces Duncan et al. 2025 (Nat Neurosci) pipeline on macOS, with bit-identical
Linux behaviour via a minimal `ubuntu:22.04` Docker image. Maps GWAS phenotypes
onto the 461-cluster Siletti adult human brain taxonomy via MAGMA's gene-property
analysis.

**First GWAS**: Watanabe et al. 2022 insomnia, UKB-only component (no DTA needed).

---

## Why this matters for our team

Once this is running, every additional GWAS is ~1 hour of compute and a config
edit. The discipline embedded in this pipeline — version-pinned MAGMA, hard
overlap gate, formal sentinel assertion, provenance logging — is what lets the
splatter-neurons writeup withstand reviewer scrutiny without us scrambling.

---

## Deviations from Watanabe 2022 (intentional)

We use **Duncan's parameters**, not Watanabe's, because the pipeline is what
Duncan validated at 461-cluster resolution. These deviations are documented
here so the methods section is clean:

| Parameter         | Watanabe 2022 | This pipeline | Why                                |
|-------------------|---------------|---------------|------------------------------------|
| MAGMA version     | v1.07         | **v1.10**     | Match Duncan exactly               |
| Annotation window | 2 / 1 kb      | **35 / 10 kb**| Match Duncan (Bryois/Skene convention) |
| Reference panel   | 1000G EUR     | 1000G EUR     | Same                               |
| Gene set source   | (no cell-type analysis in original) | Siletti 461 clusters | This is the whole point |

The version + window mismatch with the original Watanabe analysis means our
gene-level p-values won't be numerically identical to theirs. That's expected.
What we DO match: the gene-level top hits (PTPRD, MEIS1, BTBD9, LSAMP, NEGR1
all clear Bonferroni). This is enforced in `02_qc_gene_results.py`.

---

## Pipeline stages

```
01_prep_sumstats.py
    ├── MAGMA version assertion (must be v1.10)
    ├── Inspect + munge Watanabe sumstats to Duncan format
    ├── λ_GC + Manhattan + QQ
    ├── SNP-level proxy check on Jansen-2019 sentinels (weak gate)
    └── PRE-FLIGHT OVERLAP GATE (≥70% of our SNPs in g1000_eur.bim) ← hard error if fails

02_run_magma_step1and2.sh        [in Docker]
    └── Stages 1+2: annotation (35/10 kb) + gene-level analysis
        Writes run.log with command, version, timestamp

02_qc_gene_results.py
    ├── Print top-30 gene hits (eyeball vs Jansen 2019 Supp Table 8)
    └── FORMAL ASSERTION: 5/5 sentinels <2.5e-6, ≥3/5 <1e-8 ← hard error if fails

03_run_magma_celltype.sh         [in Docker]
    └── Stage 3: cell-type analysis with Siletti specificity matrix

04_celltype_results.py
    └── Parse .gsa.out, identify Bonferroni-significant cell types,
        map Cluster# → splatter / region / supercluster

[then conditional analysis stages 05+ once we have hits]
```

---

## One-time setup

```bash
# 1. Run the bootstrap (creates dirs, moves Duncan's repo, checks Docker)
cd "/Users/nickmonaco/Desktop/Science Side Projects/Allen Brain Seattle Trip"
bash magma_celltype_pipeline/setup.sh

# 2. Manually download four things (URLs in setup.sh output):
#    - MAGMA v1.10 LINUX binary       → magma_binary/magma
#    - g1000_eur reference panel     → auxfiles/g1000_eur.{bed,bim,fam,synonyms}
#    - Siletti specificity matrix    → gene-level/Siletti_l2_conti_specificity_matrix.txt
#    - Watanabe insomnia UKB GWAS    → gwas_sumstats/insomnia_watanabe2022_ukb/

# 3. Build the Docker image (uses your downloaded magma_binary/magma)
cd magma_celltype_pipeline
docker build -t magma_celltype:latest .
./run_in_docker.sh /usr/local/bin/magma --version    # smoke test → should print v1.10

# 4. Host Python (for running the notebooks in Cursor)
conda create -n magma_celltype python=3.11 -y
conda activate magma_celltype
pip install pandas numpy scipy matplotlib pyyaml jupyterlab ipykernel
python -m ipykernel install --user --name magma_celltype
```

---

## Running the pipeline

Open Cursor at the project root. Then for each step:

```bash
# Step 1: prep sumstats (host Python — fast iteration, plots, inspection)
# Open pipeline/01_prep_sumstats.py — run cells with Cmd+Enter
# First pass: run cells 1-4 to inspect file
# Update 00_config.yaml columns block
# Re-run from cell 1 to end. Stops with hard error if any gate fails.

# Step 2: run MAGMA stages 1+2 (in Docker for Linux parity)
./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh
# Takes ~45-60 min. Writes run.log to results/INS_UKB/

# Step 3: validate gene-level output (host Python)
# Open pipeline/02_qc_gene_results.py — assertion raises if sentinels fail

# Step 4: run cell-type analysis (in Docker)
./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh
# ~5 min. Output: results/INS_UKB/Siletti_l2_conti-spe_INS_UKB.gsa.out

# Step 5: interpret cell-type results (host Python)
# Open pipeline/04_celltype_results.py
```

---

## Folder layout

```
magma_celltype_pipeline/
├── duncan_repo/                         (Duncan's downloaded repo, unmodified)
├── magma_binary/magma                   (v1.10 Linux binary)
├── Dockerfile                           (ubuntu:22.04 + MAGMA + Python + R)
├── run_in_docker.sh                     (wrapper to exec commands in container)
├── setup.sh                             (one-time bootstrap)
├── auxfiles/                            (g1000_eur.{bed,bim,fam,synonyms} + gene.loc symlink)
├── gene-level/                          (Siletti specificity matrix)
├── gwas_sumstats/                       (downloaded GWAS files)
├── pipeline/
│   ├── 00_config.yaml                   (single source of truth)
│   ├── 01_prep_sumstats.py              (host Python: prep + gates)
│   ├── 02_run_magma_step1and2.sh        (Docker bash: MAGMA stages 1+2)
│   ├── 02_qc_gene_results.py            (host Python: formal assertion)
│   ├── 03_run_magma_celltype.sh         (Docker bash: MAGMA stage 3)
│   ├── 04_celltype_results.py           (host Python: interpretation)
│   ├── 05_run_conditional.R             (Docker R: Duncan script 3)
│   ├── 06_run_pairwise.sh               (Docker bash: Duncan script 4)
│   ├── 07_forward_selection.R           (Docker R: Duncan script 5)
│   └── README.md
├── results/INS_UKB/
│   ├── snploc_insomnia_ukb              (MAGMA inputs we generate)
│   ├── insomnia_ukb.no_heading
│   ├── *.step1.genes.annot              (MAGMA outputs)
│   ├── *.step2.genes.raw
│   ├── *.step2.genes.out
│   ├── Siletti_l2_conti-spe_INS_UKB.gsa.out
│   ├── top30_gene_level_hits.tsv
│   └── run.log                          (provenance: every MAGMA cmd + version + ts)
└── figures/                             (Manhattan, QQ, cell-type plots)
```

---

## Hard gates (where the pipeline halts and refuses to proceed)

| Gate                            | Location                | Threshold                           |
|---------------------------------|-------------------------|-------------------------------------|
| MAGMA version pin               | `01_prep_sumstats.py`   | Must contain `v1.10`                |
| 1000G reference panel overlap   | `01_prep_sumstats.py`   | ≥70% of our SNPs in g1000_eur.bim  |
| Gene-level sentinels (Bonferroni)| `02_qc_gene_results.py` | All 5/5 below 2.5e-6                |
| Gene-level sentinels (strict)   | `02_qc_gene_results.py` | At least 3/5 below 1e-8             |

Pre-flight checks before MAGMA runs ⇒ no wasted compute on misconfigured runs.

---

## Provenance logging

Every MAGMA invocation appends to `results/INS_UKB/run.log`:
- Timestamp (ISO 8601 UTC)
- MAGMA version string
- Exact command line with all flags
- Path of input files (with sizes)
- Path of output files (with sizes)
- Pass/fail status of any gates

This file is the audit trail. If a reviewer asks "what exact command produced
this number," it's all here.

---

## Working in Cursor — tips for this project

- **Open the project root in Cursor** (`cursor magma_celltype_pipeline/`).
  The `.py` files with `# %%` markers auto-detect as notebook cells in Cursor.
- **Cmd+Enter** runs the current cell. **Shift+Enter** runs and moves to next.
- **Cmd+J** toggles the terminal panel. Keep it open and use it for
  `./run_in_docker.sh` invocations.
- **Cmd+L** opens the AI chat panel; useful for explaining errors or quick
  syntax help. Less useful for statistical interpretation — paste output to
  Claude here instead, where the full project context is loaded.
- **Cmd+K** does inline edits (highlight code, ask for a change). Good for
  refactoring; risky for "rewrite this function from scratch" since Cursor
  doesn't know about MAGMA's quirks.
- Cursor's autocomplete will sometimes hallucinate pandas API. If something
  looks off, check official docs or ask Claude.

---

## Citations

- Pipeline: Duncan LE et al. 2025. *Nat Neurosci* 28: [in press]. DOI: 10.1038/s41593-024-01834-w
- GWAS: Watanabe K et al. 2022. *Nat Genet* 54:1125-1132. DOI: 10.1038/s41588-022-01124-w
- Sentinel anchor: Jansen PR et al. 2019. *Nat Genet* 51:394-403. DOI: 10.1038/s41588-018-0333-3
- Cell atlas: Siletti K et al. 2023. *Science* 382:eadd7046. DOI: 10.1126/science.add7046
- MAGMA: de Leeuw CA et al. 2015. *PLoS Comput Biol* 11:e1004219. DOI: 10.1371/journal.pcbi.1004219
- Specificity method: Bryois J et al. 2020. *Nat Genet* 52:482-493. DOI: 10.1038/s41588-020-0610-9
