magma_celltype_pipeline
Apply MAGMA gene-property analysis to GWAS summary statistics, mapping genetic associations onto the Siletti et al. 2023 human brain cell atlas (461 clusters). Replicates the methodology from Duncan et al. 2025 (Nat Neurosci, DOI 10.1038/s41593-024-01834-w).

Team-specific goal: resolve the heterogeneous splatter neurons (Siletti supercluster #22, mixed hypothalamus / midbrain / hindbrain / pons / medulla / thalamus) by mapping per-phenotype GWAS signal onto sub-clusters within it.

For Claude users: see CLAUDE.md before starting a new conversation.

Status
GWAS	N_eff	λ_GC	Bonferroni genes	Bonferroni cell types	Snapshot
Insomnia (Watanabe 2022 UKB-only)	314,102	1.258	26	1 (Cluster133 = ULIT_133)	results/snapshots/insomnia_ukb_v1_2026-05-13/
Schizophrenia (Trubetskoy 2022)	—	—	—	—	pending — validation run

Prerequisites
macOS (Apple Silicon or Intel) or Linux
Docker Desktop (or Docker Engine on Linux)
Python 3.11+ with: pandas numpy scipy matplotlib pyyaml openpyxl (conda or pip)
Git + a GitHub account with access to this repo
~10 GB free disk space for reference data + working files

Quick start (for new teammates)
1. Clone
git clone git@github.com:nickmonaco99/magma_celltype_pipeline.git
cd magma_celltype_pipeline

2. Get reference data
Some files are too large for git and must be downloaded separately. Either:

Option A: Use bootstrap.sh (TODO: not yet written)

./bootstrap.sh

Option B: Manual (for now)

Download these into the locations shown:

File	Source	Destination	Size
MAGMA v1.10 Linux static binary	https://ctg.cncr.nl/software/magma	magma_binary/magma	7 MB
1000G EUR reference panel (4 files)	https://ctg.cncr.nl/software/magma	auxfiles/g1000_eur.{bed,bim,fam,synonyms}	~3.5 GB
Siletti specificity matrix	Team Google Drive	gene-level/Siletti_l2_conti_specificity_matrix.txt	148 MB
(NCBI gene location file)	Symlinked from duncan_repo/Data/	auxfiles/NCBI37.3.gene.loc	tiny

Note: the Siletti cluster annotations (auxfiles/science.add7046_table_s3.xlsx) and Duncan reference code (duncan_repo/) are committed to the repo and need no separate download.

Then:

# Create the symlink for the gene location file
ln -s ../duncan_repo/Data/NCBI37.3.gene.loc auxfiles/NCBI37.3.gene.loc

# Verify
ls -la auxfiles/

3. Edit the project_root path in 00_config.yaml
The config currently hardcodes the project root to the lead's Mac:

project_root: "/Users/nickmonaco/Desktop/Science Side Projects/Allen Brain Seattle Trip/magma_celltype_pipeline"

Edit this line to your own absolute path. (This will be removed in a future refactor — scripts 04 and 05 already use dynamic path resolution.)

4. Build the Docker image
docker build -t magma_celltype:latest .

Takes ~5 min the first time. Validates that MAGMA v1.10 is correctly installed.

5. Test it works
./run_in_docker.sh /usr/local/bin/magma --version

Should print Welcome to MAGMA v1.10 (linux/s). If yes, you're set up.

6. Run the existing insomnia example end-to-end (optional sanity check)
Skip to "How to add a new GWAS" below if you have a phenotype to run. To re-run the insomnia example you'd need to download its sumstats too:

mkdir -p gwas_sumstats/insomnia_watanabe2022_ukb
cd gwas_sumstats/insomnia_watanabe2022_ukb
curl -L -o insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz \
  "https://vu.data.surfsara.nl/index.php/s/06RsHECyWqlBRwq/download?path=&files=insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz"
cd ../..
python -u pipeline/01_prep_sumstats.py 2>&1 | tee /tmp/01_prep.log
./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh
./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh
python pipeline/snapshot_results.py insomnia_ukb
python pipeline/04_annotate_results.py
python pipeline/05_make_plots.py

You should reproduce the v1 snapshot bit-for-bit (compare checksums).

Directory structure
magma_celltype_pipeline/
├── CLAUDE.md                   # ← read this before starting a Claude session
├── README.md                   # ← this file
├── Dockerfile                  # Ubuntu 22.04 + MAGMA v1.10 + Python + R
├── run_in_docker.sh            # wrapper to invoke commands inside Docker
├── setup.sh                    # legacy bootstrap (being replaced)
├── .gitignore                  # excludes large data files
│
├── pipeline/                   # all our code
│   ├── 00_config.yaml          # single source of truth for paths + GWAS metadata
│   ├── 01_prep_sumstats.py     # sumstats → MAGMA-ready, with QC and pre-flight gate
│   ├── 02_qc_gene_results.py   # formal sentinel-gene Bonferroni assertion
│   ├── 02_run_magma_step1and2.sh   # MAGMA annotation + gene-level (Duncan stages 1+2)
│   ├── 03_run_magma_celltype.sh    # MAGMA gene-property cell-type analysis (Duncan stage 3)
│   ├── snapshot_results.py     # immutable timestamped snapshot + checksums + manifest
│   ├── 04_annotate_results.py  # merges Siletti Table S3 + NCBI symbols into snapshot
│   ├── 05_make_plots.py        # Duncan-style cell-type Manhattan + forest + annotated gene Manhattan
│   └── README.md               # script-level docs
│
├── duncan_repo/                # vendored copy of Duncan's published code
│   ├── MAGMA/                  # his 5 scripts (we mirror 1+2+3 in pipeline/)
│   ├── Example_results/        # his SCZ outputs for comparison
│   ├── Preprocessing_Siletti/  # his code for generating the specificity matrix
│   └── Data/NCBI37.3.gene.loc  # the canonical gene location file (hg19)
│
├── auxfiles/                   # reference data
│   ├── science.add7046_table_s3.xlsx   # Siletti 2023 Table S3 cluster annotations (committed, ~100 KB)
│   ├── g1000_eur.{bed,bim,fam,synonyms}  # 1000G EUR LD reference panel (gitignored, ~3.5 GB)
│   └── NCBI37.3.gene.loc       # symlink to duncan_repo/Data/...
│
├── magma_binary/               # MAGMA Linux static binary (gitignored)
│   └── magma                   # v1.10
│
├── gene-level/                 # Siletti specificity matrix (gitignored)
│   └── Siletti_l2_conti_specificity_matrix.txt  # 17098 genes × 461 clusters
│
├── gwas_sumstats/              # raw GWAS files (gitignored, downloaded per-GWAS)
│   └── insomnia_watanabe2022_ukb/
│       ├── README_atlas_ukb2_sumstats.txt
│       └── insomnia_ukb2b_EUR_sumstats_...txt.gz
│
└── results/                    # MAGMA outputs and snapshots
    ├── INS_UKB/                # working files (gitignored)
    │   ├── insomnia_ukb.no_heading             # prepared sumstats
    │   ├── snploc_insomnia_ukb                 # SNP locations
    │   ├── *.step1.genes.annot                 # SNP→gene annotations
    │   ├── *.step2.genes.out                   # gene-level p-values
    │   ├── *.step2.genes.raw                   # gene-level raw stats
    │   ├── Siletti_l2_conti-spe_INS_UKB.gsa.out  # cell-type results
    │   └── run.log                             # provenance log
    └── snapshots/              # tracked in git, immutable per-run results
        └── insomnia_ukb_v1_2026-05-13/
            ├── manifest.json   # includes 'annotation' and 'plots' keys after steps 9+10
            ├── checksums.sha256
            ├── 00_config.yaml
            ├── run.log
            ├── gene_level/
            ├── cell_type/
            ├── annotated/                          # produced by 04_annotate_results.py
            │   ├── celltype_annotated.csv          # 461 clusters × Siletti metadata
            │   └── genes_annotated.csv             # 18,381 genes with HGNC symbols
            └── figures/
                ├── manhattan.png                   # original gene-level (from snapshot_results.py)
                ├── qqplot.png                      # original (from snapshot_results.py)
                ├── celltype_manhattan.png          # produced by 05_make_plots.py
                ├── celltype_forest_top20.png       # produced by 05_make_plots.py
                └── genes_manhattan_annotated.png   # produced by 05_make_plots.py

How to add a new GWAS
The end-to-end workflow (~1-2 hours of compute, mostly unattended):

Step 1: Find a public GWAS
Criteria:

Build hg19 (or liftover from hg38 yourself before continuing)
Predominantly European ancestry (or use a multi-ancestry GWAS's EUR-only stratum)
Public release (no DTA-restricted 23andMe-only data)
Reasonable sample size (effective N ≥ 50,000; less is borderline)
Pre-vetted candidates we've discussed:

Trait	Source	Public?	Splatter relevance
Schizophrenia	Trubetskoy 2022 (PGC3)	Yes, figshare	Pipeline validation target
BMI	Yengo 2018 GIANT	Yes	Hypothalamic
Chronic pain	Johnston 2019 UKB	Yes	PAG / brainstem
Migraine	FinnGen R12 G6_MIGRAINE	Yes, finngen.fi	Trigeminal / thalamic
RLS	Akçimen 2024 AoU	Yes	Basal ganglia / dopaminergic
MDD	PGC MDD2 figshare	Yes	Broad psychiatric
PD	Nalls 2019 (public, no 23andMe)	Yes	Midbrain dopamine
AD	Bellenguez 2022	Yes	Microglia control

Step 2: Download
mkdir -p gwas_sumstats/<short_name>
cd gwas_sumstats/<short_name>
curl -L -o <filename> "<url>"

Step 3: Inspect the raw file
ALWAYS look at the actual data before trusting any README. See CLAUDE.md Section 7.1 for the gotcha that bit us on insomnia (SNPID_UKB claimed to be the rsID but wasn't).

zcat gwas_sumstats/<short_name>/<filename> | head -3

Identify which columns hold: rsID, p-value, chromosome, position, alleles, MAF, INFO.

Step 4: Add a 00_config.yaml block
Copy the existing insomnia_ukb: block as a template. Fill in citation, N, column mapping, sentinel genes.

⚠️ Sentinel gene thresholds should reflect what you'd expect at the GWAS's ACTUAL sample size, not the largest published meta-analysis. See CLAUDE.md Section 7.2.

Set active_gwas: "<your_short_name>".

Step 5: Run the prep notebook
cd pipeline
python -u 01_prep_sumstats.py 2>&1 | tee /tmp/prep_<short_name>.log

Must pass:

λ_GC ∈ [1.1, 1.5] (or whatever's appropriate for the trait)
g1000 overlap ≥ 70% (the pre-flight gate hard-fails below this)
If either fails, stop and debug before continuing.

Step 6: Run MAGMA stages 1+2
cd ..
./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh

Takes 10-60 min depending on SNP count.

Step 7: Run MAGMA stage 3 (cell-type analysis)
./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh

Takes ~5 min.

Step 8: Snapshot
python pipeline/snapshot_results.py <short_name>

Creates results/snapshots/<short_name>_v1_<DATE>/. Commit to git, upload to team Drive.

Step 9: Annotate
python pipeline/04_annotate_results.py

Adds annotated/celltype_annotated.csv (with Siletti supercluster, anatomy, top markers) and annotated/genes_annotated.csv (with HGNC symbols) to the snapshot. Idempotent — safe to re-run.

Step 10: Generate plots
python pipeline/05_make_plots.py

Writes three publication-style figures to the snapshot's figures/:
- celltype_manhattan.png — 461-cluster Manhattan, Duncan/Siletti rainbow palette
- celltype_forest_top20.png — top-20 forest plot with β ± SE
- genes_manhattan_annotated.png — gene Manhattan with HGNC symbols labeled

Step 11: Interpret + share with team
Open results/snapshots/<short_name>_v1_<DATE>/annotated/celltype_annotated.csv to identify what each top cluster IS (cluster_name, supercluster, anatomical regions, top marker genes). Flag any clusters in supercluster #22 (splatter == "Splatter"). Open the three figures for visual summary. Commit the snapshot to git.

Pipeline validation status
Important: until our SCZ run reproduces Duncan's exact published cell-type list, all other GWAS results are PROVISIONAL.

Duncan's published SCZ independent cell types (the validation target):

Cluster239, Cluster278, Cluster132, Cluster233, Cluster423, Cluster404

Already in our insomnia top 20 (cross-trait pre-validation signal, with cluster identities resolved via Siletti Table S3):

Cluster132 = Misc_132 at #2 (P=1.42e-4)
Cluster278 = CGE_278 at #6 (P=4.94e-4)
Cluster404 = Misc_404 at #16 (P=1.62e-3)

Insomnia's top Bonferroni hit, Cluster133 = ULIT_133 (Upper-layer intratelencephalic, primarily V1C visual cortex, VGLUT1/2 glutamatergic), is NOT in Duncan's SCZ list — appropriately, since insomnia and SCZ have different cell-type heritability profiles.

The SCZ run will be the formal pipeline check.

Team
Nick Monaco — project lead, pipeline builder, splatter resolution
Teammate 2
Teammate 3
Teammate 4

All four members use Claude for development. See CLAUDE.md for the standard prompt to paste at the start of new conversations.

Citations
If you publish work using this pipeline, cite:

Method: Duncan LE et al. 2025. Nat Neurosci. DOI: 10.1038/s41593-024-01834-w
Siletti atlas (data + Table S3): Siletti K et al. 2023. Science 382, eadd7046.
MAGMA: de Leeuw CA et al. 2015. PLoS Comput Biol 11(4):e1004219.
1000G: 1000 Genomes Project Consortium. 2015. Nature 526:68-74.
The GWAS-specific paper for each phenotype you analyze (see 00_config.yaml citation: field per GWAS).

License
Pipeline code: MIT (or similar — team to decide).

Data licenses:

GWAS sumstats from Watanabe 2022 / Jansen 2019: CC BY-NC-SA 4.0 (non-commercial)
Siletti 2023 atlas data (incl. Table S3): per Siletti et al. publication
1000G data: open
MAGMA: free for academic use, per MAGMA license

This repository is private. Do not share without coordinating with the team about data redistribution permissions.

Last meaningful update: 2026-05-13 (annotation + plotting pipeline; Cluster133 identified as ULIT_133)
