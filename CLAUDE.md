CLAUDE.md — Onboarding for Claude conversations on this project
For teammates: paste the contents of this file at the start of any new Claude conversation about this pipeline. It contains everything Claude needs to be immediately useful instead of re-learning the project from scratch.

1. The 30-second project summary
We are running MAGMA gene-property analysis to map GWAS signals onto the Siletti et al. 2023 human brain cell atlas (461 clusters). The methodology is from Duncan et al. 2025 (Nat Neurosci, DOI 10.1038/s41593-024-01834-w). We are replicating Duncan's pipeline exactly, then applying it to additional GWAS phenotypes our team selects.

Our team's specific goal: resolve splatter neurons (Siletti supercluster #22, a heterogeneous group spanning hypothalamus, midbrain, hindbrain, pons, medulla, mammillary, thalamus) by finding which subclusters within it carry specific GWAS-phenotype enrichment.

Project lead: Nick Monaco. Team: 4 people. All using Claude.

2. Where everything lives
GitHub repo (private): https://github.com/nickmonaco99/magma_celltype_pipeline
Local working copy on lead's Mac: /Users/nickmonaco/Desktop/Science Side Projects/Allen Brain Seattle Trip/magma_celltype_pipeline
Team shared Google Drive: [link to be added by Nick] — holds the Siletti specificity matrix, 1000G reference panel, and snapshot copies
Docker image tag: magma_celltype:latest (built locally on each teammate's machine)
Siletti cluster annotations: auxfiles/science.add7046_table_s3.xlsx (Siletti 2023 Supplementary Table S3, 461 clusters × 22 columns; ~100 KB, committed to repo)

3. Architecture decisions ALREADY MADE — do not re-litigate
Claude should not propose alternatives to these unless someone explicitly asks for a change. They were chosen deliberately to match Duncan's published pipeline.

Decision	Value	Why
MAGMA version	v1.10 (Linux static binary)	Duncan ran v1.10. Watanabe used v1.07. Different versions give slightly different gene-level p-values. We pin v1.10 for exact reproducibility.
Gene annotation window	35 kb upstream, 10 kb downstream	Duncan's choice (matches FUMA default for psychiatric GWAS). Watanabe used 2/1 kb. We use Duncan's.
LD reference panel	1000 Genomes Phase 3 EUR, n=503	Duncan's choice. Watanabe used a UKB-derived reference. We use Duncan's for cross-paper comparability.
Test direction	direction=greater (one-sided positive)	Duncan's choice. We test only positive enrichment per cell type.
Specificity matrix	Siletti_l2_conti_specificity_matrix.txt (precomputed by team)	17,098 genes × 461 clusters. From Chen/Ji/Rong/Jiang preprocessing (see Chen_Ji_Rong_Jiang_Final_Report.pdf in shared materials).
Effective N per SNP	Constant (no per-SNP NEFF column)	UKB-style sumstats provide total N per SNP (NMISS), not effective N. We use constant 4/(1/Ncas + 1/Nctrl) for all SNPs. Loss of precision is negligible at MAF ≥ 1%.
Runtime environment	Docker (magma_celltype:latest)	MAGMA Linux static binary can't execute natively on macOS arm64. All MAGMA calls run inside Docker. Pipeline prep/QC notebooks run on host Python.
Cell-type Bonferroni threshold	P < 0.05 / 461 = 1.085e-4	Number of Siletti L2 clusters tested.
Cluster annotation source	Siletti 2023 Supplementary Table S3 (auxfiles/science.add7046_table_s3.xlsx)	Authoritative per-cluster metadata: name, supercluster (1-31), neurotransmitter, top three anatomical regions, top three dissections, top enriched genes. Joined into outputs by pipeline/04_annotate_results.py. The cluster number alone is biologically meaningless; always look at the annotated CSV.
Plot styling	Siletti supercluster numeric ordering (1-31) with curated rainbow palette	Matches Duncan 2025 Fig 3 conventions. Implemented in pipeline/05_make_plots.py. NOT alphabetical (which would scramble color flow and make adjacent x-positions visually unrelated).

4. Standard prompt to paste at the start of a new Claude conversation
I'm working on the magma_celltype_pipeline project. The repo is at
https://github.com/nickmonaco99/magma_celltype_pipeline. Please read
CLAUDE.md in the repo root for project context.

Today's task: <fill in: "add a new GWAS for X", "interpret cell-type results
for Y", "debug failure in step Z", etc.>

Important context:
- I'm using <macOS arm64 / Linux>
- Docker image magma_celltype:latest is built
- I'm currently in <Cursor / terminal / nothing>
If Claude doesn't have file access to the repo, also paste:

Contents of pipeline/00_config.yaml
Contents of pipeline/README.md (in pipeline/)
The error/output you're trying to interpret

5. Workflow: Adding a new GWAS
This is the most common task. Steps:

Step 1: Identify and download
Find the public sumstats URL. DO NOT use 23andMe-restricted GWAS unless the team has a current DTA. Check the GWAS metadata:

Build (must be hg19 / GRCh37 to match our reference panel; liftover from hg38 if needed)
Ancestry (must be EUR or substantially EUR; mixed-ancestry GWAS need pre-stratification)
N_cases / N_controls (compute N_eff yourself; the paper's reported N may be misleading for case-control)
Download to gwas_sumstats/<short_name>/:

mkdir -p gwas_sumstats/<short_name>
cd gwas_sumstats/<short_name>
curl -L -o <filename> "<url>"

Step 2: Inspect the file FIRST
ALWAYS look at the raw file before trusting any documentation. READMEs lie about column meanings (this happened with insomnia — SNPID_UKB was advertised as the rsID but actually contained chr:pos_allele_allele strings; the real rsID was in an undocumented RSID_UKB column):

zcat <filename> | head -3

Identify which columns hold:

rsID (must look like rs12345..., NOT 1:12345:A:T)
p-value
chromosome, base-pair position
alleles A1, A2
MAF (and which sample it's from, if multiple MAF columns exist — analysis sample > population)
INFO score (if imputed data)

Step 3: Add a GWAS block to 00_config.yaml
Copy the insomnia_ukb: block and modify. Required fields:

name, citation, download_url, license
build, ancestry, n_cases, n_controls, n_effective
raw_filename, raw_subfolder, output_folder
columns: block (verified against raw file inspection)
sentinel_genes: (5 known genes for this trait; thresholds calibrated to the GWAS's N, NOT to the meta-analysis's N — see Section 7)
Set active_gwas: "<your_new_short_name>".

Step 4: Run prep notebook
cd pipeline
python -u 01_prep_sumstats.py 2>&1 | tee /tmp/prep_<gwas>.log

Always use python -u for unbuffered output — otherwise prints don't appear until the script finishes (Python buffers stdout when piped to tee).

Two diagnostic numbers must pass:

λ_GC in [1.1, 1.5] for polygenic brain GWAS (above 1.5 suggests stratification; below 1.1 suggests weak/no signal)
g1000 overlap ≥ 70% (ideally ≥ 95%; below 70% means rsID format mismatch or build mismatch — DO NOT proceed)

Step 5: Run MAGMA stages 1+2 (gene analysis)
cd <project_root>
./run_in_docker.sh bash pipeline/02_run_magma_step1and2.sh

Takes ~10-60 minutes depending on N. Watch for two key MAGMA log lines:

read X lines from file, containing valid SNP p-values for Y SNPs in data (Z% of lines) — Z should be >95%
gene-level analysis ... processed genes: X — should reach 100% without errors

Step 6: Validate gene-level results
python pipeline/02_qc_gene_results.py     # if you've created this file
# OR quick manual check:
awk 'NR>1 && $9 < 2.57e-6' results/<FOLDER>/<file>.step2.genes.out | wc -l

Count Bonferroni-significant genes. Expectation:

500+ → strong (typical of SCZ, BMI)
100-500 → good
30-100 → moderate
10-30 → marginal (cell-type analysis still possible, expect 0-2 hits)
<10 → weak (cell-type analysis likely null)

Step 7: Run MAGMA stage 3 (cell-type analysis)
./run_in_docker.sh bash pipeline/03_run_magma_celltype.sh

Takes ~5 minutes. Output: results/<FOLDER>/Siletti_l2_conti-spe_<FOLDER>.gsa.out.

Step 8: Snapshot the results
python pipeline/snapshot_results.py <gwas_short_name>

Creates results/snapshots/<gwas>_v<N>_<DATE>/ with all key outputs + checksums + provenance manifest. Commit to git, upload to team Drive.

Step 9: Annotate the results (gene symbols + Siletti supercluster + anatomy)
python pipeline/04_annotate_results.py

Reads the latest snapshot (or pass an explicit path). Writes:

results/snapshots/<gwas>_v<N>_<DATE>/annotated/celltype_annotated.csv — every cluster joined with Siletti Table S3 (cluster_name, supercluster, class_auto, neurotransmitter, top_regions, top_dissections, top_genes, n_cells)
results/snapshots/<gwas>_v<N>_<DATE>/annotated/genes_annotated.csv — every gene joined with HGNC symbol from NCBI37.3.gene.loc column 6

Also adds an annotation key to manifest.json and regenerates checksums.sha256. Idempotent — safe to re-run.

Step 10: Generate publication-style plots
python pipeline/05_make_plots.py

Reads the annotated CSVs from the same snapshot. Writes to figures/:

celltype_manhattan.png — all 461 clusters, -log10(P) vs Siletti cluster ID, colored by supercluster in Siletti's 1-31 numeric ordering (matches Duncan 2025 Fig 3 rainbow). Three-tier alpha+size encoding so noise (P > 0.05) fades and signal (P < 0.001) pops. Y-axis floored at 5 so Bonferroni line is always visible.
celltype_forest_top20.png — horizontal forest plot of top 20 clusters, β ± SE error bars, Bonferroni-significant rows bolded red, labels include supercluster ID prefix (e.g., [#15 Upper-layer intratelencephalic]).
genes_manhattan_annotated.png — gene-level Manhattan with GWS-passing gene symbols labeled (500 kb deduplication so clustered loci like the chr17:43Mb MAPT region show just the lead).

Adds a plots key to manifest.json.

6. Schizophrenia validation (the formal pipeline test)
Status as of 2026-05-13: NOT YET RUN. Scheduled for next session.

Why it matters: Until our pipeline reproduces Duncan's exact published SCZ cell-type list, all other GWAS results are PROVISIONAL.

Expected Duncan SCZ independent cell types (from script 5 of duncan_repo/MAGMA/):

Cluster239
Cluster278   ← also appears in our insomnia top 10 (validation pre-signal); = CGE_278
Cluster132   ← also appears in our insomnia top 10 (validation pre-signal); = Misc_132
Cluster233
Cluster423   ← Splat_423, GABA, Amygdala/Thalamus/Hypothalamus (THIS IS A SPLATTER CLUSTER — direct validation that splatter contains psychiatrically-relevant cell types)
Cluster404   ← also appears in our insomnia top 20 (validation pre-signal); = Misc_404

To run:

Download Trubetskoy 2022 PGC3 SCZ wave 3 from figshare: https://figshare.com/articles/dataset/scz2022/19426775
Add schizophrenia_pgc3: block to 00_config.yaml
Set active_gwas: "schizophrenia_pgc3"
Run steps 4-10 above
After Stage 3 produces .gsa.out, write/run 08_validate_against_duncan.py to formally check that our top hits include Cluster239, 278, 132, 233, 423, 404

If validation passes: all subsequent GWAS results trustworthy. If fails: do not proceed with other GWAS until cause identified.

7. Gotchas — things Claude has failed on before
7.1 GWAS file format gotchas
READMEs can lie about column meanings. Always inspect raw data first. The insomnia README said SNPID_UKB was the rsID; it was actually chr:pos_allele_allele format. The real rsID was in RSID_UKB.
Multiple columns with similar names. UKB-style files have both MAF (analysis-sample MAF) and MAF_UKB (UKB-wide MAF). If you rename MAF_UKB→MAF while the file already has MAF, pandas ends up with two columns of the same name and boolean indexing breaks (cannot reindex on an axis with duplicate labels). Defensive code: drop any existing target column before renaming.
Build mismatches. Always verify the GWAS is hg19. If it's hg38, you need liftover before MAGMA — the gene location file and reference panel are hg19.
chr:pos:allele IDs vs rsIDs. MAGMA wants rsIDs. If you accidentally use the synthetic ID column, the pre-flight overlap gate (~0% overlap with g1000_eur) catches it.

7.2 Sentinel-gene thresholds must be GWAS-size-appropriate
Lesson learned from insomnia. Sentinel thresholds set against Watanabe full-meta (N=2.36M) expectations were too aggressive for UKB-only (N=387k). Only 1/5 sentinels passed Bonferroni at gene-level, but the underlying biology was correct (top hits: MAPT/KANSL1 cluster, TCF4, MEIS1, BTBD9 — all known insomnia genes). When picking sentinels for a new GWAS: look up gene-level results from a paper that used a SIMILAR-SIZE GWAS, not the largest available meta.

7.3 macOS / Docker / arm64 quirks
MAGMA Linux static binary doesn't execute on macOS. It's x86_64 Linux ELF. All MAGMA invocations go through Docker. 01_prep_sumstats.py has a graceful fallback that detects this and skips host-side version check.
Docker Desktop sleeps. After Mac sleep / Docker Desktop close, docker run may report "image not found" even though the cached layers are intact. Run docker build -t magma_celltype:latest . to re-register; it'll be near-instant (all CACHED).
Safari strips leading dots from downloaded filenames. .gitignore becomes gitignore in ~/Downloads/. Always rename when moving into the project.

7.4 Python stdout buffering
When piping Python through tee, output buffers until the script finishes. Use python -u for unbuffered output during long-running scripts:

python -u 01_prep_sumstats.py 2>&1 | tee /tmp/prep.log   # CORRECT
python 01_prep_sumstats.py 2>&1 | tee /tmp/prep.log      # WRONG: prints appear at end

7.5 Git / GitHub auth
HTTPS pushes to GitHub require a Personal Access Token, not account password. Generate at github.com/settings/tokens with repo scope.
Pushes of >1 MB can fail with HTTP 400 unless git config --global http.postBuffer 524288000.
SSH is more reliable long-term: ssh-keygen -t ed25519, add to GitHub at github.com/settings/keys, then git remote set-url origin git@github.com:....

7.6 SNP-level proxy ≠ formal Bonferroni
In 01_prep_sumstats.py, the sentinel proxy check looks at the single best SNP within 50kb of each gene. This is a WEAK check that systematically misses genes with distributed (rather than peaked) signal. PTPRD is a classic example: huge gene (2.3 Mb), signal spread across many SNPs, no single hit dominates. The FORMAL gene-level check (in 02_qc_gene_results.py or manual awk) uses MAGMA's LD-aware aggregation and is the authoritative test.

7.7 Cluster numbers ARE cell types (solved — was a gotcha)
Previously, the Siletti specificity matrix gave us column names like Cluster0, Cluster1, ... Cluster460 with no embedded mapping to anatomy or cell-type identity. Looking up biology required digging through Siletti supplementary tables every time.

Now solved: auxfiles/science.add7046_table_s3.xlsx is Siletti's Table S3 (committed to the repo, ~100 KB), and pipeline/04_annotate_results.py merges it into every snapshot. After running annotation, annotated/celltype_annotated.csv tells you everything per cluster: cluster_name (e.g., ULIT_133), supercluster (e.g., Upper-layer intratelencephalic), neurotransmitter (e.g., VGLUT1 VGLUT2), top three regions (e.g., V1C 63% / V2 18% / Pro 15%), top three dissections, top enriched marker genes (e.g., VIPR2, TRPC3, IGFBP2), and cell count.

When working with results, ALWAYS look at the annotated CSV, not the raw .gsa.out. The cluster number alone is biologically meaningless.

Splatter members (supercluster #22): clusters where supercluster == "Splatter" in the annotated CSV. There are 92 — the largest supercluster. They span amygdala, basal forebrain, hypothalamus, thalamus, midbrain, pons, medulla, spinal cord. Includes histaminergic TMN, dopaminergic VTA/SN, serotonergic raphe, cholinergic basal forebrain (the canonical arousal-system neurons). The team's project goal is to find which of these 92 carry trait-specific GWAS enrichment.

7.8 Plot palette ordering must follow Siletti's 1-31 hierarchy
The 31-supercluster palette in 05_make_plots.py is NOT alphabetical. It follows Siletti's published supercluster numbering (Miscellaneous=1, Microglia=2, ..., Splatter=22, ..., Midbrain-derived inhibitory=31). This matches Duncan 2025 Fig 3 and is what makes the rainbow flow visually coherent across the x-axis. If a future Claude proposes "let me sort the legend alphabetically for readability" — push back. The ordering is the convention shared with the published figure, not an arbitrary choice.

8. Current state — as of 2026-05-13

Completed
Local + Docker setup on M3 Pro Mac
All reference data downloaded (1000G EUR, MAGMA v1.10, Siletti matrix, NCBI gene loc, Siletti Table S3 annotations committed to repo)
Pipeline scripts (10 stages total):
  01_prep_sumstats.py — sumstats prep + QC + pre-flight overlap gate
  02_run_magma_step1and2.sh — MAGMA annotation + gene-level
  02_qc_gene_results.py — formal sentinel Bonferroni assertion (file exists, untested)
  03_run_magma_celltype.sh — MAGMA cell-type analysis
  snapshot_results.py — immutable timestamped snapshot + checksums + manifest
  04_annotate_results.py — merges Siletti Table S3 (supercluster, anatomy, markers) + NCBI Entrez → HGNC symbols into snapshot
  05_make_plots.py — cell-type Manhattan (Duncan Fig 3 style with Siletti 1-31 rainbow palette), top-20 forest plot, annotated gene Manhattan

First end-to-end run: insomnia (Watanabe 2022 UKB-only)
  7.73M SNPs after QC
  λ_GC = 1.258 ✓
  g1000 overlap = 98.58% ✓
  26 Bonferroni-significant genes (top: MAPT/KANSL1 chr17 cluster, TCF4, MEIS1)
  1 Bonferroni-significant cell type: Cluster133 = ULIT_133 (Upper-layer intratelencephalic, V1C cortex 63%, VGLUT1/2 glutamatergic; β=24.0, P=6.14e-05). Top markers: VIPR2 (circadian-relevant), TRPC3, IGFBP2.
  Cross-trait validation: 3 of Duncan's 6 SCZ-significant clusters appear in our insomnia top 20: Misc_132 (#2, P=1.42e-4), CGE_278 (#6, P=4.94e-4), Misc_404 (#16, P=1.62e-3). Suggestive cross-trait signal, NOT formal SCZ replication — that requires running the actual SCZ GWAS.
  Top-20 biological pattern: cortical bias (6× CGE interneurons, 3× DLIT, 2× ULIT, 2× Misc) with 2× Mammillary body (#29, sleep/arousal-relevant subcortical) and 1× Upper rhombic lip (#26, cerebellar) as the only non-cortical signals. NO Splatter (#22) clusters in the top 20 — striking given splatter contains the canonical arousal system (histaminergic TMN, dopaminergic VTA, serotonergic raphe). Either insomnia is genuinely cortical, or UKB-only is underpowered for subcortical signals.

Git repository initialized and pushed to private GitHub

Pending — high priority
SCZ validation run (Section 6 above) — formal proof pipeline reproduces Duncan. Special interest: Cluster423 (Splat_423) — if it replicates, that's direct validation that splatter contains psychiatrically-relevant cell types.

Pending — medium priority
Multi-GWAS refactor (the short_name parameterization deferred from initial build — several scripts still hardcode INS_UKB paths)
bootstrap.sh — one-command setup for new teammates
02_qc_gene_results.py formal Bonferroni assertion script (file exists, untested)
Conditional / forward-selection analysis (Duncan scripts 4, 5) — identify independent cell types when multiple are significant
Multi-trait comparison plot in Duncan Fig 3 style (black background, stacked panels per phenotype). Requires ≥2 GWAS runs. Will become the headline figure once SCZ + ≥2 other GWAS are done.
Make project_root portable in 00_config.yaml — currently hardcoded to lead's Mac path

Pending — additional GWAS to consider
Curated list, public, no DTA, splatter-relevant:

Chronic pain (Johnston 2019 UKB) — PAG/brainstem signal
BMI (Yengo 2018 GIANT) — hypothalamus signal
Migraine (FinnGen R12) — trigeminal/thalamic signal
RLS (Akçimen 2024 AoU) — basal ganglia/dopaminergic signal
MDD (PGC MDD2 from figshare) — broad psychiatric signal
PD (Nalls 2019 public, no 23andMe) — midbrain dopamine signal
AD (Bellenguez 2022) — microglia control / negative-control

9. Citation context (for paper writing)
If asked about citations, the correct references are:

Method: Duncan LE et al. 2025. "Linking cell types to genetic risk for schizophrenia and other phenotypes." Nat Neurosci. DOI 10.1038/s41593-024-01834-w.
Siletti atlas (data + Table S3 annotations): Siletti K et al. 2023. "Transcriptomic diversity of cell types across the adult human brain." Science 382, eadd7046.
MAGMA: de Leeuw CA et al. 2015. "MAGMA: Generalized gene-set analysis of GWAS data." PLoS Comput Biol 11(4):e1004219.
Insomnia GWAS: Watanabe K et al. 2022. Nat Genet 54:1125-1132 (UKB component originally from Jansen PR et al. 2019. Nat Genet 51:394-403). DOI 10.1038/s41588-022-01124-w.
1000 Genomes: 1000 Genomes Project Consortium. 2015. Nature 526:68-74.

10. When to escalate to humans (do not over-promise)
Claude should NOT:

Promise specific biological interpretations of single cluster results without supporting evidence beyond the Siletti Table S3 annotations (e.g., "Cluster133 being a V1C neuron means insomnia is fundamentally a visual cortex disease" — V1C is what Siletti annotated, but that's an anatomical label, not a functional mechanism)
Modify the MAGMA invocation parameters (version, window, model direction) without explicit user approval — those match Duncan and changing them invalidates the methodological comparability
Modify the plot palette ordering away from Siletti's 1-31 numeric scheme without explicit user approval — that ordering is the convention shared with the Duncan paper
Skip the SCZ validation step. If a teammate wants to publish insomnia results before SCZ replicates, push back.
Claim results match Duncan's without actual side-by-side comparison

Claude SHOULD:

Be skeptical when results look too clean or too messy
Always recommend snapshotting before any analysis is "done"
Encourage running new GWAS through the whole pipeline (steps 4-10) before drawing cross-trait conclusions
Flag when sentinel-gene thresholds were over-aggressive (we made this mistake on insomnia)
Always cross-reference cluster numbers with annotated/celltype_annotated.csv before discussing biology — the cluster number alone is meaningless

Last updated: 2026-05-13, after annotation + plotting pipeline build (insomnia v1 fully interpreted: Cluster133 = ULIT_133).
