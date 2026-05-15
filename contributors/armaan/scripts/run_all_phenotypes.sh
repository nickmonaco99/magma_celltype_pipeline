#!/bin/bash


MAGMA=/mnt/c/PhD_Data/allen_transcriptomics/migraine_celltypes/magma/magma
GENELOC=/mnt/c/PhD_Data/allen_transcriptomics/migraine_celltypes/magma/NCBI37.3.gene.loc
BFILE=/mnt/c/PhD_Data/allen_transcriptomics/migraine_celltypes/magma/aux/g1000_eur
COVAR=/mnt/c/PhD_Data/allen_transcriptomics/migraine_celltypes/data/Siletti_l2_conti_specificity_matrix.txt
DATA=/mnt/c/PhD_Data/allen_transcriptomics/migraine_celltypes/data
RESULTS=/mnt/c/PhD_Data/allen_transcriptomics/migraine_celltypes/results

declare -A CODES
CODES[sleep_apnea]=G6_SLEEPAPNO
CODES[sleep_combined]=SLEEP
CODES[menst_excessive]=N14_MESNRUIRREG
CODES[menorrhagia]=N14_MENORRHAGIA
CODES[menst_irregular]=N14_MENSIRREG

declare -A NVALS
NVALS[sleep_apnea]=453733
NVALS[sleep_combined]=453733
NVALS[menst_excessive]=293218
NVALS[menorrhagia]=293218
NVALS[menst_irregular]=293218

for PHENO in sleep_apnea sleep_combined menst_excessive menorrhagia menst_irregular; do
    CODE=${CODES[$PHENO]}
    N=${NVALS[$PHENO]}
    GWAS=${DATA}/finngen_R12_${CODE}_magma.txt
    OUT=${RESULTS}/${PHENO}

    echo "================================================"
    echo "Starting: ${PHENO} (${CODE}) N=${N}"
    echo "$(date)"
    echo "================================================"

    echo "[Step 1] Annotating..."
    ${MAGMA} --annotate window=35,10 \
             --snp-loc ${GWAS} \
             --gene-loc ${GENELOC} \
             --out ${OUT}_annotation

    echo "[Step 2] Gene-based analysis..."
    ${MAGMA} --bfile ${BFILE} \
             --pval ${GWAS} use=SNP,P N=${N} \
             --gene-annot ${OUT}_annotation.genes.annot \
             --out ${OUT}_genes

    echo "[Step 3] Cell type enrichment..."
    ${MAGMA} --gene-results ${OUT}_genes.genes.raw \
             --gene-covar ${COVAR} \
             --model direction=greater \
             --out ${OUT}_celltypes

    echo "DONE: ${PHENO} at $(date)"
done

echo "ALL PHENOTYPES COMPLETE"