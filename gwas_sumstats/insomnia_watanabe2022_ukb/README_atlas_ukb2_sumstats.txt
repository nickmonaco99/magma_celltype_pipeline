# GWAS summary statistics for UK Biobank release 2
# README file initialized on 21 Sep 2018
# updated on 14 Jan 2019

# File name
This archive file contains 600 GWAS summary statistics.
Each file contains phenotype ID from UK Biobank and "res" or "logistic" for linear or logistic analyses, respectively.

# Columns
SNP: unique ID of the SNP consists of chromosome, position and alphabetically ordered alleles
CHR: chromosome
BP: base pair position on GRCh37
A1: effect allele
TEST: Type of test (ADD for all files)
NMISS: Number of non-missing genotypes
BETA/OR: Regression coefficient or odds ratio
SE: Standard error (for binary traits, SE in logOR scale)
L95: Lower bound on confidence interval for CMH odds ratio
U95: Upper bound on confidence interval for CMH odds ratio
STAT: Coefficient t-statistics
P: P-value
A2: non effect allele
MAF: Minor allele frequency
NCHROBS: Number of allele observation
SNPID_UKB: rsID provided by UK Biobank
A1_UKB: A1 allele in UK Biobank
A2_UKB: A2 allele in UK Biobank
INFO_UKB: Info score provided by UK Biobank
MAF_UKB: MAF of entire UK Boiobank samples

# Phenotype information
gwasDB_ukb2_sumstats.xlsx: The file is the database used for the atlas of GWAS summary statistics website but only 600 UKB summary staitistcs are extracted.
ukb2_gwas_pheno_info.xlsx: The file contains information of phenotype definition with exception. Details can be found in the publication listed below.

# Phenotypes with highly skewed distribusion and long tail
There are multiple phenotypes with highlye skewed distribusion and long tail which resulted in inflation of summary statistics. Those phenotype need a caution and advised to exclude for any systematic analysis.
Atlas ID for those phenotypes:
3193
3198
3199
3239
3243
3266
3275
3306
3316
3326
3381
3382
3383
3393
3408
3623
3723
3738

# Citation
Watanabe et al. A global view of pleiotropy and genetic architecture in complex traits. bioRxiv 2018
