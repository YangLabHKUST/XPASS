# XPASS
The XPASS package implement the XPASS approach for generating PRS in a target population by integrating multi-ethnic datasets.

# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/XPASS")
```

# Quick start

We illustrate the usage of XPASS using the GWAS summary statistics of height from UKB and BBJ. For convenience, we use the 1000 Genomes project genotypes as reference panels, which may not achieve optimal prediction accuracy due to the limited sample size. In practice, it is suggested to use larger datasets as reference panels (n>2000). The datasets involved in the following example can be downloaded from [here](https://www.dropbox.com/sh/i7rhnko69974dje/AACfcDXz0cmwshbli8q7PZA5a?dl=0).

The GWAS summary statistics should be prepared in the following format:
```
head height_bbj_3M_format.txt

SNP	N	Z	A1	A2
rs117086422	159095	-1.20413423957762	T	C
rs28612348	159095	-1.25827042089233	T	C
rs4475691	159095	-1.19842287777303	T	C
rs950122	159095	-1.2014434974188	C	G
rs3905286	159095	-1.27046106136441	T	C
rs28407778	159095	-1.26746342605063	A	G
rs4246505	159095	-1.24706211285128	A	G
rs4626817	159095	-1.26366297625074	A	G
rs11507767	159095	-1.28611566069053	G	A
```

```
head height_ukb_3M_format.txt

SNP	N	Z	A1	A2
rs117086422	429312	1.42436004338939	T	C
rs28612348	429312	1.48706291417224	T	C
rs4475691	429312	1.53977372135067	T	C
rs950122	429312	1.37958155329171	C	G
rs3905286	429312	1.77045946243262	T	C
rs28407778	429312	1.9908370573435	A	G
rs4246505	429312	1.90922505355565	A	G
rs4626817	429312	1.53216668392479	A	G
rs11507767	429312	1.55873328059033	G	A
```
Once the summary files are formatted, XPASS will automatically process the datasets, including SNPs overlapping and allele matching.
Run XPASS with the following comand:
```{r}
# library(devtools)
# install_github("https://github.com/YangLabHKUST/XPASS")
library(XPASS)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(30)

# reference genotypes for EAS (prefix of plink file bim/bed/fam)
ref_EAS <- "1000G.EAS.QC.hm3.ind"

# covariates of EAS reference genotypes
cov_EAS <- "1000G.EAS.QC.hm3.ind.pc5.txt"

# reference genotypes for EUR (prefix of plink file bim/bed/fam)
ref_EUR <- "1000G.EUR.QC.hm3.ind"

# covariates of EUR reference genotypes
cov_EUR <- "1000G.EUR.QC.hm3.ind.pc20.txt"

# genotype file of test data (plink prefix). 
# Note: for demonstration, we assume that the genotypes of prediction target are used as the reference panel of target population. 
# In practice, one can also use genotypes from other sources as reference panel.
height_test <- "1000G.EAS.QC.hm3.ind"


# sumstats of height
height_bbj <- "height_bbj_3M_format.txt"  # target
height_ukb <- "height_ukb_3M_format.txt"  # auxiliary

fit_bbj <-XPASS(file_z1 = height_bbj,file_z2 = height_ukb,file_ref1 = ref_EAS,file_ref2 = ref_EUR,
                file_cov1 = cov_EAS,file_cov2 = cov_EUR,file_predGeno = height_test,
                pop = "EAS",sd_method="LD_block",compPosMean = T,
                file_out = "height_bbj_ukb_ref_TGP")

Summary statistics file 1: height_bbj_3M_format.txt
Summary statistics file 2: height_ukb_3M_format.txt
Reference file 1: 1000G.EAS.QC.hm3.ind
Reference file 2: 1000G.EUR.QC.hm3.ind
Covariates file 1: 1000G.EAS.QC.hm3.ind.pc5.txt
Covariates file 2: 1000G.EUR.QC.hm3.ind.pc20.txt
Test genotype file : 1000G.EAS.QC.hm3.ind
Reading data from summary statisitcs...
3621503 and 3621503 SNPs found in summary statistics files 1 and 2.
There are two reference panels. Assume two phenotypes are from different populations.
Reading SNP info from reference panels...
1209411 and 1313833 SNPs found in reference panel 1 and 2.
Reading SNP info from test genotype file...
1209411 SNPs found in test file.
754616 SNPs are matched in all files.
0 SNPs are removed because of ambiguity; 754616 SNPs remained.
Calculating kinship matrix from the both reference panels...
128169 SNPs in the second reference panel are alligned for alleles according to the first.
14410 SNPs have different minor alleles in phenotype 1, z-scores are corrected according to reference panel.
14410 SNPs have different minor alleles in phenotype 1, z-scores are corrected according to reference panel.
Assigning SNPs to LD Blocks...
Calculate PVE...
Compute posterior mean from  1443  blocks ...
...
Predicting PRS from test genotypes...
Done.
             h1         h2        h12        rho
[1,] 0.43474307 0.63559207 0.37876061 0.72054189
[2,] 0.02087136 0.03403142 0.02031872 0.01727539
```

XPASS returns a list of results:
```{r}
# H: a table of estimated heritabilities, co-heritability and genetic correlation (first row)
# and their corresponding standard erros (second row).
> fit_bbj$H
             h1         h2        h12        rho
[1,] 0.43474307 0.63559207 0.37876061 0.72054189
[2,] 0.02087136 0.03403142 0.02031872 0.01727539

# mu: a data frame storing the posterior means computed by LDpred-inf using only the target dataset (mu1) and
# only the auxiliary dataset (mu2), and the posterior mean computed by XPASS (mu_XPASS).
> head(fit_bbj$mu)
        SNP           mu1           mu2      mu_XPASS
1 rs4475691 -0.0008581511  0.0004678541 -0.0004164065
2 rs7537756 -0.0017717779  0.0010456091 -0.0004210284
3 rs7523549  0.0021641953  0.0003773637  0.0018222099
4 rs3748592  0.0011018723  0.0008842415  0.0012244474
5 rs3748593  0.0022335309 -0.0001147579  0.0014497337
6 rs2272756  0.0010033806  0.0011124711  0.0019509467

# PRS (with file_predGeno provided): a data frame storing the PRS generated using mu1, mu2 and mu_XPASS, respectively.
> head(fit_bbj$PRS)
      FID     IID        PRS1       PRS2   PRS_XPASS
1 HG00403 HG00403 -1.23431226  0.5992838 -0.29932176
2 HG00404 HG00404 -0.28765365 -0.8408158 -0.02486850
3 HG00406 HG00406  0.20704538 -1.8499626  0.03852178
4 HG00407 HG00407 -0.26639864 -1.2527320  0.00212844
5 HG00409 HG00409 -0.06152861  0.9371440  0.32425772
6 HG00410 HG00410 -0.14512120 -2.5842302 -0.75355484

```
XPASS will also write above outputs into the files with `file_out` prefix, if provided.


 # Usage
 
To fit XPASS using your own datasets, follow the steps below:

 Step 1: Download GWAS summary-level data of both target and auxiliary populations from public resources
 
 Step 2: Prepare summary statistics files in the XPASS format (check data format in the [example data](https://www.dropbox.com/sh/i7rhnko69974dje/AACfcDXz0cmwshbli8q7PZA5a?dl=0))
 
 Step 3: Prepare reference genotypes and assiciated covariates from target and auxiliary populations
 
 Step 4: Fit XPASS



# Development
The XPASS package is developed by Mingxuan Cai (mcaiad@ust.hk).

# Contact information

Please contact Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.
