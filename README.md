# XPASS
The XPASS package implement the XPASS approach for generating PRS in a target population by integrating multi-ethnic datasets.

# Installation 
```{r}
#install.packages("devtools")
devtools::install_github("YangLabHKUST/XPASS")
```

# Quick start

We illustrate the usage of XPASS using the GWAS summary statistics of height from UKB and BBJ. For convenience, we use the 1000 Genomes project genotypes as reference panels, which may not achieve optimal prediction accuracy due to the limited sample size. In practice, it is suggested to use larger datasets as reference panels (n>2000). The datasets involved in the following example can be downloaded from [here](https://www.dropbox.com/sh/i7rhnko69974dje/AACfcDXz0cmwshbli8q7PZA5a?dl=0).

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

# genotype file of test data (plink prefix)
height_test <- "1000G.EAS.QC.hm3.ind"

# sumstats of height
height_bbj <- "height_bbj_3M_format.txt"  # target
height_ukb <- "height_ukb_3M_format.txt"  # auxiliary

fit_bbj <-XPASS(file_z1 = height_bbj,file_z2 = height_ukb,file_ref1 = ref_EAS,file_ref2 = ref_EUR,
                file_cov1 = cov_EAS,file_cov2 = cov_EUR,file_predGeno = height_test,
                pop = "EAS",sd_method="LD_block",compPosMean = T,
                file_out = "height_bbj_ukb_ref_TGP")

```

XPASS returns a list of results:
```{r}
# H: a table of estimated heritabilities, co-heritability and genetic correlation (first row)
# and their corresponding standard erros (second row).
> fit_bbj$H
             h1         h2        h12        rho
[1,] 0.43474307 0.63559207 0.37876061 0.72054189
[2,] 0.02087136 0.03403142 0.02031872 0.01727539

# mu: a data frame storing the posterior means obtained by LDpred-inf using only the target dataset (mu1) and
# only the auxiliary dataset (mu2), and the posterior mean obtained by XPASS (mu_XPASS).
> head(fit_bbj$mu)
        SNP           mu1           mu2      mu_XPASS
1 rs4475691 -0.0008581511  0.0004678541 -0.0004164065
2 rs7537756 -0.0017717779  0.0010456091 -0.0004210284
3 rs7523549  0.0021641953  0.0003773637  0.0018222099
4 rs3748592  0.0011018723  0.0008842415  0.0012244474
5 rs3748593  0.0022335309 -0.0001147579  0.0014497337
6 rs2272756  0.0010033806  0.0011124711  0.0019509467

# PRS (if file_predGeno is provided): a data frame storing the PRS generated using mu1, mu2 and mu_XPASS, respectively.
> head(fit_bbj$PRS)
      FID     IID        PRS1       PRS2   PRS_XPASS
1 HG00403 HG00403 -1.23431226  0.5992838 -0.29932176
2 HG00404 HG00404 -0.28765365 -0.8408158 -0.02486850
3 HG00406 HG00406  0.20704538 -1.8499626  0.03852178
4 HG00407 HG00407 -0.26639864 -1.2527320  0.00212844
5 HG00409 HG00409 -0.06152861  0.9371440  0.32425772
6 HG00410 HG00410 -0.14512120 -2.5842302 -0.75355484

```
XPASS will also write above outputs in the files with `file_out` prefix, if provided.


 # Usage
 
To fit XPASS using your own datasets, follow the steps below:

 Step 1: Download GWAS summary-level data of both target and auxiliary populations from public resources
 
 Step 2: Prepare summary statistics files in the XPASS format (check data format in the [example data](https://www.dropbox.com/sh/i7rhnko69974dje/AACfcDXz0cmwshbli8q7PZA5a?dl=0))
 
 Step 3: Prepare reference genotypes and assiciated covariates from target and auxiliary populations
 
 Step 4: Fit XPASS



# Development
The XPASS package is developed by Mingxuan Cai (mcaiad@ust.hk).

# Contact information

Please feel free to contact Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
