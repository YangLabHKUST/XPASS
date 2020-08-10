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
