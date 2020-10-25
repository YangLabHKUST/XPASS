# construct PRS with new genotype file
predict_XPASS <- function(pm,file_predgeno){
  geno <- read_data(file_predgeno,"zero")
  X <- geno$X
  rm(geno)
  ref <- fread(paste0(file_predgeno,".bim"),data.table=F)
  fam <- fread(paste0(file_predgeno,".fam"),data.table=F)

  # find common SNPs in test and training data
  snps <- intersect(pm$SNP,ref$V2)
  if(length(snps)<nrow(ref)){
    warning("Test data and training data have different SNPs. This may reduce prediction accuracy. To avoid this problem, pecify the test genotype file in the `file_predgeno` argument in XPASS function to pre-mathc SNPs.")
  }

  # match SNPs in both datasets
  idx_snps <- match(snps,ref$V2)
  X <- X[,idx_snps]
  ref <- ref[idx_snps,]

  idx_snps <- match(snps,pm$SNP)
  pm <- pm[idx_snps,]

  # flip alleles
  idx_flip1 <- which(ref$V5%in%c("A","T")&pm$A1%in%c("C","G"))
  idx_flip2 <- which(ref$V5%in%c("C","G")&pm$A1%in%c("A","T"))
  idx_flip <- c(idx_flip1,idx_flip2)

  X[,idx_flip] <- 2-X[idx_flip]

  # compute PRS
  PRS <- X %*% data.matrix(pm[,c("mu1","mu2","mu_XPASS")])
  colnames(PRS) <- c("PRS1","PRS2","PRS_XPASS")
  PRS <- data.frame(FID=fam$V1,IID=fam$V2,PRS)

}


# evalueate R2 using external summary statistics
evalR2_XPASS <- function(pm,file_z_pred,file_predgeno){
  geno <- read_data(file_predgeno,"zero")
  X <- geno$X
  rm(geno)
  ref <- fread(paste0(file_predgeno,".bim"),data.table=F)

  zf <- fread(file_z_pred,data.table=F)

  # find common SNPs in pm, test sumstats and test geno
  snps <- intersect(pm$SNP,ref$V2)
  snps <- intersect(snps,zf$SNP)
  if(length(snps)<nrow(ref)){
    warning("Test data and training data have different SNPs. This may reduce prediction accuracy. To avoid this problem, pecify the test genotype file in the `file_predgeno` argument in XPASS function to pre-match SNPs.")
  }

  # match SNPs in pm, test sumstats and test geno
  idx_snps <- match(snps,ref$V2)
  X <- X[,idx_snps]
  ref <- ref[idx_snps,]

  idx_snps <- match(snps,pm$SNP)
  pm <- pm[idx_snps,]

  idx_snps <- match(snps,zf$SNP)
  zf <- zf[idx_snps,]

  # flip alleles in test geno
  idx_flip1 <- which(ref$V5%in%c("A","T")&pm$A1%in%c("C","G"))
  idx_flip2 <- which(ref$V5%in%c("C","G")&pm$A1%in%c("A","T"))
  idx_flip <- c(idx_flip1,idx_flip2)

  X[,idx_flip] <- 2-X[idx_flip]

  # flip alleles in test sumstats
  zs <- zf$Z
  idx_flip1 <- which(zf$V5%in%c("A","T")&pm$A1%in%c("C","G"))
  idx_flip2 <- which(zf$V5%in%c("C","G")&pm$A1%in%c("A","T"))
  idx_flip <- c(idx_flip1,idx_flip2)

  zs[idx_flip] <- -zs[idx_flip]

  # compute R2
  betas <- data.matrix(pm[,c("mu1","mu2","mu_XPASS")])
  denom <- sqrt(colMeans((X%*%betas)^2))

  xsd <- apply(X,2,sd)
  R2 <- (colSums(zs*betas*xsd/sqrt(median(zf$N)))/denom)^2
  names(R2) <- c("PRS1","PRS2","PRS_XPASS")

  return(R2)
}
