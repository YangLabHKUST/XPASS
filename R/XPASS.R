#'@title  Main function for XPASS
#'@description XPASS: a unified framework for cross-population trait prediction with GWAS summary statistics.
#'This function provides both XPASS and XPASS+ implementation
#'
#' @param file_z1  Z-score file for population one. Five columns are required with column names SNP, N, Z, A1, A2. A1 is the effect allele.
#' @param file_z2  Z-score file for population two. Same format with file_z1. Note that file_z1 and file_z2 are exchangeable, both can be target or auxiliary dataset.
#' @param file_ref1 Prefix of file name/path of reference genotypes corresponding to file_z1 in bed format. Exclude the .bim/.bed/.fam. extensions in the path name.
#' @param file_ref2 Prefix of file name/path of reference genotypes corresponding to file_z2 in bed format.
#' @param file_cov1 Covariates file associated with the reference genotypes file_ref1. The covariates files should not include row names and column names. The rows should be exactly corresponding to the individuals in the .fam file of reference genotypes. Each column corresponds to one covariate. A column of one should not be included in the file.
#' @param file_cov2 Covariates file associated with the reference genotypes file_ref2. Same format with file_cov2.
#' @param file_predGeno File name/path of the genotypes of prediction dataset.
#' @param K1 n1 by n1 kingship matrix computed from file_ref1. Can improve computational efficiency if provided. If not provided, K1 will be computed from file_ref1.
#' @param K2 n2 by n2 kingship matrix computed from file_ref2. Can improve computational efficiency if provided. If not provided, K2 will be computed from file_ref2.
#' @param K12 n1 by n2 Cross-population kingship matrix computed from file_ref1 and file_ref2. Can improve computational efficiency if provided. If not provided, K12 will be computed from file_ref1 and file_ref2.
#' @param X1 n1 by p Raw genotype matrix from file_ref1. Can improve computational efficiency if provided. If not provided, X1 will be read from file_ref1.
#' @param X2 n2 by p Raw genotype matrix from file_ref2. Can improve computational efficiency if provided. If not provided, X2 will be read from file_ref2.
#' @param snps_fe1 A vector containing the SNP rs id to be included as population one-specific fixed effects in XPASS+.
#' @param snps_fe2 A vector containing the SNP rs id to be included as population two-specific fixed effects in XPASS+.
#' @param snps_list A vector containing the SNP rs id to be included in the analysis. The intersection of this list and SNPs from file_z1, file_z2, file_ref1 and file_ref2 will be used in the XPASS analysis.
#' @param sd_method A character string indicating the method to obtain the Jackknife standard error of parameter estimates. Should be one of "Chromosome', "LD_block" and "Jackknife", corresponding to chromosome based, LD-block (of population one) based, and SNP based Jackknife.
#' @param pop A character string indicating from which population the LD partions should be used. Currently support "EUR" and "EAS". Default is "EUR".
#' @param compPosMean TRUE/FALSE indicating whether to compute the posterior mean of effect sizes. Default is TRUE.
#' @param use_CG TRUE/FALSE indicating whether to use conjugate gradient to solve for the posterior mean. Default is TRUE.
#' @param compPRS TRUE/FALSE indicating whether to compute the PRS when file_predGeno is not NULL. If FALSE, one can use the function predict_XPASS to obtain the PRS later.
#' @param file_out A character string of the prefix of output file name.
#'
#'
#' @return a list with the following elements:
#' \describe{
#' \item{H: }{4 by 2 matrix of estimated heritabilities, co-heritability and genetic correlation (first row) and their corresponding standard erros (second row)}
#' \item{mu: }{A data frame storing the posterior means computed by LDpred-inf using only the target dataset (mu1) and only the auxiliary dataset (mu2), and the posterior means of the target population and the auxiliary population computed by XPASS (mu_XPASS1 and mu_XPASS2). In eff_type1 and eff_type2 columns, "RE" indicates the SNP effect is random effect and "FE" indicates the effect is treated as population-specific fixed effect in the corresponding populations. SNPs information is also returned: A1 is the effect allele, A2 is the other allele.}
#' \item{PRS (if file_predGeno provided and compPRS=T): }{A data frame storing the PRS generated using mu1, mu2, mu_XPASS1, and mu_XPASS2, respectively.}
#' \item{snps: }{Final SNPs included in the analysis.}
#' \item{fit: }{An object returned by the parameter estimation function corr_ss.}
#' }
#'
#' @examples
#' See vignette.
#' @export

XPASS <- function(file_z1,file_z2,file_ref1,file_ref2=NULL,file_cov1=NULL,file_cov2=NULL,file_predGeno=NULL,K1=NULL,K2=NULL,K12=NULL,X1=NULL,X2=NULL,
                  snps_fe1=NULL,snps_fe2=NULL,snp_list=NULL,
                  sd_method="Chromosome",pop="EUR",compPosMean=T,use_CG=T,compPRS=F,file_out=""){
  if(nchar(file_out)>0){
    cat("Writing to log file: ",file_out,".log\n",sep="")
    sink(paste0(file_out,".log"),append=F,split=T)
  }
  cat("Summary statistics file 1: ",file_z1,"\n",sep="")
  cat("Summary statistics file 2: ",file_z2,"\n",sep="")
  cat("Reference file 1: ",file_ref1,"\n",sep="")
  cat("Reference file 2: ",file_ref2,"\n",sep="")
  if(! is.null(file_cov1)) cat("Covariates file 1: ",file_cov1,"\n",sep="")
  if(! is.null(file_cov2)) cat("Covariates file 2: ",file_cov2,"\n",sep="")
  if(! is.null(file_predGeno)) cat("Test genotype file : ",file_predGeno,"\n",sep="")

  cat("Reading data from summary statisitcs...\n")
  zf1 <- fread(file_z1,data.table = F,header = T,stringsAsFactors = F)
  zf2 <- fread(file_z2,data.table = F,header = T,stringsAsFactors = F)
  # colnames(zf1) <- colnames(zf2) <- c("SNP","N","Z","A1","A2")
  cat(nrow(zf1)," and ",nrow(zf2)," SNPs found in summary statistics files 1 and 2.\n",sep = "")

  cat("Reading SNP info from reference panels...\n")
  ref1_info <- fread(paste(file_ref1,".bim",sep=""),data.table = F,stringsAsFactors = F)
  snps_info <- ref1_info
  ref2_info <- fread(paste(file_ref2,".bim",sep=""),data.table = F,stringsAsFactors = F)
  colnames(ref1_info) <- colnames(ref2_info) <- c("CHR","SNP","POS","BP","A1","A2")
  cat(nrow(ref1_info)," and ",nrow(ref2_info)," SNPs found in reference panel 1 and 2.\n",sep = "")

  # matching SNPs
  snps <- intersect(zf1$SNP,ref1_info$SNP)
  snps <- intersect(snps,zf2$SNP)
  snps <- intersect(snps,ref2_info$SNP)
  if(!is.null(snp_list)) snps <- snps[snps%in%snp_list]    # extract SNPs in the provided list if any
  if(!is.null(file_predGeno)&compPosMean){                 # match SNPs in test genotype if provided
    cat("Reading SNP info from test genotype file...\n")
    test_info <- fread(paste(file_predGeno,".bim",sep=""),data.table = F,stringsAsFactors = F)
    colnames(test_info) <-  c("CHR","SNP","POS","BP","A1","A2")
    cat(nrow(test_info)," SNPs found in test file.\n",sep = "")

    snps <- intersect(test_info$SNP,snps)
  }

  zf1 <- zf1[match(snps,zf1$SNP),]
  zf2 <- zf2[match(snps,zf2$SNP),]
  idx_ref1 <- match(snps,ref1_info$SNP)
  ref1_info <- ref1_info[idx_ref1,]
  snps_info <- snps_info[idx_ref1,]
  idx_ref2 <- match(snps,ref2_info$SNP)
  ref2_info <- ref2_info[idx_ref2,]
  if(!is.null(file_predGeno)&compPosMean){
    idx_test <- match(snps,test_info$SNP)
    test_info <- test_info[idx_test,]
  }
  cat(length(snps)," SNPs are matched in all files.\n",sep = "")

  # replace T with A, replace G with C; A=1, C=2
  zf1$A1[zf1$A1=="T"] <- "A"
  zf1$A2[zf1$A2=="T"] <- "A"
  zf1$A1[zf1$A1=="G"] <- "C"
  zf1$A2[zf1$A2=="G"] <- "C"
  zf1_A1 <- ifelse(zf1$A1=="A",1,2)
  zf1_A2 <- ifelse(zf1$A2=="A",1,2)

  ref1_info$A1[ref1_info$A1=="T"] <- "A"
  ref1_info$A2[ref1_info$A2=="T"] <- "A"
  ref1_info$A1[ref1_info$A1=="G"] <- "C"
  ref1_info$A2[ref1_info$A2=="G"] <- "C"
  ref1_A1 <- ifelse(ref1_info$A1=="A",1,2)
  ref1_A2 <- ifelse(ref1_info$A2=="A",1,2)


  zf2$A1[zf2$A1=="T"] <- "A"
  zf2$A2[zf2$A2=="T"] <- "A"
  zf2$A1[zf2$A1=="G"] <- "C"
  zf2$A2[zf2$A2=="G"] <- "C"
  zf2_A1 <- ifelse(zf2$A1=="A",1,2)
  zf2_A2 <- ifelse(zf2$A2=="A",1,2)

  ref2_info$A1[ref2_info$A1=="T"] <- "A"
  ref2_info$A2[ref2_info$A2=="T"] <- "A"
  ref2_info$A1[ref2_info$A1=="G"] <- "C"
  ref2_info$A2[ref2_info$S2=="G"] <- "C"
  ref2_A1 <- ifelse(ref2_info$A1=="A",1,2)
  ref2_A2 <- ifelse(ref2_info$A2=="A",1,2)

  # remove ambiguous SNPs
  snps_rm <- (ref1_A1+ref1_A2)!=(zf1_A1+zf1_A2) | (ref1_A1+ref1_A2)!=(zf2_A1+zf2_A2) | (ref1_A1+ref1_A2)!=(ref2_A1+ref2_A2) | (ref2_A1+ref2_A2)!=(zf1_A1+zf1_A2) | (ref2_A1+ref2_A2)!=(zf2_A1+zf2_A2) | (zf2_A1+zf2_A2)!=(zf1_A1+zf1_A2)

  if(!is.null(file_predGeno)&compPosMean){
    # replace T with A, replace G with C; A=1, C=2
    test_info$A1[test_info$A1=="T"] <- "A"
    test_info$A2[test_info$A2=="T"] <- "A"
    test_info$A1[test_info$A1=="G"] <- "C"
    test_info$A2[test_info$A2=="G"] <- "C"
    test_A1 <- ifelse(test_info$A1=="A",1,2)
    test_A2 <- ifelse(test_info$A2=="A",1,2)
    snps_rm <- snps_rm | (test_A1+test_A2)!=(zf1_A1+zf1_A2)
  }

  snps <- snps[!snps_rm]
  idx_ref1 <- idx_ref1[!snps_rm]
  idx_ref2 <- idx_ref2[!snps_rm]
  ref1_info <- ref1_info[!snps_rm,]
  ref2_info <- ref2_info[!snps_rm,]
  snps_info <- snps_info[!snps_rm,]

  if(!is.null(file_predGeno)&compPosMean){
    idx_test <- idx_test[!snps_rm]
    test_info <- test_info[!snps_rm,]
  }

  zf1 <- zf1[!snps_rm,]
  zf2 <- zf2[!snps_rm,]

  cat(sum(snps_rm)," SNPs are removed because of ambiguity; ",sum(!snps_rm)," SNPs remained.\n",sep = "")

  # if compute posterior mean, need raw genotypes
  if(compPosMean){
    if(is.null(X1)){
      ref1 <- read_data(file_ref1,fillMiss = "zero")
      X1 <- ref1$X[,idx_ref1]
    } else {
      X1 <- X1[,idx_ref1]
    }

    if(is.null(X2)){
      ref2 <- read_data(file_ref2,fillMiss = "zero")
      X2 <- ref2$X[,idx_ref2]
    } else {
      X2 <- X2[,idx_ref2]
    }
  }


  # calculate kinship matrix if needed
  if(is.null(K1)&is.null(K2)){
    cat("Calculating kinship matrix from the both reference panels...\n")
    if(is.null(X1)){
      ref1 <- read_data(file_ref1,fillMiss = "zero")
      X1 <- ref1$X[,idx_ref1]
    }

    if(is.null(X2)){
      ref2 <- read_data(file_ref2,fillMiss = "zero")
      X2 <- ref2$X[,idx_ref2]
    }

    no_var <- apply(X1,2,function(x) all(diff(x)==0)) | apply(X2,2,function(x) all(diff(x)==0))
    if(sum(no_var)!=0){
      snps <- snps[!no_var]
      X1 <- X1[,!no_var]
      X2 <- X2[,!no_var]
      ref1_info <- ref1_info[!no_var,]
      ref2_info <- ref2_info[!no_var,]
      snps_info <- snps_info[!no_var,]
      zf1 <- zf1[!no_var,]
      zf2 <- zf2[!no_var,]
      cat(sum(no_var)," SNPs are removed because of no variation; ",length(snps)," SNPs remained.\n",sep = "")
    }

    # align alleles in two reference panels based on the first ref
    ind_ref <- ref1_info$A1!=ref2_info$A1
    if(sum(ind_ref)>0){
      X2[,ind_ref] <- 1-(X2[,ind_ref]-1)
      cat(sum(ind_ref)," SNPs in the second reference panel are alligned for alleles according to the first.\n",sep="")
    }

    # X1sd <- apply(X1,2,sd)
    # X1 <- scale(X1)/sqrt(ncol(X1))
    # K1 <- X1%*%t(X1)
    # X2sd <- apply(X2,2,sd)
    # X2 <- scale(X2)/sqrt(ncol(X2))
    # K2 <- X2%*%t(X2)

    tmp1 <- scaleC(X1)
    X1sd <- c(tmp1$Xs)
    X1 <- tmp1$X
    K1 <- X1%*%t(X1)

    tmp2 <- scaleC(X2)
    X2sd <- c(tmp2$Xs)
    X2 <- tmp2$X
    K2 <- X2%*%t(X2)

    K12 <- X1%*%t(X2)

  } else if (is.null(K1)){
    cat("Calculating kinship matrix from the first reference panel...\n")
    if(is.null(X1)){
      ref1 <- read_data(file_ref1,fillMiss = "zero")
      X1 <- ref1$X[,idx_ref1]
    }

    no_var <- apply(X1,2,function(x) all(diff(x)==0))
    if(sum(no_var)!=0){
      snps <- snps[!no_var]
      X1 <- X1[,!no_var]
      ref1_info <- ref1_info[!no_var,]
      ref2_info <- ref2_info[!no_var,]
      snps_info <- snps_info[!no_var,]
      zf1 <- zf1[!no_var,]
      zf2 <- zf2[!no_var,]
      cat(sum(no_var)," SNPs are removed because of no variation; ",ncol(X1)," SNPs remained.\n",sep = "")
    }
    # X1sd <- apply(X1,2,sd)
    # X1 <- scale(X1)/sqrt(ncol(X1))
    # K1 <- X1%*%t(X1)

    tmp1 <- scaleC(X1)
    X1sd <- c(tmp1$Xs)
    X1 <- tmp1$X
    K1 <- X1%*%t(X1)
  } else if (is.null(K2)){
    cat("Calculating kinship matrix from the second reference panel...\n")
    if(is.null(X2)){
      ref2 <- read_data(file_ref2,fillMiss = "zero")
      X2 <- ref2$X[,idx_ref2]
    }

    no_var <- apply(X2,2,function(x) all(diff(x)==0))
    if(sum(no_var)!=0){
      snps <- snps[!no_var]
      X2 <- X2[,!no_var]
      ref1_info <- ref1_info[!no_var,]
      ref2_info <- ref2_info[!no_var,]
      snps_info <- snps_info[!no_var,]
      zf1 <- zf1[!no_var,]
      zf2 <- zf2[!no_var,]
      cat(sum(no_var)," SNPs are removed because of no variation; ",ncol(X2)," SNPs remained.\n",sep = "")
    }
    # X2sd <- apply(X2,2,sd)
    # X2 <- scale(X2)/sqrt(ncol(X2))
    # K2 <- X2%*%t(X2)

    tmp2 <- scaleC(X2)
    X2sd <- c(tmp2$Xs)
    X2 <- tmp2$X
    K2 <- X2%*%t(X2)
  } else {
    tmp1 <- scaleC(X1)
    X1sd <- c(tmp1$Xs)
    X1 <- tmp1$X

    tmp2 <- scaleC(X2)
    X2sd <- c(tmp2$Xs)
    X2 <- tmp2$X

    no_var <- rep(F,nrow(zf1))
  }


  # based on the minor allel of reference panel, correct for the z score in summary statistics
  ind1 <- zf1$A1 != ref1_info$A1
  ind2 <- zf2$A1 != ref1_info$A1
  z_score1 <- zf1$Z
  z_score2 <- zf2$Z
  z_score1[ind1] <- -z_score1[ind1]
  z_score2[ind2] <- -z_score2[ind2]
  cat(sum(ind1)," SNPs have different minor alleles in population 1, z-scores are corrected according to reference panel.\n",sep = "")
  cat(sum(ind2)," SNPs have different minor alleles in population 2, z-scores are corrected according to reference panel.\n",sep = "")

  if(sd_method=="LD_block"|compPosMean){
    cat("Assigning SNPs to LD Blocks...\n")
    block <- read.table(system.file("extdata", paste0(pop,"_fourier_ls-all.bed"), package = "XPASS"),header = T)
    group <- rep(0,nrow(zf1))
    idx_group <- 1
    for(i in 1:22){
      block_i <- block[block$chr==paste("chr",i,sep=""),]
      # chr_i <- snp_ref[snp_ref$CHR==i,]
      n_block <- nrow(block_i)

      for(j in 1:n_block){
        tmp <- with(ref1_info,CHR==i & BP>=block_i$start[j] & BP<block_i$stop[j])
        if(sum(tmp!=0)){
          group[tmp] <- idx_group
          idx_group <- idx_group+1
        }
      }
    }
  } else if(sd_method=="Jackknife"){
    group <- NULL
  } else if(sd_method=="Chromosome"){
    group <- ref1_info$CHR
  }

  if(is.null(file_cov1)){
    cov1 <- NULL
  } else {
    cov1 <- fread(file_cov1,data.table = F)
    cov1 <- data.matrix(cov1)
  }
  if(is.null(file_cov2)){
    cov2 <- NULL
  } else {
    cov2 <- fread(file_cov2,data.table = F)
    cov2 <- data.matrix(cov2)
  }

  cat("Calculate PVE...\n")
  isfe1 <- snps_info$V2%in%snps_fe1
  isfe2 <- snps_info$V2%in%snps_fe2
  # isre <- !(isfe1|isfe2)
  # fit <- corr_ss(z_score1[isre],z_score2[isre],K1,K2,K12,zf1$N[isre],zf2$N[isre],Z1=cbind(cov1,X1[,isfe1]),Z2=cbind(cov2,X2[,isfe2]),group = group)
  # fit <- corr_ss(z_score1,z_score2,K1,K2,K12,zf1$N,zf2$N,Z1=cbind(cov1,X1[,isfe1]),Z2=cbind(cov2,X2[,isfe2]),group = group)
  fit <- corr_ss(z_score1,z_score2,K1,K2,K12,zf1$N,zf2$N,Z1=cov1,Z2=cov2,group = group)

  print(fit$H)
  ret <- c(fit,snps=list(snps))

  if(compPosMean){
    # Set rho=0.99 (-0.99) if it exceed 1 (-1).
    # h12 <- ifelse(fit$H[1,4]>1,fit$H[1,1]*fit$H[1,2],fit$H[1,3])
    if(abs(fit$H[1,4])>1){
      cat("The MoM estimate of genetic correlation rho is ",fit$H[1,4],", which exceeds boundary. Set rho=0.99 (-0.99) for computing posterior means. \n")
      h12 <- sign(fit$H[1,4])*0.99*sqrt(fit$H[1,1]*fit$H[1,2])
    } else {
      h12 <- fit$H[1,3]
    }
    ngroup <- length(unique(group))
    mu <- matrix(0,ncol(X1),4)
    eff_type <- matrix("RE",ncol(X1),2)
    cat("Compute posterior mean from ",ngroup," blocks ...\n")
    start_j <- 1
    for(j in 1:ngroup){
      end_j <- start_j+sum(group==j)-1
      idx_j <- start_j:end_j

      isfe1 <- snps_info$V2[idx_j]%in%snps_fe1
      isfe2 <- snps_info$V2[idx_j]%in%snps_fe2

      LDMs1=t(X1[,idx_j])%*%X1[,idx_j]/nrow(X1)
      LDMs2=t(X2[,idx_j])%*%X2[,idx_j]/nrow(X2)
      if(sum(isfe1>0)|sum(isfe2)>0){
        # if some snps effects are fixed effects

        idx_fe1 <- idx_j[isfe1]
        idx_fe2 <- idx_j[isfe2]

        LDMl1=t(X1[,idx_fe1])%*%X1[,idx_fe1]/nrow(X1)
        LDMsl1=t(X1[,idx_j])%*%X1[,idx_fe1]/nrow(X1)

        LDMl2=t(X2[,idx_fe2])%*%X2[,idx_fe2]/nrow(X2)
        LDMsl2=t(X2[,idx_j])%*%X2[,idx_fe2]/nrow(X2)

        out_fe <- compute_fe(zs1 = z_score1[idx_j],zs2 = z_score2[idx_j],
                             zl1 = z_score1[idx_fe1],zl2 = z_score2[idx_fe2],
                             LDMs1=LDMs1,LDMl1=LDMl1,LDMsl1=LDMsl1,
                             LDMs2=LDMs2,LDMl2=LDMl2,LDMsl2=LDMsl2,
                             h1 = fit$H[1,1],h2 = fit$H[1,2],h12 = h12,n1 = median(zf1$N),n2 = median(zf2$N),p=ncol(X1),use_CG=use_CG)

        LDMslB1 <- LDMsl1 %*% out_fe$beta_XPASS1
        LDMslB2 <- LDMsl2 %*% out_fe$beta_XPASS2

        mu[idx_fe1,1] <- out_fe$beta1
        mu[idx_fe2,2] <- out_fe$beta2
        mu[idx_fe1,3] <- out_fe$beta_XPASS1
        mu[idx_fe2,4] <- out_fe$beta_XPASS2

        eff_type[idx_fe1,1] <- "FE"
        eff_type[idx_fe2,2] <- "FE"
      } else{
        # if no fixed effects, all snps are random effects
        LDMslB1 <- NULL
        LDMslB2 <- NULL
      }
      mu[idx_j,] <- mu[idx_j,] + compute_pm2(z1 = z_score1[idx_j],z2 = z_score2[idx_j],
                                             LDM1 = LDMs1,LDM2 = LDMs2,
                                             LDMslB1=LDMslB1,LDMslB2=LDMslB2,
                                             h1 = fit$H[1,1],h2 = fit$H[1,2],h12 = h12,n1 = median(zf1$N),n2 = median(zf2$N),p=ncol(X1),use_CG=use_CG)

      start_j <- end_j + 1
      cat("Block",j,"/",ngroup," finished. \n")
    }
    mu[,c(1,3)] <- mu[,c(1,3)]/sqrt(ncol(X1))/X1sd
    mu[,c(2,4)] <- mu[,c(2,4)]/sqrt(ncol(X2))/X2sd

    if(!is.null(file_predGeno)&compPRS){
      cat("Predicting PRS from test genotypes...\n")
      fam_test <- fread(paste0(file_predGeno,".fam"))
      test_geno <- read_data(file_predGeno,fillMiss = "zero")
      Xtest <- test_geno$X[,idx_test]

      # remove SNPs shhowing no variation in the training genotypes
      Xtest <- Xtest[,!no_var]

      # align alleles in test genotype based on the first ref
      ind_ref <- ref1_info$A1!=test_info$A1
      if(sum(ind_ref)>0){
        Xtest[,ind_ref] <- 1-(Xtest[,ind_ref]-1)
        cat(sum(ind_ref)," SNPs in the test genotypes are alligned for alleles according to the first.\n",sep="")
      }
      yh <- Xtest %*% mu
      yh <- data.frame(FID=fam_test$V1,IID=fam_test$V2,PRS1=yh[,1],PRS2=yh[,2],PRS_XPASS=yh[,3])
      write.table(yh,file=paste0(file_out,"_PRS.txt"),col.names = T,row.names = F,quote=F,sep="\t")
      ret <- c(ret,PRS=list(yh))
    }

    mu <- data.frame(CHR=snps_info$V1,SNP=snps_info$V2,POS=snps_info$V4,A1=snps_info$V5,A2=snps_info$V6,
                     mu1=mu[,1],mu2=mu[,2],mu_XPASS1=mu[,3],mu_XPASS2=mu[,4],eff_type1=eff_type[,1],eff_type2=eff_type[,2])
    write.table(mu,file=paste0(file_out,"_PosteriorMean.txt"),col.names = T,row.names = F,quote=F,sep="\t")
    ret <- c(ret,mu=list(mu))
  }

  cat("Done.\n")
  print(fit$H)

  if(nchar(file_out)>0) sink()
  return(ret)
}

