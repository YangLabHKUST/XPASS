corr_ss <- function(z1,z2,K1,K2,K12,n1,n2,Z1=NULL,Z2=NULL,group=NULL){
  m1 <- nrow(K1)
  m2 <- nrow(K2)
  p <- length(z1)
  ngroup <- length(unique(group))

  # calculate h1^2 for y1
  if(is.null(Z1)){
    Z1 <- matrix(1,m1,1)
    M1 <- diag(m1) - matrix(1/m1,m1,m1)
    MK1 <- K1
    MK12 <- K12
  } else{
    Z1 <- cbind(1,Z1)
    M1 <- diag(m1) - Z1%*%solve(t(Z1)%*%Z1)%*%t(Z1)
    MK1 <- M1%*%K1
    MK12 <- M1%*%K12
  }
  q1 <- ncol(Z1)

  trK1 <- sum(diag(MK1))
  trSQK1 <- trK1^2
  trK1SQ <- sum(MK1^2)

  S1 <- (trK1SQ-trSQK1/(m1-q1)) / (m1-q1)^2
  zz1 <- z1^2/n1 - 1/n1
  c1 <- sum(zz1) / p

  h1 <- c1/S1

  # calculate h2^2 for y2
  if(is.null(Z2)){
    Z2 <- matrix(1,m2,1)
    M2 <- diag(m2) - matrix(1/m2,m2,m2)
    MK2 <- K2
    K12M <- K12
  } else{
    Z2 <- cbind(1,Z2)
    M2 <- diag(m2) - Z2%*%solve(t(Z2)%*%Z2)%*%t(Z2)
    MK2 <- M2%*%K2
    K12M <- K12%*%M2
  }
  q2 <- ncol(Z2)

  trK2 <- sum(diag(MK2))
  trSQK2 <- trK2^2
  trK2SQ <- sum(MK2^2)

  S2 <- (trK2SQ-trSQK2/(m2-q2)) / (m2-q2)^2
  zz2 <- z2^2/n2 - 1/n2
  c2 <- sum(zz2) / p

  h2 <- c2/S2

  # calculate co-heritability
  trK12SQ <- sum(MK12 *K12M)# sum(diag(MK12%*%t(K12M)))

  if(is.logical(all.equal(K1,K2))){
    S3 <- S1
    # S3 <- trK12SQ/(m1-q1)/(m2-q2) - 1/sqrt(m1-q1)/sqrt(m2-q2)# + 1/sqrt(n1)/sqrt(n2)
  } else{
    S3 <- trK12SQ/(m1-q1)/(m2-q2)# - 1/sqrt(m1-q1)/sqrt(m2-q2) + 1/sqrt(n1)/sqrt(n2)
  }
  zz12 <- z1*z2/sqrt(n1)/sqrt(n2)
  c3 <- sum(zz12)/p

  h12 <- c3/S3

  rho <- h12/sqrt(h1*h2)

  if(is.null(group)){
    # compute standard error by Jackknife
    # zj <- (zs/sqrt(n))^2
    # c_jf <- (zz-zj)/(p-1) - 1/(median(n))
    c1_jf <- (sum(zz1)-zz1)/(p-1)
    var_h1 <- var(c1_jf)/S1/S1*(p-1)

    c2_jf <- (sum(zz2)-zz2)/(p-1)
    var_h2 <- var(c2_jf)/S2/S2*(p-1)

    c3_jf <- (sum(zz12)-zz12)/(p-1)
    var_h12 <- var(c3_jf)/S3/S3*(p-1)

    var_rho <- var(c3_jf/sqrt(c1_jf*c2_jf)) * S1*S2/S3/S3 * (p-1)
  } else{
    # zj <- sapply(1:ngroup,function(j){
    #
    #   tmp <- sum((zs[group==j]/sqrt(n[group==j]))^2)
    #   c(tmp,sum(group==j))
    # })
    # c_jf <- (zz-zj[1,])/(p-zj[2,]) - 1/(median(n))
    #
    # var_h <- var(c_jf)/S/S*(ngroup-1)
    zj <- sapply(1:ngroup,function(j){

      tmp1 <- sum(zz1[group==j])
      tmp2 <- sum(zz2[group==j])
      tmp3 <- sum(zz12[group==j])
      c(tmp1,tmp2,tmp3,sum(group==j))
    })
    c1_jf <- (sum(zz1)-zj[1,])/(p-zj[4,])
    var_h1 <- var(c1_jf)/S1/S1*(ngroup-1)

    c2_jf <- (sum(zz2)-zj[2,])/(p-zj[4,])
    var_h2 <- var(c2_jf)/S2/S2*(ngroup-1)

    c3_jf <- (sum(zz12)-zj[3,])/(p-zj[4,])
    var_h12 <- var(c3_jf)/S3/S3*(ngroup-1)

    var_rho <- var(c3_jf/sqrt(c1_jf*c2_jf)) * S1*S2/S3/S3 * (ngroup-1)
  }

  H <- matrix(0,2,4)
  H[1,] <- c(h1,h2,h12,rho)
  H[2,] <- sqrt(c(var_h1,var_h2,var_h12,var_rho))

  colnames(H) <- c("h1","h2","h12","rho")

  ret <- list(H=H,K1=K1,K2=K2,K12=K12,S1=S1,S2=S2,S3=S3,c1=c2,c2=c2,c3=c3,trK1SQ=trK1SQ,trK2SQ=trK2SQ,trK12SQ=trK12SQ,trSQK1=trSQK1,trSQK2=trSQK2,m1=m1,q1=q1,m2=m2,q2=q2)
}
