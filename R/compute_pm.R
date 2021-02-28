#  X1 and X2 are standardized with mean 0 and variance 1/(p*m)

# compute_pm uses schur complement. This function is NOT USED any more sicne it is slower
compute_pm <- function(z1,z2=NULL,X1=NULL,X2=NULL,h1,h2=NULL,h12=NULL,n1,n2=NULL,p,LDM1=t(X1)%*%X1,LDM2=t(X2)%*%X2,use_CG=T){

  # compute mu from first dataset
  S1 <- LDM1
  diag(S1) <- diag(S1) + (1-h1)/h1/n1

  if(use_CG){
    S1invz1 <- conjugate_gradient(S1,z1,verbose=F)
  } else {
    S1invz1 <- chol2inv(chol(S1))%*%z1
  }
  mu1 <- sqrt(n1/p) * S1invz1/n1

  ret <- matrix(mu1)
  if(is.null(z2)|is.null(h2)|is.null(h12)|is.null(n2)){
    cat("Compute posterior mean from one dataset.\n")
  } else {
    cat("Compute posterior mean from two datasets.\n")

    denom <- h1*h2-h12^2

    S2 <- LDM2
    diag(S2) <- diag(S2) + (1-h2)*h1/denom/n2

    S1 <- LDM1
    diag(S1) <- diag(S1) + (1-h1)*h2/denom/n1

    invS2 <- chol2inv(chol(S2))
    if(use_CG){
      S2invz2 <- conjugate_gradient(S2,z2,verbose=F)
    } else {
      S2invz2 <- invS2%*%z2
    }

    # compute mu from second dataset
    mu2 <- sqrt(n2/p) * S2invz2/n2

    # compute mu from both datasets using XPASS
    invA1 <- S1 - (1-h1)*(1-h2)*h12^2/denom^2/n1/n2 * invS2

    if(use_CG){
      A1z1 <- conjugate_gradient(invA1,z1,verbose=F)
      A1S2invz2 <- conjugate_gradient(invA1,S2invz2,verbose=F)
    } else{
      A1 <- chol2inv(chol(invA1))
      A1z1 <- A1%*%z1
      A1S2invz2 <- A1%*%S2invz2
    }

    mu_XPASS <- sqrt(n1/p) * A1z1/n1 + sqrt(n2/p) * (1-h1)*h12/denom/n2 * A1S2invz2/n1

    ret <- cbind(mu1,mu2,mu_XPASS)
  }
  return(ret)
}


#  X1 and X2 are standardized with mean 0 and variance 1/(p*m)
# compute_pm2 directly solves 2p by 2p linear system. It is faster.
compute_pm2 <- function(z1,z2=NULL,X1=NULL,X2=NULL,h1,h2=NULL,h12=NULL,n1,n2=NULL,p,
                        LDM1=t(X1)%*%X1,LDM2=t(X2)%*%X2,LDMslB1=NULL,LDMslB2=NULL,
                        use_CG=T){
  p_block <- ncol(LDM1)

  # compute mu1 from the first dataset
  S1 <- n1*LDM1
  diag(S1) <- diag(S1) + (1-h1)/h1

  tmp1 <- z1*sqrt(n1/p)
  if(!is.null(LDMslB1)){
    tmp1 <- tmp1-n1*LDMslB1
  }

  if(use_CG){
    S1invz1 <- conjugate_gradient(S1,tmp1,verbose=F)
  } else {
    S1invz1 <- chol2inv(chol(S1))%*%tmp1
  }
  mu1 <- S1invz1

  ret <- matrix(mu1)
  if(is.null(z2)|is.null(h2)|is.null(h12)|is.null(n2)){
    cat("Compute posterior mean from one dataset.\n")
  } else {
    cat("Compute posterior mean from two datasets.\n")

    # compute mu2 from the second dataset
    S2 <- n2*LDM2
    diag(S2) <- diag(S2) + (1-h2)/h2

    tmp2 <- z2*sqrt(n2/p)
    if(!is.null(LDMslB2)){
      tmp2 <- tmp2-n2*LDMslB2
    }

    if(use_CG){
      S2invz2 <- conjugate_gradient(S2,tmp2,verbose=F)
    } else {
      S2invz2 <- chol2inv(chol(S2))%*%tmp2
    }
    mu2 <- S2invz2

    # compute mu_XPASS from both datasets
    Delta <- matrix(c(h1/(1-h1),h12/(1-h1),h12/(1-h2),h2/(1-h2)),2,2)
    invDelta <- solve(Delta)

    S1 <- n1*LDM1
    diag(S1) <- diag(S1) + invDelta[1,1]

    S2 <- n2*LDM2
    diag(S2) <- diag(S2) + invDelta[2,2]

    S <- matrix(0,2*p_block,2*p_block)
    S[1:p_block,1:p_block] <- S1
    S[(p_block+1):(2*p_block),(p_block+1):(2*p_block)] <- S2
    S[(p_block+1):(2*p_block),1:p_block] <- diag(invDelta[2,1],p_block)
    S[1:p_block,(p_block+1):(2*p_block)] <- diag(invDelta[1,2],p_block)

    Sinvz <- solve(S,c(tmp1,tmp2))

    mu_XPASS1 <- Sinvz[1:p_block]
    mu_XPASS2 <- Sinvz[(p_block+1):(2*p_block)]

    ret <- cbind(mu1,mu2,mu_XPASS1,mu_XPASS2)
    # ret <- list(out=cbind(mu1,mu2,mu_XPASS1,mu_XPASS2),Delta=Delta,invDelta=invDelta,S=S,tmp1,tmp2=tmp2)
  }
  return(ret)
}




#  X1 and X2 are standardized with mean 0 and variance 1/(p*m)
# compute_pm2 directly solves 2p by 2p linear system. It is faster.
compute_fe <- function(zs1,zl1=NULL,zs2=NULL,zl2=NULL,Xs1=NULL,Xs2=NULL,Xl1=NULL,Xl2=NULL,h1,h2=NULL,h12=NULL,n1,n2=NULL,p,
                       LDMs1=t(Xs1)%*%Xs1,LDMl1=t(Xl1)%*%Xl1,LDMsl1=t(Xs1)%*%Xl1,
                       LDMs2=t(Xs2)%*%Xs2,LDMl2=t(Xl2)%*%Xl2,LDMsl2=t(Xs2)%*%Xl2,use_CG=T){
  p_small <- ncol(LDMs1)
  p_large1 <- ncol(LDMl1)
  p_large2 <- ncol(LDMl2)

  if(length(zl1>0)){
    # compute mu1 from the first dataset
    S1 <- n1*LDMs1
    diag(S1) <- diag(S1) + (1-h1)/h1

    tmp1 <- sqrt(n1/p) * zs1

    if(use_CG){
      Sinvz1 <- conjugate_gradient(S1,tmp1,verbose=F)
    } else {
      invS1 <- chol2inv(chol(S1))
      Sinvz1 <- invS1%*%tmp1
    }

    b1 <- sqrt(n1/p) * zl1-t(n1*LDMsl1)%*%Sinvz1

    if(use_CG){
      SinvLDMsl1 <- apply(n1*LDMsl1,2,function(b) conjugate_gradient(A=S1,b=b,verbose=F))
    } else {
      SinvLDMsl1 <- invS1%*%(n1*LDMsl1)
    }

    A1 <- n1*LDMl1 - t(n1*LDMsl1) %*% SinvLDMsl1

    beta1 <- conjugate_gradient(A1,b1,verbose = F)
  } else {
    beta1 <- numeric(0)
  }

  if(length(zl2>0)){
    # compute mu2 from the second dataset
    S2 <- n2*LDMs2
    diag(S2) <- diag(S2) + (1-h2)/h2

    tmp2 <- sqrt(n2/p) * zs2

    if(use_CG){
      Sinvz2 <- conjugate_gradient(S2,tmp2,verbose=F)
    } else {
      invS2 <- chol2inv(chol(S2))
      Sinvz2 <- invS2%*%tmp2
    }

    b2 <- sqrt(n2/p) * zl2-t(n2*LDMsl2)%*%Sinvz2

    if(use_CG){
      SinvLDMsl2 <- apply(n2*LDMsl2,2,function(b) conjugate_gradient(A=S2,b=b,verbose=F))
    } else {
      SinvLDMsl2 <- invS2%*%(n2*LDMsl2)
    }

    A2 <- n2*LDMl2 - t(n2*LDMsl2) %*% SinvLDMsl2

    beta2 <- conjugate_gradient(A2,b2,verbose = F)
  } else {
    beta2 <- numeric(0)
  }

  # compute mu_XPASS from both datasets
  Delta <- matrix(c(h1/(1-h1),h12/(1-h1),h12/(1-h2),h2/(1-h2)),2,2)
  invDelta <- solve(Delta)

  S1 <- n1*LDMs1
  diag(S1) <- diag(S1) + invDelta[1,1]

  S2 <- n2*LDMs2
  diag(S2) <- diag(S2) + invDelta[2,2]

  S <- matrix(0,2*p_small,2*p_small)
  S[1:p_small,1:p_small] <- S1
  S[(p_small+1):(2*p_small),(p_small+1):(2*p_small)] <- S2
  S[(p_small+1):(2*p_small),1:p_small] <- diag(invDelta[2,1],p_small)
  S[1:p_small,(p_small+1):(2*p_small)] <- diag(invDelta[1,2],p_small)

  tmp1 <- sqrt(n1/p) * zs1
  tmp2 <- sqrt(n2/p) * zs2

  Sinvz <- solve(S,c(tmp1,tmp2))

  b1 <- sqrt(n1/p)*zl1 - t(n1*LDMsl1)%*%Sinvz[1:p_small]
  b2 <- sqrt(n2/p)*zl2 - t(n2*LDMsl2)%*%Sinvz[(p_small+1):(2*p_small)]

  LDMsl <- adiag(n1*LDMsl1,n2*LDMsl2)
  if(use_CG){
    SinvLDMsl <- apply(LDMsl,2,function(b) solve(a=S,b=b))
  } else {
    SinvLDMsl <- invS%*%LDMsl
  }

  A <- adiag(n1*LDMl1,n2*LDMl2) - t(LDMsl) %*% SinvLDMsl

  beta_XPASS <- conjugate_gradient(A,c(b1,b2),verbose = F)

  # beta_XPASS1 <- LDMl1%*%tmp1 - t(LDMsl1)%*%Sinvtmp[1:p_small]/n1
  # beta_XPASS2 <- LDMl2%*%tmp2 - t(LDMsl2)%*%Sinvtmp[(p_small+1):(2*p_small)]/n2

  if(length(zl1)>0){
    beta_XPASS1 <- beta_XPASS[1:p_large1]
  } else{
    beta_XPASS1 <- numeric(0)
  }

  if(length(zl2)>0){
    beta_XPASS2 <- beta_XPASS[(p_large1+1):(p_large1+p_large2)]
  } else{
    beta_XPASS2 <- numeric(0)
  }

  ret <- list(beta1=beta1,beta2=beta2,beta_XPASS1=beta_XPASS1,beta_XPASS2=beta_XPASS2)

  return(ret)
}

