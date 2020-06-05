#  X1 and X2 are standardized with mean 0 and variance 1/(p*m)
compute_pm <- function(z1,z2=NULL,X1=NULL,X2=NULL,h1,h2=NULL,h12=NULL,n1,n2=NULL,p,LDM1=t(X1)%*%X1,LDM2=t(X2)%*%X2){

  # compute mu from first dataset
  S1 <- LDM1
  diag(S1) <- diag(S1) + (1-h1)/h1/n1
  mu1 <- sqrt(n1/p) * (chol2inv(chol(S1))%*%z1)/n1
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

    # compute mu from second dataset
    mu2 <- sqrt(n2/p) * (invS2%*%z2)/n2

    # compute mu from both datasets using XPASS
    invA1 <- S1 - (1-h1)*(1-h2)*h12^2/denom^2/n1/n2 * invS2
    A1 <- chol2inv(chol(invA1))/n1

    mu_XPASS <- sqrt(n1/p) * A1%*%z1 + sqrt(n2/p) * (1-h1)*h12/denom/n2 * A1%*%(invS2%*%z2)

    ret <- cbind(mu1,mu2,mu_XPASS)
  }
  return(ret)
}
