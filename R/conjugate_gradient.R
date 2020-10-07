# A^(-1)b=x by conjugate gradient
conjugate_gradient <- function(A,b,x=NULL,maxIter=500,verbose=T) {
  if(is.null(x)) x <- rep(0,ncol(A))
  r <- b-A%*%x
  p <- r
  rsold <- sum(r^2)
  for(iter in 1:maxIter){
    Ap <- A%*%p
    alpha <- sum(r^2)/sum(p*Ap)
    x <- x + alpha*p
    r <- r - alpha*Ap
    rsnew <- sum(r^2)
    if(verbose) cat(iter,"-th iteration: residual:",sqrt(rsnew),"\n")
    if(sqrt(rsnew)<1e-6) break
    beta <- rsnew/rsold
    p <- r + beta*p
    rsold <- rsnew
  }
  return(x)
}
