################################################################################
#          Function to find the variance of beta_hat                           #
################################################################################

var_beta <- function(beta,Z, weights){

  p <- length(beta)
  A <- matrix(0, nrow=p, ncol=p)
  B <- matrix(0, nrow=p, ncol=p)
  
  for(i in 1:p){
    for(j in 1:p){
      B[i,j] <- sum(Z[,i]*Z[,j]*weights*exp(-Z%*%beta))
      A[i,j] <- sum(Z[,i]*Z[,j]*(weights^2)*exp(-2*Z%*%beta))
    }
  }

  vbeta=solve(B)%*%A%*%solve(B)
   
  return(vbeta)
} 
