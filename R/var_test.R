var_test <-
function(Z, lambda, beta, weights, counts, mu, caseid, controlid){

  #first calculate the variance of the statistic when the true parameter
  #estimate is being used

  V1 <- 0
  t1 <- 0
  t2 <- 0
  n <- dim(mu)[1]
  p <- dim(mu)[2]

  for(i in 1:n){
      id <- which(caseid==controlid[i])
      r <- lambda[id]
    
      T1 <- 0
      leng=length(id)
      if(leng>0){
        for(j in 1:leng){
            T1 <- T1+sum(1/(r^2))/(r[j]^2)-1/(r[j]^4)    
        }
        
        T2 <- counts[i]^2*sum(1/(r^2))

        V1 <- V1+2*T1+4*T2    
        t1 <- t1+T1
        t2 <- t2+T2
      }
  }

  #now estimate the variance due to an estimated regression parameter
  A <- matrix(0, nrow=p, ncol=p)
  B <- matrix(0, nrow=p, ncol=p)
  
  for(i in 1:p){
    for(j in 1:p){
      B[i,j] <- sum(Z[,i]*Z[,j]*weights*exp(-Z%*%beta))
      A[i,j] <- sum(Z[,i]*Z[,j]*(weights^2)*exp(-2*Z%*%beta))
    }
  }

  vbeta=solve(B)%*%A%*%solve(B)


  T3 <- rep(0,p)
  for(i in 1:n){
      T3 <- T3+counts[i]*mu[i,]
  }
  V2 <- 4*t(T3)%*%vbeta%*%T3    
    
  #now estimate the covariance between the two terms in the Taylor series
  #expansion

  T4 <- rep(0,p)
  for(i in 1:n){
    id <- which(caseid==controlid[i])
    z <- matrix(Z[id,],nc=p)
    r <- lambda[id]
    w <- weights[id]
    
    leng=length(id)
    if(leng>0){
      if(leng==1){
          T4 <- T4+w*counts[i]*z/(r^2)
      }
      if(leng>1){
          for(k in 1:p){
            T4[k] <- T4[k]+counts[i]*sum(w*z[,k]/(r^2))
          }
      }
    }
  }
  
  V3 <- -8*t(T3)%*%solve(B)%*%t(T4)

  #here is the final estimate
  
  return(V=V1+V2+V3)
}

