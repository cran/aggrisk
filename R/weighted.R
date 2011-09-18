`weighted` <-
function(Z, mu_t, weight_all, weight_data,beta0){

    p <- dim(Z)[2]
    n <- dim(mu_t)[1]/p

  # Create the objective function:

  objfun_w <- function(beta){

    # Sum up the m values of exp(-Zbeta) multiplied by weights
    # to obtain a scalar
    term1 <- sum(exp(-Z%*%beta)*weight_data)

    term2 <- 0

    for(i in 1:p){
        X=mu_t[((i-1)*n+1):(i*n),]
        term2 <- term2 + sum(apply(X*weight_all,2,sum)*beta[i])
    }

    out <-term1+term2
    return(out)
  }


 

  # Use optim to minimise the objfun:
  beta_w <- nlminb(beta0, objfun_w, control=list(iter.max=1500,eval.max=1500),
                  lower=-100,upper=100)


  return(beta_w)
}

