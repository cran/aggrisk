`unweighted` <-
function(Z, mu, beta0){

  # Create the objective function:

  objfun <- function(beta){

    # Sum up the m values of exp(-Zbeta) multiplied
    # to obtain a scalar
    term1 <- sum(exp(-Z%*%beta))

    # Sum the columns of mu to get a 1*p vector, and multiply
    # by beta to get a scalar.
    term2 <- apply(mu,2,sum)%*%beta

    # Sum the two terms
    out <- term1+term2

    return(out)
  }

  # Need to propose initial values for the vector beta
  # Specify number of cases
  ncases <- dim(Z)[1]

  # Specify number of parameters
  p <- dim(Z)[2]


  # Use optim to minimise the objfun:
  beta_uw <- nlminb(beta0, objfun, control=list(iter.max=1500,eval.max=1500),
                  lower=-100,upper=100)


  return(beta_uw)
}

