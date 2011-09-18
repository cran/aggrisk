################################################################################
#               Function to create spline terms using cubic splines           #
################################################################################


sterms <- function(caseT, N, startT, endT, knots,pop){

  K <- length(knots)
  n <- length(caseT)


  # Obtain the spline terms for each of the cases
  spline_case <- matrix(0, nrow=n,ncol=K)

  for(i in 1:K){
    for(j in 1:n){
      if(caseT[j]>=knots[i]){
        spline_case[j,i] <- (caseT[j]-knots[i])^3
      }
    }
  }

  # Obtain fixed terms for the cases
  fixed_case <- cbind((caseT-startT), (caseT-startT)^2, (caseT-startT)^3)
  

  time <- startT:endT
  TIME <- endT-startT+1

  # Obtain spline terms for the controls
  spline_control <- NULL

  for(i in 1:K){
    cmat <- matrix(0,nrow=N, ncol=TIME)
    for(j in 1:N){
      for(k in 1:TIME){
        if(time[k]>=knots[i]){
          cmat[j,k] <- pop[j,k]*(time[k]-knots[i])^3
        }
      }
    }
    spline_control <- rbind(spline_control, cmat)
  }

  # fixed terms for controls

  fixed_control <- NULL

  for(i in 1:3){
    fixed <- matrix(0,nrow=N, ncol=TIME)
    for(j in 1:TIME){
      fixed[,j] <- pop[,j]*((j-1)^i)
    }
    fixed_control <- rbind(fixed_control,fixed)
  }

  # Controls summed over time

  sumfixed <- apply(fixed_control,1,sum)
  fixed_otime <- matrix(sumfixed,nrow=N,ncol=3, byrow=FALSE)  

  sumyr <- apply(spline_control,1,sum)
  spline_otime <- matrix(sumyr,nrow=N,ncol=K, byrow=FALSE)
  
    
  return(list(cases=cbind(fixed_case,spline_case), 
              control=rbind(fixed_control,spline_control),
              sumcontrol=cbind(fixed_otime,spline_otime)))
}  
