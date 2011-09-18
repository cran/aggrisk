#####################################################################
#   Function to calculate the weights for the Poisson process case  #
#####################################################################


# Variables needed to calculate the weights are:
# beta_uw: 	the values of the coefficients calculated using equal weights
# mu_t:	A matrix of aggregated control data in the format n*p by t
#		where n is the number of small areas, p is the number of
#		parameters in the model, and t is the time element.
# areaID	A vector small areas being used (controls)
# caseID	A vector of the small areas corresponding to the cases
# caseT	A vector of the t in which cases were diagnosed
# startT	Value of start time t
# endT	Value of end time t
# pop		Size of population of each small area 



find_weights <- function(beta_uw, mu_t, areaID, caseID, caseT, startT,endT,pop){

  p <- length(beta_uw)
  n <- dim(mu_t)[1]/p
  Tlength <- endT-startT+1 

  weight_all=matrix(0,nrow=n,ncol=Tlength)

  for(i in 1:n){
    for(j in 1:Tlength){
        if(pop[i,j]>0){
            weight_all[i,j] <- exp(beta_uw%*%mu_t[seq(i,n*p,by=n),j]/pop[i,j])
        }
    }   
  }

  ncases <- length(caseID)
  
  weight_data=rep(0,ncases)
  for(i in 1:ncases){
    id_region=which(areaID==caseID[i])
    id_time=(caseT[i]-startT)+1
    weight_data[i]=weight_all[id_region,id_time]  
  }

  return(list(weight_all=weight_all, weight_data=weight_data))
}



