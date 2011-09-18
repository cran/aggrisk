teststat <-
function(controlid, caseid, lambda, counts){

  T1_w=rep(0,length(controlid))
  for(j in 1:length(controlid)){
      id_j=which(caseid==controlid[j])    
      R_j=lambda[id_j]
      if(length(id_j)>0){
        for(k in 1:length(id_j)){
            T1_w[j]=T1_w[j]+1/R_j[k]*(sum(1/R_j)-1/R_j[k])
        }
      }
  }

  return(sum(T1_w)-sum(counts^2))
}

