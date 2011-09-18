`testf` <-
function(controlid, caseid, counts, lambda, lags, xcoord, ycoord, h){

  # Match the small areas in the case data to the small areas in the control
  # data
  postid <- NULL
  for(i in 1:length(caseid)){
    postid[i]=which(controlid==caseid[i])
  }

  # Get a seqence of 1:no of small areas
  seqOA <- 1:length(controlid)


  # Calculate the denominator term
  den=rep(0,length(lags))
  for(i in 1:length(controlid)){
    centre_x=xcoord[i]
    centre_y=ycoord[i]
    distances=sqrt((centre_x-xcoord)^2+(centre_y-ycoord)^2)
    
    for(j in 1:length(lags)){
        lag=lags[j]
        ids=which(distances<=lag+h & distances>=lag-h)
        den[j]=den[j]+sum(counts[seqOA[i]]*counts[seqOA[ids]])
    }
  }

  # Calculate the numerator term
  num=rep(0,length(lags))

  xcoord=xcoord[postid]
  ycoord=ycoord[postid]
  for(i in 1:length(postid)){
    centre_x=xcoord[i]
    centre_y=ycoord[i]
    distances=sqrt((centre_x-xcoord)^2+(centre_y-ycoord)^2)
    
    for(j in 1:length(lags)){
        lag=lags[j]
        ids=which(distances<=lag+h & distances>=lag-h)
        G=sum(1/lambda[ids]/lambda[i])-1/lambda[i]^2
        num[j]=num[j]+G          
    }
  }

  pcf=num/den
  
  return(pcf)
}

