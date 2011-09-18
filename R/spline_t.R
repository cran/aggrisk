################################################################################
# Function for calculating the fixed terms and spline terms of the temporal    #
# trend once beta has been calculated                                          #
################################################################################

# beta:		Values of coefficients
# startT:		Value of t at start of the study
# endT:		Value of t at end of the study
# knots:		Value of t at each of the n knots
# ref1:		Position of the coefficients of the fixed terms (length 3)
# ref2:		Position of the coefficients of the spline terms (length n)
		 	

spline_t <- function(beta,startT,endT,knots,ref1, ref2){

  Tseq <- startT:endT
  Tref <- 0:(endT-startT)
  
  # Calculate fixed terms: t*beta[ref1[1]] + t^2*beta[ref1[2]] + t*beta[ref1[3]]

  fixedterms <- beta[ref1[1]]*Tref + beta[ref1[2]]*Tref^2 + beta[ref1[3]]*Tref^3

  splineterm <- rep(0, length(Tseq))
  term <- NULL

  for(i in 1:length(knots)){
    for(j in 1:length(Tseq)){
      if(Tseq[j]>=knots[i]){
        term[j] <- beta[ref2[i]]*(Tseq[j]-knots[i])^3}
      else{
        term[j] <- 0}
      }
    splineterm <- term+splineterm 
  }

  return(list(fixed=fixedterms, spline=splineterm))
}
