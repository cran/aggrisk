\name{weighted}
\alias{weighted}
\title{ This function uses estimating equations with weights selected via data-driven methods to estimate regression parameters of an inhomogeneous Poisson process }
\description{
  This function estimates the effects of covariates on individual-level risk using estimating equations with weights selected via data-driven methods in the situation where covariate information
 for controls is known only at a spatially aggregated level. Temporal trends in risk are estimated using a cubic spline. Optimisation
 is performed using the function \code{ nlminb()}, i.e. constrained optimisation using PORT routines. 
}
\usage{
weighted(Z, mu_t, weight_all, weight_data, beta0)
}

\arguments{
  \item{Z}{ A matrix of covariates, (dimension \code{n} by \code{p}) for the cases.  }
  \item{mu_t}{ A matrix of spatially aggregated covariates of dimension (\code{N} by \code{p}) by \code{T}. }
  \item{weight_all}{A \code{N} by \code{T} matrix of weights obtained from the function \code{find_weights()}}
  \item{weight_data}{ A vector of weights of length \code{n} obtained from the function \code{find_weights()} }
  \item{beta0}{ Initial values of the coefficients which is required to initialise the optimisation algorithm }
}
\details{
  Let \code{n} be the number of individual cases, \code{N} be the number of small areas within the study region, and \code{T} be the number of time points under consideration.
  The format of the matrix \code{mu_t} is such that the first N rows correspond to the spatially aggregated values of the first covariate for small areas \code{1} to \code{N}, and times \code{1} to \code{T}, and rows \code{N+1} to \code{2N} correspond to
  the spatially aggregated values of the second covariate for small areas \code{1} to \code{N}, and times \code{1} to \code{T} etc.
}
\value{
  Returns the list of components associated with the optimisation algorithm \code{ nlminb() }:
  \item{par }{The best set of parameters found.}
  \item{objective }{The value of the objective function corresponding to \code{par}.}
  \item{convergence}{An integer code. \code{0} indicates successful convergence.}
  
  \item{message}{A character string giving any additional information returned by the optimiser}
  
  \item{iterations}{Number of iterations performed}
  
  \item{evaluations}{Number of objective function evaluations}
}
\references{ The PORT documentation is at \url{ http://netlib.bell-labs.com/cm/cs/cstr/153.pdf }. 

Diggle, P.J., Guan, Y., Hart, C.A., Paize, F. 
and Stanton, M.C. (2010) Estimating 
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402
}
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}}
\note{ The \code{nlminb} routine is performed using \code{start=beta0, gradient=NULL, hessian=NULL, control=list(iter.max=1500,eval.max=1500),
                  lower=-100,upper=100}

}
\seealso{ \code{\link{nlminb}}, \code{\link{optim}}, \code{\link{unweighted}}, \code{\link{find_weights}}}
\examples{
data(alderhey, package="aggrisk")
 

attach(alderhey)

# Select Liverpool only

cases <- cases[which(cases$LA==2),]
OA.liv <- which(control[[5]][,1]!=0)


# Select age and sex as risk factors
# Assume there is no temporal trend in risk

X <- cbind(cases$age, cases$sex)


mu_t <- rbind(control$sumage[OA.liv,], control$sumgender[OA.liv,])

mu <- NULL
select <- 9:10
for(i in 1:2){
  sum.var <- apply(control[[select[i]]][OA.liv,],1,sum)
  mu <- cbind(mu, sum.var)
}

ncases <- dim(X)[1]
beta0 <- c(0,0)
optim_uw <- unweighted(X,mu, beta0)
beta_uw <- optim_uw$par


weights <- find_weights(beta_uw, mu_t, areaID=control$OAid[OA.liv],
caseID=cases$OA, caseT=cases$year, startT=1981,
endT=2007, pop=control$countLA2[OA.liv,])


optim_w <- weighted(X,mu_t,weights$weight_all,
+ weights$weight_data,beta0)
beta_w <- optim_w$par



\dontrun{
## This step may be time consuming ##
## Assume a temporal trend in risk ##

spline <- sterms(cases$year,length(OA.liv),1981,2007,
c(1987,1994,2000),control$countLA2[OA.liv,])


# Select age as a risk factor
X <-matrix(cases$age, nc=1)
Xtilde <- cbind(X, spline$cases)


mu_t <- control$sumage[OA.liv,]
mu_t.tilde <- rbind(mu_t, spline$control)

mu <- NULL
select <- 9
for(i in 1:1){
  sum.var <- apply(control[[select[i]]][OA.liv,],1,sum)
  mu <- cbind(mu, sum.var)
}

mu.tilde <- cbind(mu,spline$sumcontrol)
ncases <- dim(Xtilde)[1]
p <- dim(Xtilde)[2]
beta0 <- c(log(ncases)/sum(mu.tilde[,1]), rep(0.1,p-1))
optim_uw <- unweighted(Xtilde,mu.tilde, beta0)
beta_uw <- optim_uw$par


weights <- find_weights(beta_uw, mu_t.tilde, areaID=control$OAid[OA.liv],
caseID=cases$OA, caseT=cases$year, startT=1981,
endT=2007, pop=control$countLA2[OA.liv,])


optim_w <- weighted(Xtilde,mu_t.tilde,weights$weight_all,
+ weights$weight_data,beta0)
beta_w <- optim_w$par


}


}

\keyword{ spatial }
\keyword{ htest }
