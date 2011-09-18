\name{unweighted}
\alias{unweighted}
\title{This function uses estimating equations with equal weights to estimate regression parameters of an inhomogeneous Poisson process}
\description{
 This function estimates the effects of covariates on individual-level risk using unweighted estimating equations in the situation where covariate information
 for controls is known only at a spatially aggregated level. Temporal trends in risk are estimated using a cubic spline. Optimisation
 is performed using the function \code{ nlminb()}, i.e. constrained optimisation using PORT routines. 
}
\usage{
unweighted(Z, mu, beta0)
}

\arguments{
  \item{Z}{ A matrix of covariates, (dimension \code{n} by \code{p}) for the cases. }
  \item{mu}{ A matrix of spatially aggregated covariates summed over time corresponding to each of the small areas within the study region (dimension \code{N} by \code{p}). }
  \item{beta0}{ Initial values of the coefficients which is required to initialise the optimisation algorithm }
}
\details{
  Let \code{n} be the number of individual cases, \code{N} be the number of small areas within the study region and \code{p} be the number of covariates (including those related to the cubic spline).
}
\value{
  Returns the list of components associated with the optimisation algorithm \code{ nlminb()}:
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
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402 }
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}}

\note{
  The \code{nlminb} routine is performed using \code{start=beta0, gradient=NULL, hessian=NULL, control=list(iter.max=1500,eval.max=1500),
                  lower=-100,upper=100}
}
\seealso{ \code{\link{nlminb}}, \code{\link{optim}}, \code{\link{weighted}}, \code{\link{find_weights}} }
\examples{
  data(alderhey, package="aggrisk")
 

  attach(alderhey)  

# Select Liverpool only

cases <- cases[which(cases$LA==2),]
OA.liv <- which(control[[5]][,1]!=0)

# Consider age and sex as risk factors
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



\dontrun{
## This next step may be time consuming ##
## Assume a temporal trend in risk ##

spline <- sterms(cases$year,length(OA.liv),1981,2007,
c(1987,1994,2000),control$countLA2[OA.liv,])


# Consider age only as a risk factor
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


}
}


\keyword{ spatial }
\keyword{ htest }
