\name{find_weights}
\alias{find_weights}

\title{ A function to find the weights for the estimating equations used to estimate the
        regression parameters of an inhomogeneous Poisson process}
\description{
  This function implements data-driven methods to select the weights for a set of estimating equations.
  These estimating equations will be used to estimate the effects of covariates on individual-level risk
  in a situation where covariates for controls are known only at the spatially
  aggregated level.
}
\usage{
find_weights(beta_uw, mu_t, areaID, caseID, caseT, startT, endT, pop)
}

\arguments{
  \item{beta_uw}{ Estimates of the regression coefficients obtained using unweighted estimating equations }
  \item{mu_t}{ A matrix of spatially aggregated covariates of dimension (\code{N} by \code{p}) by \code{T}. }
  \item{areaID}{ A vector of identification variables for each small area in the study }
  \item{caseID}{ A vector of small area identification variables corresponding to each case }
  \item{caseT}{ A vector of times corresponding to each case }
  \item{startT}{ The start time of the study }
  \item{endT}{ The end time of the study }
  \item{pop}{ An \code{N} by \code{T} matrix of population counts }
}
\details{
   Let \code{N} be the number of small areas within the study region, \code{T} the number of time points under consideration, and \code{p} the number of covariates (see example).
  The format of the matrix \code{mu_t} is such that the first \code{N} rows correspond to the spatially aggregated values of the first covariate for small areas \code{1} to \code{N}, and times \code{1} to \code{T}, and rows \code{N+1} to \code{2N} correspond to
  the spatially aggregated values of the second covariate for small areas \code{1} to \code{N}, and times \code{1} to \code{T} etc.
}
\value{
  Returns a list of the following components:
  \item{weight_all }{A \code{N} by \code{T} matrix of weights for each small area at each time point}
  \item{weight_data }{A vector of \code{n} weights corresponding to each case at time of event}
}
\references{Diggle, P.J., Guan, Y., Hart, C.A., Paize, F. 
and Stanton, M.C. (2010) Estimating 
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402 }
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}  }

\seealso{ \code{\link{unweighted}, \link{weighted}} }
\examples{
data(alderhey, package="aggrisk")
 

attach(alderhey)

# Select Liverpool only

cases <- cases[which(cases$LA==2),]
OA.liv <- which(control[[5]][,1]!=0)

# Consider age and sex as risk factors. Note we are not considering
# a temporal trend in risk in this example
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


\dontrun{

## May be time consuming ##

# Include a temporal trend in risk
spline <- sterms(cases$year,length(OA.liv),1981,2007,
c(1987,1994,2000),control$countLA2[OA.liv,])

# Consider age as a risk factor only
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

}

}

\keyword{ spatial }
\keyword{ htest }
