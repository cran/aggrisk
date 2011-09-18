\name{var_beta}
\alias{var_beta}
\title{ Function to estimate the variance of the regression coefficients}
\description{
  This function estimates the variance of the regression coefficients in the intensity function of an 
inhomogeneous spatial Poisson process, using individual-level case data spatially aggregated data for 
the population at risk. The regression coefficients are estimated using estimating equations.
}
\usage{
var_beta(beta, Z, weights)
}
\arguments{
  \item{beta}{ A vector of estimates of the regression coefficients }
  \item{Z}{ An \code{n} by \code{p} matrix of covariates relating to the cases }
  \item{weights}{ A vector of length \code{n} of weights used to fit the model }
}
\details{
  \code{n} is the number of cases in the study, and \code{p} is the number of covariates (including those related to the
  cubic spline).
}
\value{
The \code{p} by \code{p} covariance matrix corresponding to the regression coefficients
}
\references{Diggle, P.J., Guan, Y., Hart, C.A., Paize, F. 
and Stanton, M.C. (2010) Estimating 
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402}
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}  }

\seealso{ \code{\link{unweighted}}, \code{\link{weighted}} }
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

weights <- find_weights(beta_uw, mu_t, areaID=control$OAid[OA.liv],
caseID=cases$OA, caseT=cases$year, startT=1981,
endT=2007, pop=control$countLA2[OA.liv,])


optim_w <- weighted(X,mu_t,weights$weight_all,
+ weights$weight_data,beta0)
beta_w <- optim_w$par


var_bw <- var_beta(beta_w, X, weights=weights$weight_data)


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


var_bw <- var_beta(beta_w, Xtilde, weights=weights$weight_data)

}
}

\keyword{ spatial }
\keyword{ htest }
