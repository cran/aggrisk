\name{teststat}
\alias{teststat}

\title{A formal test statistic for checking the assumption that the location of cases form an inhomogeneous Poisson process

}
\description{
This function calculates a formal test statistic which can be used to check that the assumption 
that the location of the cases form an inhomogeneous Poisson process is valid. Assuming asymptomatic 
Normality of this test statistic, a Wald-type test can be applied.
}
\usage{
teststat(controlid, caseid, lambda, counts)
}
\arguments{
  \item{controlid}{ A vector of length \code{N} of small area identification variables }
  \item{caseid}{ A vector of length \code{n} of small area identification variables for the case data }
  \item{lambda}{ The parametric part of the intensity function i.e. \code{exp(Z\%*\%beta)} }
  \item{counts}{ Population of each small area summed over time }
}
\details{
\code{N} is the number of small areas in the study region, and \code{n} is the number of cases in the study.
}
\value{
Returns a scalar corresponding to the value of the formal test statistic
}
\references{
Diggle, P.J., Guan, Y., Hart, C.A., Paize, F. 
and Stanton, M.C. (2010) Estimating 
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402}
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}  }




\seealso{
\code{\link{var.test}}
}
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

test <- teststat(controlid=control$OAid[OA.liv], caseid=cases$OA,
lambda=exp(X\%*\%beta_w),counts=control$countLA2[OA.liv,])

var.test <- var_test(Z=X, lambda=exp(X\%*\%beta_w), beta=beta_w,
weights=weights$weight_data, counts=control$countLA2[OA.liv], mu=mu,
controlid=control$OAid[OA.liv], caseid=cases$OA)
p.val <- 2*pnorm(test/sqrt(var.test))


\dontrun{
## This step may be time consuming ##
## Include a temporal trend in risk ##

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


optim_w <- weighted(Xtilde,mu_t.tilde,weights$weight_all,
+ weights$weight_data,beta0)
beta_w <- optim_w$par

test <- teststat(controlid=control$OAid[OA.liv], caseid=cases$OA,
lambda=exp(Xtilde\%*\%beta_w),counts=control$countLA2[OA.liv,])

var.test <- var_test(Z=Xtilde, lambda=exp(Xtilde\%*\%beta_w), beta=beta_w,
weights=weights$weight_data, counts=control$countLA2[OA.liv], mu=mu.tilde,
controlid=control$OAid[OA.liv], caseid=cases$OA)
p.val <- 2*pnorm(test/sqrt(var.test))

}
}


\keyword{ spatial }
\keyword{ htest }
