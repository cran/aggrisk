\name{testf}
\alias{testf}

\title{ A function to test the assumption that cases form an inhomogeneous Poisson process }
\description{
  This function calculates a test statistic with the purpose of assessing whether the assumption that the cases
form an inhomogeneous Poisson process is valid. 
}
\usage{
testf(controlid, caseid, counts, lambda, lags, xcoord, ycoord, h)
}
\arguments{
  \item{controlid}{ A vector of length \code{N} of small area identification variables }
  \item{caseid}{ A vector of length \code{n} of small area identification variables for the case data }
  \item{counts}{ Population of each small area summed over time }
  \item{lambda}{ The parametric part of the intensity function i.e. \code{exp(Z\%*\%beta)} }
  \item{lags}{ A vector of distances at which to calculate the test statistic }
  \item{xcoord}{ The x-coordinates of the centroids of each small area }
  \item{ycoord}{ The y-coordinates of the centroids of each small area }
  \item{h}{ bandwidth}
}
\details{
  \code{N} is the number of small areas in the study region, and \code{n} is the number of cases in the study.
}
\value{
 Returns a vector of the same length as the lags variable which corresponds to the value of the test statistic
at each of the lags.
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


lags=seq(1000,10000, by=1000)
R_w=exp(X\%*\%beta_w)
CENTRE_X=cent[OA.liv,2]
CENTRE_Y=cent[OA.liv,3]
count <- apply(control$countLA2[OA.liv,],1,sum)
pcf_w <- testf(controlid=control$OAid[OA.liv], caseid=cases$OA, counts=count,
lambda=R_w, lags=lags, xcoord=CENTRE_X, ycoord=CENTRE_Y, h=1000) # 9.89
plot(lags,pcf_w, type="l")



\dontrun{
## This step may be time consuming ##
## Include a temporal trend in risk ##

spline <- sterms(cases$year,length(OA.liv),1981,2007,
c(1987,1994,2000),control$countLA2[OA.liv,])

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


lags=seq(1000,10000, by=1000)
R_w=exp(Xtilde\%*\%beta_w)
CENTRE_X=cent[OA.liv,2]
CENTRE_Y=cent[OA.liv,3]
count <- apply(control$countLA2[OA.liv,],1,sum)
pcf_w <- testf(controlid=control$OAid[OA.liv], caseid=cases$OA, counts=count,
lambda=R_w, lags=lags, xcoord=CENTRE_X, ycoord=CENTRE_Y, h=1000) # 9.89
plot(lags,pcf_w, type="l")


}
}

\keyword{ spatial }
\keyword{ htest }
