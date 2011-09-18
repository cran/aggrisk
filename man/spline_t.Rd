\name{spline_t}
\alias{spline_t}
\title{ A function to estimate the temporal variation in risk }
\description{
  This function obtains an estimate of the temporal trend in risk using regression
  coefficients obtained from estimating equations. The temporal trend is modelled
  using a cublic spline.
}
\usage{
spline_t(beta, startT, endT, knots, ref1, ref2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ Estimates of the regression coefficients }
  \item{startT}{ The start time of the study}
  \item{endT}{ The end time of the study }
  \item{knots}{ The locations of the knots used to fit the cubic spline }
  \item{ref1}{ The position of the coefficients of the fixed time terms in the vector \code{beta} }
  \item{ref2}{ The position of the coefficients of the spline terms in the vector \code{beta}  }
}
\details{
  There should be \code{3} fixed time term coefficients since a cubic spline is being
  used and \code{K} spline term coefficients, where \code{K} is the number of knots being used.
}
\value{
  Returns a list of the following components:
  \item{fixed}{A vector of the value of the fixed terms at each time point}
  \item{spline }{A vector of the value of the spline terms at each time point}
}
\references{Diggle, P.J., Guan, Y., Hart, C.A., Paize, F. 
and Stanton, M.C. (2010) Estimating 
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402}
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}  }

\seealso{ \code{\link{unweighted}, \link{weighted}} }
\examples{
data(alderhey, package="aggrisk")
 

attach(alderhey)

# Select Liverpool only

cases <- cases[which(cases$LA==2),]
OA.liv <- which(control[[5]][,1]!=0)


# Consider age and sex as risk factors only
# Include a temporal trend in risk

spline <- sterms(cases$year,length(OA.liv),1981,2007,
c(1987,1994,2000),control$countLA2[OA.liv,])

X <- cbind(cases$age, cases$sex)
Xtilde <- cbind(X, spline$cases)

mu_t <- rbind(control$sumage[OA.liv,], control$sumgender[OA.liv,])
mu_t.tilde <- rbind(mu_t, spline$control)

mu <- NULL
select <- 9:10
for(i in 1:2){
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


sp_fix <- spline_t(beta_w, startT=1981, endT=2007,knots=c(1987,1994,2000),
                   ref1=c(2,3,4), ref2=c(5,6,7))

plot(1981:2007,exp(sp_fix$fixed+sp_fix$spline), type="l", xlab="", ylab="exp[h(t)]")




}

\keyword{ spatial }
\keyword{ htest }
