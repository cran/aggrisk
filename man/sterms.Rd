\name{sterms}
\alias{sterms}
\title{ A function to obtain the covariates relating to the cubic spline }
\description{
  This function obtains the fixed terms and the spline terms relating to a cubic spline for both the individual-level
case data and spatially aggregated control data 
}
\usage{
sterms(caseT, N, startT, endT, knots, pop)
}

\arguments{
  \item{caseT}{ A vector of length \code{n} of the times of event for each case }
  \item{N}{ The number of small areas within the study region }
  \item{startT}{The start time of the study}
  \item{endT}{ The end time of the study}
  \item{knots}{ The location of each of the \code{K} knots (\code{K} determined by user) }
  \item{pop}{A vector of the population of each of the small areas aggregated over time}
}
\details{
  Let \code{n} be the number of individual cases
}
\value{
  Returns the following list of components:
  \item{cases }{An \code{n} by (\code{K+3}) matrix, where the first \code{3} columns are the fixed components and the following \code{K} columns
                are the spline components of the cubic spline relating to the cases}
  \item{control }{An \code{N} by (\code{K+3}) by \code{T} matrix, where the first (\code{3} by \code{N}) rows are the fixed components and the following (\code{K} by \code{N}) rows
                are the spline components of the cubic spline relating to each small area over time}
  \item{sumcontrol}{An \code{N} by (\code{K+3}) matrix where first \code{3} columns are the fixed components and the following \code{K} columns
                are the spline components of the cubic spline relating to each small area summed over time}
}  
\references{Diggle, P.J., Guan, Y., Hart, C.A., Paize, F. 
and Stanton, M.C. (2010) Estimating 
individual-level risk in spatial epidemiology using spatially aggregated information on the population at risk. \emph{Journal of the American Statistical Association}. \bold{105}:1394--1402}
\author{ Michelle Stanton \email{m.stanton@lancaster.ac.uk}, Yongtao Guan \email{yongtao.guan@yale.edu}  }

\seealso{ \code{\link{unweighted}}}
\examples{
data(alderhey, package="aggrisk")
 

attach(alderhey)

# Calculate the spline terms

spline <- sterms(cases$year,dim(control$sumall)[1],1981,2007,c(1987,1994,2000),
                 control$sumall)


}

\keyword{ spatial }
\keyword{ htest }
