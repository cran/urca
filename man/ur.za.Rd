\name{ur.za}
\alias{ur.za}
\encoding{latin1}
\title{Zivot \& Andrews Unit Root Test}
\description{
  Performs the Zivot \& Andrews unit root test, which allows a
  break at an unknown point in either the intercept, the linear
  trend or in both.
}
\usage{
ur.za(y, model = c("intercept", "trend", "both"), lag=NULL)
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{model}{Specification if the potential break occured in either the intercept, the linear trend or in both.}
  \item{lag}{The highest number of lagged endogenous differenced variables to be included in the test regression}
}
\details{
  This test is based upon the recursive estimation of a test
  regression. The test statistic is defined as the minimum t-statistic
  of the coeffcient of the lagged endogenous variable. 
}
\value{
  An object of class \code{ur.za}.
}
\references{
  Zivot, E. and Andrews, Donald W.K. (1992), Further Evidence on the
  Great Crash, the Oil-Price Shock, and the Unit-Root Hypothesis,
  \emph{Journal of Business \& Economic Statistics}, \bold{10(3)},
  251--270.

  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{\code{\link{ur.za-class}}}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
za.gnp <- ur.za(gnp, model="both", lag=2)
summary(za.gnp)
}
\author{Bernhard Pfaff}
\keyword{regression}
