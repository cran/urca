\name{ur.pp}
\alias{ur.pp}
\title{Phillips \& Perron Unit Root Test}
\description{
  Performs the Phillips \& Perron unit root test. Beside the Z
  statistics Z-alpha and Z-tau, the Z statistics for the deterministic
  part of the test regression are computed, too.
}
\usage{
ur.pp(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"),
      lags = c("short", "long"), use.lag = NULL)
}
\arguments{
  \item{x}{Vector to be tested for a unit root.}
  \item{type}{Test type, either \code{"Z-alpha"} or \code{"Z-tau"}.}
  \item{model}{Determines the deterministic part in the test regression.}
  \item{lags}{Lags used for correction of error term.}
  \item{use.lag}{Use of a different lag number, specified by the user.}
}
\details{
  The function \code{ur.pp()} computes the Phillips \& Perron test. For
  correction of the error term a Bartlett window is used. 
}
\value{
  An object of class \code{`ur.pp'}.
}
\references{
  Phillips, P.C.B. and Perron, P. (1988), Testing for a unit root in
  time series regression, \emph{Biometrika},  \bold{75(2)}, 335--346.

  MacKinnon, J.G. (1991), Critical Values for Cointegration Tests,
  \emph{Long-Run Economic Relationships}, eds. R.F. Engle and
  C.W.J. Granger, London, Oxford, 267--276.
  
  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)' and
  \url{http://www.econ.ucsd.edu/papers/files/90-4.pdf}.
}
\seealso{\code{\link{ur.pp-class}}.}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
pp.gnp <- ur.pp(gnp, type="Z-tau", model="trend", lags="short")
summary(pp.gnp)
}
\author{Bernhard Pfaff}
\keyword{regression}
