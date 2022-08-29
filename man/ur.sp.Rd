\name{ur.sp}
\alias{ur.sp}
\encoding{latin1}
\title{Schmidt and Phillips Unit Root Test}
\description{
  Performs the Schmidt and Phillips unit root test, where under the Null
  and Alternative Hypothesis the coefficients of the deterministic
  variables are included.
}
\usage{
ur.sp(y, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4),
      signif = c(0.01, 0.05, 0.1))
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{type}{Test type, either \code{tau} or \code{rho} test.}
  \item{pol.deg}{Degree of polynomial in the test regression.}
  \item{signif}{Significance level for the critical value of the test statistic.}
}
\details{
  Under the Null and the Alternative hypothesis the coefficients of the
  deterministic part of the test regression are included. Two test types
  are available: the \code{rho}-test and the \code{tau}-test.
  Both test are extracted from the LM principle.
}
\value{
  An object of class \code{"ur.sp"}.
}
\references{
  Schmidt, P. and Phillips, P.C.B. (1992), LM Test for a Unit Root in
  the Presence of Deterministic Trends, \emph{Oxford Bulletin of Economics and
  Statistics}, \bold{54(3)}, 257--287.

  Download possible at: \url{https://cowles.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{\code{\link{ur.sp-class}}}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
sp.gnp <- ur.sp(gnp, type="tau", pol.deg=1, signif=0.01)
summary(sp.gnp)
}
\author{Bernhard Pfaff}
\keyword{regression}
