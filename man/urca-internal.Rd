\name{.spcv}
\alias{.spcv}
\title{Critical values for Schmidt \& Phillips Unit Root Test}
\description{
  This function is an internal function and is called by
  \code{ur.sp}. It computes the critical value of the Schmidt \&
  Phillips test, given a level of significance, the polynomial degree of
  the test regression, the test type and the sample size.
}
\usage{
.spcv(obs, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4),
      signif = c(0.01, 0.025, 0.05, 0.1))
}
\arguments{
  \item{obs}{The sample size.}
  \item{type}{The test type.}
  \item{pol.deg}{The polynomial degree.}
  \item{signif}{The significance level.}
}
\value{
  The critical value of the test.
}
\references{
  Schmidt, P. and Phillips, P.C.B. (1992), LM Test for a Unit Root in
  the Presence of Deterministic Trends, \emph{Oxford Bulletin of Economics and
  Statistics}, \bold{54(3)}, 257--287.

  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{\code{\link{ur.sp}}.}
\author{Bernhard Pfaff}
\keyword{regression}
