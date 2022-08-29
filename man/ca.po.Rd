\name{ca.po}
\alias{ca.po}
\encoding{latin1}
\title{Phillips and Ouliaris Cointegration Test}
\description{
  Performs the Phillips and Ouliaris \code{"Pu"} and \code{"Pz"}
  cointegration test.
}
\usage{
ca.po(z, demean = c("none", "constant", "trend"),
      lag = c("short", "long"), type = c("Pu", "Pz"), tol = NULL)
}
\arguments{
  \item{z}{Data matrix to be investigated for cointegration.}
  \item{demean}{The method for detrending the series, either
    \code{"none"}, \code{"constant"} or \code{"trend"}.}
  \item{lag}{Either a short or long lag number used for
    variance/covariance correction.}
  \item{type}{The test type, either \code{"Pu"} or \code{"Pz"}.}
  \item{tol}{Numeric, this argument is passed to \code{solve()} in
    \code{ca.po()}.}
}
\details{
  The test \code{"Pz"}, compared to the test \code{"Pu"}, has the
  advantage that it is invariant to the normalization of the
  cointegration vector, \emph{i.e.} it does not matter which variable
  is on the left hand side of the equation. In case convergence
  problems are encountered by matrix inversion, one can pass a higher
  tolerance level \emph{via} \code{"tol=..."} to the \code{solve()}-function.
}
\value{
  An object of class \code{ca.po}.
}
\references{
  Phillips, P.C.B. and Ouliaris, S. (1990), Asymptotic Properties of
  Residual Based Tests for Cointegration, \emph{Econometrica},
  \bold{Vol. 58, No. 1}, 165--193.
}
  
\seealso{\code{\link{ca.po-class}}}
\examples{
data(ecb)
m3.real <- ecb[,"m3"]/ecb[,"gdp.defl"]
gdp.real <- ecb[,"gdp.nom"]/ecb[,"gdp.defl"]
rl <- ecb[,"rl"]
ecb.data <- cbind(m3.real, gdp.real, rl)
m3d.po <- ca.po(ecb.data, type="Pz")
summary(m3d.po)
}
\author{Bernhard Pfaff}
\keyword{regression}
