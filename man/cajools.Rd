\name{cajools}
\alias{cajools}
\encoding{latin1}
\title{OLS regression of VECM}
\description{
  This function returns the OLS regressions of an unrestricted VECM,
  \emph{i.e.} it returns an object of class \code{lm}. The user can provide a
  certain number of which equation in the VECM should be estimated and
  reported, or if \code{"reg.number=NULL"} each equation in the VECM
  will be estimated and its results are reported.
}
\usage{
cajools(z, reg.number = NULL)
}
\arguments{
  \item{z}{An object of class \code{ca.jo} or \code{cajo.test}.}
  \item{reg.number}{The number of the equation in the VECM that should
    be estimated or if set to \code{NULL} (the default), all equations
    within the VECM are estimated.}
}
\value{
  Returns an object of class \code{lm}.
}
\references{

  Johansen, S. (1988), Statistical Analysis of Cointegration Vectors,
  \emph{Journal of Economic Dynamics and Control}, \bold{12}, 231--254.

  Johansen, S. and Juselius, K. (1990), Maximum Likelihood Estimation and
  Inference on Cointegration -- with Applications to the Demand for
  Money, \emph{Oxford Bulletin of Economics and Statistics}, \bold{52,
    2}, 169--210.

  Johansen, S. (1991), Estimation and Hypothesis Testing of
  Cointegration Vectors in Gaussian Vector Autoregressive Models,
  \emph{Econometrica}, \bold{Vol. 59, No. 6}, 1551--1580.

}
\seealso{
  \code{\link{ca.jo}}, \code{\link{cajorls}}, \code{\link[base]{lm}},
  \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm1 <- ca.jo(sjd, ecdet = "const", type="eigen", K=2, spec="longrun",
season=4)
sjd.vecm2 <- ca.jo(sjd, ecdet = "const", type="eigen", K=2, spec="transitory",
season=4)
sjd.vecm.ols1 <- cajools(sjd.vecm1)
sjd.vecm.ols2 <- cajools(sjd.vecm2)
summary(sjd.vecm.ols1)
summary(sjd.vecm.ols2)
}
\author{Bernhard Pfaff}
\keyword{regression}
\concept{VECM OLS Johansen Juselius Cointegration Co-integration}
