\name{cajools}
\alias{cajools}
\title{OLS regression of VECM}
\description{
  This function returns the OLS regressions of an unrestricted VECM,
  \emph{i.e.} it returns an object of class `lm'. The user can provide a
  certain number of which equation in the VECM should be estimated and
  reported, or if \code{"reg.number=NULL"} each equation in the VECM
  will be estimated and its results are reported.
}
\usage{
cajools(z, reg.number = NULL)
}
\arguments{
  \item{z}{An object of class `ca.jo' or `cajo.test'.}
  \item{reg.number}{The number of the equation in the VECM that should
    be estimated or if set to \code{NULL} (the default), all equations
    within the VECM are estimated.}
}
\value{
  Returns an object of class \code{`lm'}.
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
  \code{\link{ca.jo}}, \code{\link[base]{lm}}, \code{\link{ca.jo-class}}
  and \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, constant=TRUE, type="eigen", K=2, season=4)
HD1 <- matrix(c(1, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1), c(5,3))
DA <- matrix(c(1,0,0,0, 0, 1, 0, 0, 0, 0, 0, 1), c(4,3))
sjd.vecm.res <- ablrtest(sjd.vecm, H=HD1, A=DA, r=1)
sjd.vecm.ols <- cajools(sjd.vecm.res)
summary(sjd.vecm.ols)
sjd.vecm.ols1 <- cajools(sjd.vecm.res, reg.number=1)
summary(sjd.vecm.ols1)
}
\keyword{regression}
