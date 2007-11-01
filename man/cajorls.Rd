\name{cajorls}
\alias{cajorls}
\encoding{latin1}
\title{OLS regression of VECM}
\description{
  This function returns the OLS regressions of a restricted VECM,
  \emph{i.e.} it returns a list object with elements of class `lm'
  containing the restricted VECM and a matrix object with the normalised
  cointegrating relationships. The user can provide a certain number of
  which equation in the VECM should be estimated and reported, or if
  \code{"reg.number = NULL"} each equation in the VECM will be estimated
  and its results are reported. Furthermore, the cointegratioon rank has
  to be supplied too.
}
\usage{
cajorls(z, r = 1, reg.number = NULL)
}
\arguments{
  \item{z}{An object of class `ca.jo' or `cajo.test'.}
  \item{r}{An integer, signifiying the cointegration rank.}
  \item{reg.number}{The number of the equation in the VECM that should
    be estimated or if set to \code{NULL} (the default), all equations
    within the VECM are estimated.}
}
\details{
 The cointegration space is normalised as \eqn{\bold{\beta}_c =
   \bold{\beta}(S'\bold{\beta})^{-1}}, with \eqn{S' = (I_r, 0)}. 
}
\value{
  Returns a list object with elements of class \code{`lm'} for the
  restricted VECM and a matrix object with the normalised cointegrating
  vectors. 
}
\references{

  Johansen, S. (1995), \emph{Likelihood-Based Inference in Cointegrated Vector
    Autoregressive Models}, Oxford University Press, Oxford.

  Lütkepohl, H. (2006), \emph{New Introduction to Multiple Time Series
    Analysis}, Springer, New York.

}
\seealso{
  \code{\link{ca.jo}}, \code{\link{cajools}}, \code{\link[base]{lm}},
  \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(finland)
sjf <- finland
sjf.vecm <- ca.jo(sjf, ecdet = "none", type = "eigen", K = 2,
spec = "longrun", season = 4)
sjf.vecm.rls <- cajorls(sjf.vecm, r = 2)
summary(sjf.vecm.rls$rlm)
sjf.vecm.rls$beta
}
\author{Bernhard Pfaff}
\keyword{regression}
\concept{VECM OLS Johansen Juselius Cointegration Co-integration}
