\name{alphaols}
\alias{alphaols}
\title{OLS regression of VECM weighting matrix}
\description{
  This functions estimates the \eqn{\bold{\alpha}} matrix of a VECM.
  The following OLS regression of the R-form of the VECM is hereby
  utilised:
  \deqn{\bold{R}_{0t} = \bold{\alpha}\bold{\beta}\prime \bold{R}_{kt} +
    \bold{\varepsilon}_t \qquad t=1, \dots, T}
}
\usage{
alphaols(z, reg.number = NULL)
}
\arguments{
  \item{z}{An object of class `ca.jo'.}
  \item{reg.number}{The number of the equation in the R-form that should
    be estimated or if set to \code{NULL} (the default), all equations
    within the R-form are estimated.}
}
\details{
  The cointegrating relations, \emph{i.e.} \eqn{\bold{R}_{kt}\prime
    \bold{\beta}} are calculated by using \code{z@RK} and \code{z@V}.
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
sjd.vecm1 <- ca.jo(sjd, constant=TRUE, type="eigen", K=2, spec="longrun",
season=4)
summary(alphaols(sjd.vecm1))
summary(alphaols(sjd.vecm1, reg.number=1))
}
\author{Bernhard Pfaff}
\keyword{regression}
\concept{VECM OLS Loading Johansen Juselius Cointegration Co-integration}
