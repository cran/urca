\name{lttest}
\alias{lttest}
\title{Likelihood ratio test for no linear trend in VAR}
\description{
  Conducts a likelihood ratio test for no inclusion of a linear trend in a
  VAR. That is, the Null hypothesis is for not including a linear trend
  and is assigned as 'H2*(r)'. The test statistic is distributed as
  \eqn{\chi^2} square with \eqn{(p-r)} degrees of freedom.
}
\usage{
lttest(z, r)
}
\arguments{
  \item{z}{An object of class `ca.jo'.}
  \item{r}{The count of cointegrating relationships.}
}
\details{
  The count of cointegrating relations should be given as integer and
  should be in the interval \eqn{1 \leq r < P}.
}
\value{
  \item{lttest}{Matrix containing the value of the test statistic and its
    p-value.}
}
\references{

  Johansen, S. and Juselius, K. (1990), Maximum Likelihood Estimation and
  Inference on Cointegration -- with Applications to the Demand for
  Money, \emph{Oxford Bulletin of Economics and Statistics}, \bold{52,
    2}, 169--210.

  Johansen, S. (1991), Estimation and Hypothesis Testing of
  Cointegration Vectors in Gaussian Vector Autoregressive Models,
  \emph{Econometrica}, \bold{Vol. 59, No. 6}, 1551--1580.

}
\seealso{
  \code{\link{ca.jo}} and \code{\link{ca.jo-class}}.
}
\examples{
data(denmark)
sjd <- as.matrix(denmark[, c("LRM", "LRY", "IBO", "IDE")])
sjd.vecm <- ca.jo(sjd, constant=TRUE, type="eigen", K=2, spec="longrun",
season=4)
lttest(sjd.vecm, r=1)
#
data(finland)
sjf <- as.matrix(finland)
sjf.vecm <- ca.jo(sjf, constant=FALSE, type="eigen", K=2,
spec="longrun", season=4)
lttest(sjf.vecm, r=3)
}
\author{Bernhard Pfaff}
\keyword{regression}
\concept{VECM Test Linear Trend Johansen Juselius}
