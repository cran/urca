\name{blrtest}
\alias{blrtest}
\encoding{latin1}
\title{Likelihood ratio test for restrictions on beta}
\description{
  This function estimates a restricted VAR, where the restrictions are
  base upon \eqn{\bold{\beta}}, \emph{i.e.} the cointegration vectors. The test
  statistic is distributed as \eqn{\chi^2} with \eqn{r(p-s)} degrees of
  freedom, with \eqn{s} equal to the columns of the restricting matrix
  \eqn{\bold{H}}.
}
\usage{
blrtest(z, H, r)
}
\arguments{
  \item{z}{An object of class \code{ca.jo}.}
  \item{H}{The \eqn{(p \times s)} matrix containing the restrictions on
    \eqn{\bold{\beta}}.}
  \item{r}{The count of cointegrating relationships; \cr
    inferred from \code{summary(ca.jo-object)}.}
}
\details{
  Please note, that in the case of nested hypothesis, the reported
  p-value should be adjusted to \eqn{r(s1-s2)} (see Johansen, S. and
  K. Juselius (1990)).
}
\value{
  An object of class \code{cajo.test}.
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
  \code{\link{ca.jo}}, \code{\link{alrtest}}, \code{\link{ablrtest}},
  \code{\link{bh5lrtest}}, \code{\link{bh6lrtest}}, \code{\link{cajo.test-class}},
  \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, ecdet="const", type="eigen", K=2, spec="longrun",
season=4)
HD0 <- matrix(c(-1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), c(5,4))
summary(blrtest(sjd.vecm, H=HD0, r=1))
}
\author{Bernhard Pfaff}
\keyword{regression}
