\name{ablrtest}
\alias{ablrtest}
\encoding{latin1}
\title{Likelihood ratio test for restrictions on alpha and beta}
\description{
  This function estimates a restricted VAR, where the restrictions are
  based upon \eqn{\bold{\alpha}}, \emph{i.e.} the loading vectors and
  \eqn{\bold{\beta}}, \emph{i.e} the matrix of cointegration vectors. The test
  statistic is distributed as \eqn{\chi^2} with \eqn{(p-m)r + (p-s)r} degrees of
  freedom, with \eqn{m} equal to the columns of the restricting matrix
  \eqn{\bold{A}}, \eqn{s} equal to the columns of the restricting matrix
  \eqn{\bold{H}} and \eqn{p} the order of the VAR.
}
\usage{
ablrtest(z, H, A, r)
}
\arguments{
  \item{z}{An object of class \code{ca.jo}.}
  \item{H}{The \eqn{(p \times s)} matrix containing the restrictions on
    \eqn{\bold{\beta}}.}
  \item{A}{The \eqn{(p \times m)} matrix containing the restrictions on
    \eqn{\bold{\alpha}}.}
  \item{r}{The count of cointegrating relationships; \cr
    inferred from \code{summary(ca.jo-object)}.}
}
\details{
  The restricted \eqn{\bold{\alpha}} matrix, as well as \eqn{\bold{\beta}} is
  normalised with respect to the first variable.
}
\value{
  An object of class \code{cajo.test}.
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
  \code{\link{ca.jo}}, \code{\link{alrtest}}, \code{\link{blrtest}},
  \code{\link{cajo.test-class}}, \code{\link{ca.jo-class}} and
  \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, ecdet = "const", type="eigen", K=2, spec="longrun",
season=4)
HD1 <- matrix(c(1, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1), c(5,3))
DA <- matrix(c(1,0,0,0, 0, 1, 0, 0, 0, 0, 0, 1), c(4,3))
summary(ablrtest(sjd.vecm, H=HD1, A=DA, r=1))
}
\author{Bernhard Pfaff}
\keyword{regression}
