\name{alrtest}
\alias{alrtest}
\title{Likelihood ratio test for restrictions on alpha}
\description{
  This function estimates a restricted VAR, where the restrictions are
  base upon \eqn{\bold{\alpha}}, \emph{i.e.} the loading vectors. The test
  statistic is distributed as \eqn{\chi^2} with \eqn{r(p-m)} degrees of
  freedom, with \eqn{m} equal to the columns of the restricting matrix
  \eqn{\bold{A}}.
}
\usage{
alrtest(z, A, r)
}
\arguments{
  \item{z}{An object of class `ca.jo'.}
  \item{A}{The \eqn{(p \times m)} matrix containing the restrictions on
    \eqn{\bold{\alpha}}.}
  \item{r}{The count of cointegration relationships; \cr
    inferred from \code{summary(`ca.jo-object')}.}
}
\details{
  The orthogonal matrix to \eqn{\bold{A}} can be accessed as
  \code{object@B}. The restricted \eqn{\bold{\alpha}} matrix is
  normalised with respect to the first variable.
}
\value{
  An object of class \code{`cajo.test'}.
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
  \code{\link{ca.jo}}, \code{\link{blrtest}}, \code{\link{ablrtest}},
  \code{\link{cajo.test-class}}, \code{\link{ca.jo-class}} and
  \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, ecdet = "const", type="eigen", K=2, spec="longrun",
season=4)
DA <- matrix(c(1,0,0,0), c(4,1))
summary(alrtest(sjd.vecm, A=DA, r=1))
}
\author{Bernhard Pfaff}
\keyword{regression}
