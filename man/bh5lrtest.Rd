\name{bh5lrtest}
\alias{bh5lrtest}
\title{Likelihood ratio test for restrictions under partly known beta}
\description{
  This function estimates a restricted VAR, where some of the
  cointegration vectors are known. The known cointegration relationships
  have to be provided in an \eqn{p x r1} matrix \eqn{\bold{H}}. The test
  statistic is distributed as \eqn{\chi^2} with \eqn{(p-r)r1} degrees of
  freedom, with \eqn{r} equal to total number of cointegration relations.
}
\usage{
bh5lrtest(z, H, r)
}
\arguments{
  \item{z}{An object of class `ca.jo'.}
  \item{H}{The \eqn{(p \times r1)} matrix containing the known
    cointegration relations.}
  \item{r}{The count of cointegrating relationships; \cr
    inferred from \code{summary(`ca.jo-object')}.}
}
\details{
  Please note, that the number of columns of \eqn{\bold{H}} must be
  smaller than the count of cointegration relations \eqn{r}.
}
\value{
  An object of class \code{`cajo.test'}.
}
\references{

Johansen, S. (1995), \emph{Likelihood-Based Inference in Cointegrated Vector
Autoregressive Models}, Oxford University Press, Oxford.

Johansen, S. and Juselius, K. (1992), Testing structural hypotheses in a
multivariate cointegration analysis of the PPP and the UIP for UK,
\emph{Journal of Econometrics}, \bold{53}, 211--244.

}
\seealso{
  \code{\link{ca.jo}}, \code{\link{alrtest}}, \code{\link{ablrtest}},
  \code{\link{blrtest}}, \code{\link{bh6lrtest}}, \code{\link{cajo.test-class}},
  \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(UKpppuip)
attach(UKpppuip)
dat1 <- cbind(p1, p2, e12, i1, i2)
dat2 <- cbind(doilp0, doilp1)
H1 <- ca.jo(dat1, type='trace', K=2, season=4, dumvar=dat2)
H51 <- c(1, -1, -1, 0, 0)
H52 <- c(0, 0, 0, 1, -1)
summary(bh5lrtest(H1, H=H51, r=2))
summary(bh5lrtest(H1, H=H52, r=2))
}
\author{Bernhard Pfaff}
\keyword{regression}
