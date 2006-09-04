\name{bh6lrtest}
\alias{bh6lrtest}
\title{Likelihood ratio test for restrictions under partly known beta in
a subspace}
\description{
  This function estimates a restricted VAR, where some restrictions are
  placed on \eqn{r1} cointegrating relations which are chosen in the
  space of the matrix H. The test statistic is distributed as
  \eqn{\chi^2} with \eqn{(p-s-r2)r1} degrees of freedom, with \eqn{s}
  equal to the number of columns of \eqn{\bold{H}}, \eqn{r1} the number
  of cointegrating relations in the first partition and \eqn{r2} the
  number of cointegrating relations in the second partition which will
  be estimated without any restrictions.
}
\usage{
bh6lrtest(z, H, r, r1, conv.val=0.0001, max.iter=50)
}
\arguments{
  \item{z}{An object of class `ca.jo'.}
  \item{H}{The \eqn{(p \times s)} matrix containing the known
    cointegration relations.}
  \item{r}{The count of cointegrating relationships; \cr
    inferred from \code{summary(`ca.jo-object')}.\cr}
  \item{r1}{The count of cointegrating relationships in the first
    partition of the cointegration space; \cr}
  \item{conv.val}{The convergence value of the algorithm. (see details); \cr}
  \item{max.iter}{The maximal number of iterations.}
    
  
}
\details{
  Please note, that the following ordering of the dimensions should be
  obeyed: \eqn{r1 \leq s \leq p - r2}. A two-step algorithm is used to
  determine the eigen values of the restricted model. Convergence is
  achieved if the quadratic norm of the eigen values is smaller than
  \code{conv.val}.
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
  \code{\link{blrtest}}, \code{\link{bh5lrtest}}, \code{\link{cajo.test-class}},
  \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(UKpppuip)
attach(UKpppuip)
dat1 <- cbind(p1, p2, e12, i1, i2)
dat2 <- cbind(doilp0, doilp1)
H1 <- ca.jo(dat1, type='trace', K=2, season=4, dumvar=dat2, ctable='A1')
H6 <- matrix(c(1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0), c(5,3))
bh6lrtest(z=H1, H=H6, r=2, r1=1, conv.val=0.0001, max.iter=50)
}
\author{Bernhard Pfaff}
\keyword{regression}
