\name{plotres}
\alias{plotres}
\encoding{latin1}
\title{Graphical inspection of VECM residuals}
\description{
  The function \code{plotres} should be used for graphical inspection
  of the VAR residuals, \emph{i.e.} the estimated specification as
  elaborated in the `Details' section of \code{ca.jo}. It displays the
  residuals for each equation within a VAR and their acf's and pacf's.
}
\usage{
plotres(x)
}
\arguments{
  \item{x}{Object of class `ca.jo'.}
}
\references{
  Johansen, S. and Juselius, K. (1990), Maximum Likelihood Estimation and
  Inference on Cointegration -- with Applications to the Demand for
  Money, \emph{Oxford Bulletin of Economics and Statistics}, \bold{52,
    2}, 169--210.
  
}
\seealso{
  \code{\link{ca.jo}} and \code{\link{ca.jo-class}}.
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, ecdet="const", type="eigen", K=2, spec="longrun",
season=4)
plotres(sjd.vecm)
}
\author{Bernhard Pfaff}
\keyword{regression}
\concept{VECM Residuals Johansen Juselius}
