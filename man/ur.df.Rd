\name{ur.df}
\alias{ur.df}
\title{Augmented-Dickey-Fuller Unit Root Test}
\description{
  Performs the augmented Dickey-Fuller unit root test. 
}
\usage{
ur.df <- function (y, type = c("none", "drift", "trend"), lags = 1) 
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{type}{Test type, either \code{"none"}, \code{"drift"}  or
    \code{"trend"}.}
  \item{lags}{Number of lags for endogenous variable to be included.}
}
\details{
  The function \code{ur.df()} computes the augmented Dickey-Fuller
  test. If type is set to \code{"none"} neither an intercept nor a trend
  is included in the test regression. If it is set to \code{"drift"} an
  intercept is added and if it is set to \code{"trend"} both an intercept
    and a trend is added. The critical values are taken from Hamilton (1994).
  }
\value{
  An object of class \code{`ur.df'}.
}
\references{
  Dickey, D. A. and Fuller, W. A. (1979), Distributions of the
  Estimators For Autoregressive Time Series with a Unit Root,
  \emph{Journal of the American Statistical Association},  \bold{75}, 427--431.

  Hamilton (1994), \emph{Time Series Analysis}, Princeton University Press.
}
\seealso{\code{\link{ur.df-class}}.}
\examples{
data(Raotbl3)
attach(Raotbl3)
lc.df <- ur.df(y=lc, lags=3, type='trend')
summary(lc.df)
}
\author{Bernhard Pfaff}
\keyword{regression}
