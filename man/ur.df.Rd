\name{ur.df}
\alias{ur.df}
\encoding{latin1}
\title{Augmented-Dickey-Fuller Unit Root Test}
\description{
  Performs the augmented Dickey-Fuller unit root test. 
}
\usage{
ur.df(y, type = c("none", "drift", "trend"), lags = 1,
      selectlags = c("Fixed", "AIC", "BIC")) 
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{type}{Test type, either \code{"none"}, \code{"drift"}  or
    \code{"trend"}.}
  \item{lags}{Number of lags for endogenous variable to be included.}
  \item{selectlags}{Lag selection can be achieved according to the
    Akaike \code{"AIC"} or the Bayes \code{"BIC"} information
    criteria. The maximum number of lags considered is set by
    \code{lags}. The default is to use a \code{"fixed"} lag length set
    by \code{lags}.} 
}
\details{
  The function \code{ur.df()} computes the augmented Dickey-Fuller
  test. If type is set to \code{"none"} neither an intercept nor a trend
  is included in the test regression. If it is set to \code{"drift"} an
  intercept is added and if it is set to \code{"trend"} both an intercept
    and a trend is added. The critical values are taken from Hamilton
    (1994) and Dickey and Fuller(1981).
  }
\value{
  An object of class \code{ur.df}.
}
\references{
  Dickey, D. A. and Fuller, W. A. (1979), Distributions of the
  Estimators For Autoregressive Time Series with a Unit Root,
  \emph{Journal of the American Statistical Association},  \bold{75},
  427--431.

  Dickey, D. A. and Fuller, W. A. (1981), Likelihood Ratio Statistics
  for Autoregressive Time Series with a Unit Root, \emph{Econometrica},
  \bold{49}, 1057--1072.
  
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
