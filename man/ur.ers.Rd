\name{ur.ers}
\alias{ur.ers}
\title{Elliott, Rothenberg \& Stock Unit Root Test}
\description{
  Performs the Elliott, Rothenberg \& Stock unit root test. 
}
\usage{
ur.ers(y, type = c("DF-GLS", "P-test"), model = c("constant", "trend"),
       lag.max = 4)
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{type}{Test type, either \code{"DF-GLS"} (default), or \code{"P-test"}.}
  \item{model}{The deterministic model used for detrending.}
  \item{lag.max}{The maximum numbers of lags used for testing of a
    decent lag truncation for the \code{"P-test"} (\code{BIC} used),
    or the maximum number of lagged differences to be included in the test
    regression for \code{"DF-GLS"}.}
}
\details{
  To improve the power of the unit root test, Elliot, Rothenberg \& Stock
  proposed a local to unity detrending of the time series. ERS developed
  a feasible point optimal test, \code{"P-test"}, which takes serial
  correlation of the error term into account. The second test type is
  the \code{"DF-GLS"} test, which is an ADF-type test applied to the
  detrended data without intercept. Critical values for this test are
  taken from MacKinnon in case of \code{model="constant"} and else from
  Table 1 of Elliot, Rothenberg \& Stock.
}
\value{
  An object of class \code{"ur.ers"}.
}

\references{
  Elliott, G., Rothenberg, T.J. and Stock, J.H. (1996), Efficient Tests
  for an Autoregressive Unit Root, \emph{Econometrica},
  \bold{Vol. 64, No. 4}, 813--836.

  MacKinnon, J.G. (1991), Critical Values for Cointegration Tests,
  \emph{Long-Run Economic Relationships}, eds. R.F. Engle and
  C.W.J. Granger, London, Oxford, 267--276.

  Download possible at: \url{http://www.econ.ucsd.edu/papers/files/90-4.pdf}.
}

\seealso{\code{\link{ur.ers-class}}}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
ers.gnp <- ur.ers(gnp, type="DF-GLS", model="const", lag.max=4)
summary(ers.gnp)
}
\author{Bernhard Pfaff}
\keyword{regression}
