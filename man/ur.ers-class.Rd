\name{ur.ers-class}
\docType{class}
\alias{ur.ers-class}
\title{Representation of class `ur.ers'}
\description{
  This class contains the relevant information by applying the Elliott,
  Rothenberg \& Stock unit root test.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"vector"}: The time series to
      be tested.}
    \item{\code{yd}:}{Object of class \code{"vector"}: The detrended
      time series.}
    \item{\code{type}:}{Object of class \code{"character"}: Test type,
      either \code{"DF-GLS"} (default), or \code{"P-test"}.}
    \item{\code{model}:}{Object of class \code{"character"}: The
      deterministic model used for detrending, either intercept only, or
      intercept with linear trend.}
    \item{\code{lag}:}{Object of class \code{"integer"}: The number of
      lags used in the test/auxiliary regression.}
    \item{\code{cval}:}{Object of class \code{"matrix"}: The critical
      values of the test at the 1\%, 5\% and 10\% level of significance.}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: The value
      of the test statistic.}
    \item{\code{testreg}:}{Object of class \code{"ANY"}: The test
      regression, only set for \code{"DF-GLS"}.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e.} `Elliott, Rothenberg \& Stock'.}
  }
}
\section{Extends}{
Class \code{`urca'}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ur.ers")} at the R prompt for a
  complete list of methods which are available for this class.
  
  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic.}
    \item{\code{summary}:}{like show, but test type, test regression (\code{type="DF-GLS"}) and critical values added.}
    \item{\code{plot}:}{Diagram of fit, residual plot and their acfs'
      and pacfs' for \code{type="DF-GLS"}.}
  }

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
\seealso{
  \code{\link{ur.ers}} and \code{\link{urca-class}}.
}
\author{Bernhard Pfaff}
\keyword{classes}
