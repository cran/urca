\name{ur.kpss-class}
\docType{class}
\alias{ur.kpss-class}
\encoding{latin1}
\title{Representation of class ur.kpss}
\description{
  This class contains the relevant information by applying the
  Kwiatkowski, Phillips, Schmidt \& Shin unit root test to a time series.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"vector"}: The time series to
      be tested.}
    \item{\code{type}:}{Object of class \code{"character"}: Test type,
      \code{"mu"} or \code{"tau"} depending on the deterministic part.}
    \item{\code{lag}:}{Object of class \code{"integer"}: Number of lags
      for error term correction.}
    \item{\code{cval}:}{Object of class \code{"matrix"}: Critical value
      of test.}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: Value of
      test statistic.}
    \item{\code{res}:}{Object of class \code{"vector"}: Residuals of
      test regression.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e.} `KPSS'.}
  }
}
\section{Extends}{
Class \code{urca}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ur.kpss")} at the R prompt for a
  complete list of methods which are available for this class.
  
  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic.}
    \item{\code{summary}:}{like show, but critical values, lags and test
      type added.}
    \item{\code{plot}:}{Residual plot and their acfs' and pacfs'.}
  }

}
\references{
  Kwiatkowski, D., Phillips, P.C.B., Schmidt, P. and Shin, Y., (1992),
  Testing the Null Hypothesis of Stationarity Against the Alternative of
  a Unit Root: How Sure Are We That Economic Time Series Have a Unit
  Root?, \emph{Journal of Econometrics}, \bold{54}, 159--178.

  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{
  \code{\link{ur.kpss}} and \code{\link{urca-class}}.
}
\author{Bernhard Pfaff}
\keyword{classes}
