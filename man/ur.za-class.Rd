\name{ur.za-class}
\docType{class}
\alias{ur.za-class}

\title{Representation of class `ur.za'}
\description{
  This class contains the relevant information by applying the Zivot \& Andrews
  unit root test to a time series.}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"vector"}: The time series to
      be tested.}
    \item{\code{model}:}{Object of class \code{"character"}: The model
      to be used, \emph{i.e.} intercept, trend or both}
    \item{\code{lag}:}{Object of class \code{"integer"}: The highest
      number of lags to include in the test regression.}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: The t-statistic.}
    \item{\code{cval}:}{Object of class \code{"vector"}: Critical values
      at the 1\%, 5\% and 10\% level of significance.}
    \item{\code{bpoint}:}{Object of class \code{"integer"}: The
      potential break point.}
    \item{\code{tstats}:}{Object of class \code{"vector"} The
      t-statistics of the rolling regression.}
    \item{\code{res}:}{Object of class \code{"vector"} The residuals of
      the test regression.}
    \item{\code{test.name}:}{Object of class \code{"character"} The name
      of the test, \emph{i.e.} `Zivot \& Andrews'.}
    \item{\code{testreg}:}{Object of class \code{"ANY"} The summary
      output of the test regression.}
  }
}
\section{Extends}{
Class \code{`urca'}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ur.za")} at the R prompt for a
  complete list of methods which are available for this class.
  
  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic and critical values.}
    \item{\code{summary}:}{like show, but summary of test regression
      added.}
    \item{\code{plot}:}{plot of recursive t-statistics.}
  }
}
\references{
  Zivot, E. and Andrews, Donald W.K. (1992), Further Evidence on the
  Great Crash, the Oil-Price Shock, and the Unit-Root Hypothesis,
  \emph{Journal of Business \& Economic Statistics}, \bold{10(3)},
  251--270.

  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{
  \code{\link{ur.za}} and \code{\link{urca-class}}. 
}
\author{Bernhard Pfaff}
\keyword{classes}
