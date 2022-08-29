\name{ur.pp-class}
\docType{class}
\alias{ur.pp-class}
\encoding{latin1}
\title{Representation of class ur.pp}
\description{
  This class contains the relevant information by applying the Phillips
  and Perron unit root test to a time series.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"vector"}: The time series to
      be tested.}
    \item{\code{type}:}{Object of class \code{"character"}: Test type of
      Z statistic, either \code{"Z-alpha"} or \code{"Z-tau"}.}
    \item{\code{model}:}{Object of class \code{"character"}: The type of
      the deterministic part, either \code{"constant"} or
      \code{"trend"}. The latter includes a constant term, too.}
    \item{\code{lag}:}{Object of class \code{"integer"}: Number of lags
      for error correction.}
    \item{\code{cval}:}{Object of class \code{"matrix"}: Critical values
      at the 1\%, 5\% and 10\% level of significance.}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: Value of
      the test statistic.}
    \item{\code{testreg}:}{Object of class \code{"ANY"}: The summary
      output of the test regression.}
    \item{\code{auxstat}:}{Object of class \code{"matrix"}: Test
      statistic(s) of the deterministic part.}
    \item{\code{res}:}{Object of class \code{"vector"}: The residuals of
      the test regression.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e} `Phillips-Perron'.}
  }
}
\section{Extends}{
Class \code{urca}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ur.pp")} at the R prompt for a
  complete list of methods which are available for this class.
  
  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic.}
    \item{\code{summary}:}{like show, but critical value and summary of
      test regression added.}
    \item{\code{plot}:}{Diagram of fit plot, residual plot and their
      acfs' and pacfs'.}
  }
}
\references{
  Phillips, P.C.B. and Perron, P. (1988), Testing for a unit root in
  time series regression, \emph{Biometrika},  \bold{75(2)}, 335--346.

  MacKinnon, J.G. (1991), Critical Values for Cointegration Tests,
  \emph{Long-Run Economic Relationships}, eds. R.F. Engle and
  C.W.J. Granger, London, Oxford, 267--276.
  
  Download possible at: \url{https://cowles.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{
  \code{\link{ur.pp}} and \code{\link{urca-class}}
}
\author{Bernhard Pfaff}
\keyword{classes}
