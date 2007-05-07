\name{ur.sp-class}
\docType{class}
\alias{ur.sp-class}

\title{Representation of class `ur.sp'}
\description{
  This class contains the relevant information by applying the Schmidt
  \& Phillips unit root test to a time series.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"vector"}: The time series to
      be tested.}
    \item{\code{type}:}{Object of class \code{"character"}: Test type,
      \code{"rho"} or \code{"tau"} test statistic.}
    \item{\code{polynomial}:}{Object of class \code{"integer"}:
      Deterministic trend specification}
    \item{\code{signif}:}{Object of class \code{"numeric"}: Critical values.}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: Value of
      the test statistic.}
    \item{\code{cval}:}{Object of class \code{"numeric"}: The critical
      values, depending on \code{"signif"}, \code{"polynomial"} and the
      sample size.}
    \item{\code{res}:}{Object of class \code{"vector"}: The residuals of
      the test regression.}
    \item{\code{testreg}:}{Object of class \code{"ANY"}: The summary
      output of the test regression.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e.} `"Schmidt \& Phillips'.}
  }
}
\section{Extends}{
Class \code{`"urca'}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ur.sp")} at the R prompt for a
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
  Schmidt, P. and Phillips, P.C.B. (1992), LM Test for a Unit Root in
  the Presence of Deterministic Trends, \emph{Oxford Bulletin of Economics and
  Statistics}, \bold{54(3)}, 257--287.

  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{
  \code{\link{ur.sp}} and \code{\link{urca-class}}. 
}
\author{Bernhard Pfaff}
\keyword{classes}
