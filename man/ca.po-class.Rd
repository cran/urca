\name{ca.po-class}
\docType{class}
\alias{ca.po-class}

\title{Representation of class `ca.po'}
\description{
  This class contains the relevant information by applying the Phillips
  \& Ouliaris cointegration test to a data matrix.
}
\section{Slots}{
  \describe{
    \item{\code{z}:}{Object of class \code{"ANY"}: A data matrix, or an
      object that can be coerced to it.}
    \item{\code{type}:}{Object of class \code{"character"}: The type of
      the test, either the \code{"Pu"}-test or the normalisation
	invariant \code{"Pz"}-test.}
    \item{\code{model}:}{Object of class \code{"character"}: Determines
      how the series should be detrended.}
    \item{\code{lag}:}{Object of class \code{"integer"}: The lags used
      for variance/covariance correction.}
    \item{\code{cval}:}{Object of class \code{"matrix"}: The critical
      values of the test at the 1\%, 5\% and 10\% level of significance.}
    \item{\code{res}:}{Object of class \code{"matrix"}: The residuals of
      the the cointegration regression(s).}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: The value
      of the test statistic.}
    \item{\code{testreg}:}{Object of class \code{"ANY"}: The summary
      output of the cointegration regression(s).}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e.} `Phillips \& Ouliaris'.}
  }
}
\section{Extends}{
Class \code{`urca'}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ca.po")} at the R prompt for a
  complete list of methods which are available for this class.
  
  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic.}
    \item{\code{summary}:}{like show, but critical value and summary of
      test regression(s) added.}
    \item{\code{plot}:}{Residual plot(s) and their acfs' and pacfs'.}
  }
}
\references{
  Phillips, P.C.B. and Ouliaris, S. (1990), Asymptotic Properties of
  Residual Based Tests for Cointegration, \emph{Econometrica},
  \bold{Vol. 58, No. 1}, 165--193.
}
  
\seealso{
  \code{\link{ca.po}} and \code{\link{urca-class}}. 
}
\author{Bernhard Pfaff}
\keyword{classes}
