\name{ur.df-class}
\docType{class}
\alias{ur.df-class}
\title{Representation of class `ur.df'}
\description{
  This class contains the relevant information by applying the augmented
  Dickey-Fuller unit root test to a time series.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"vector"}: The time series to
      be tested.}
    \item{\code{model}:}{Object of class \code{"character"}: The type of
      the deterministic part, either \code{"none"}, \code{"drift"} or
      \code{"trend"}. The latter includes a constant term, too.}
    \item{\code{lags}:}{Object of class \code{"integer"}: Number of lags
      for error correction.}
    \item{\code{cval}:}{Object of class \code{"matrix"}: Critical values
      at the 1\%, 5\% and 10\% level of significance.}
    \item{\code{teststat}:}{Object of class \code{"matrix"}: Value of
      the test statistic.}
    \item{\code{testreg}:}{Object of class \code{"ANY"}: The summary
      output of the test regression.}
    \item{\code{res}:}{Object of class \code{"vector"}: The residuals of
      the test regression.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e} `Augmented-Dickey-Fuller Test'.}
  }
}
\section{Extends}{
Class \code{`urca'}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ur.df")} at the R prompt for a
  complete list of methods which are available for this class.
  
  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic.}
    \item{\code{summary}:}{like show, but critical value and summary of
      test regression added.}
    \item{\code{plot}:}{Residual plot, acfs' and pacfs'.}
  }
}
\references{
  Dickey, D. A. and Fuller, W. A. (1979), Distributions of the
  Estimators For Autoregressive Time Series with a Unit Root,
  \emph{Journal of the American Statistical Association},  \bold{75}, 427--431.
  
  Dickey, D. A. and Fuller, W. A. (1981), Likelihood Ratio Statistics
  for Autoregressive Time Series with a Unit Root, \emph{Econometrica},
  \bold{49}, 1057--1072.

  Hamilton (1994), \emph{Time Series Analysis}, Princeton University Press.
}
\seealso{
  \code{\link{ur.df}} and \code{\link{urca-class}}
}
\author{Bernhard Pfaff}
\keyword{classes}
