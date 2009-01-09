\name{ca.jo-class}
\docType{class}
\alias{ca.jo-class}
\encoding{latin1}
\title{Representation of class ca.jo}
\description{
  This class contains the relevant information by applying the Johansen
  procedure to a matrix of time series data.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"ANY"}: A data matrix, or an
      object that can be coerced to it.}
    \item{\code{Z0}:}{Object of class \code{"matrix"}: The matrix of the
      differenced series.}
    \item{\code{Z1}:}{Object of class \code{"matrix"}: The regressor
      matrix, except for the lagged variables in levels.}
    \item{\code{ZK}:}{Object of class \code{"matrix"}: The matrix of the
      lagged variables in levels.}
    \item{\code{type}:}{Object of class \code{"character"}: The type of the
      test, either \code{"trace"} or \code{"eigen"}.}
    \item{\code{model}:}{Object of class \code{"character"}: The model
      description in prose, with respect to the inclusion of a linear
      trend.}
    \item{\code{ecdet}:}{Object of class \code{"character"}: Specifies
      the deterministic term to be included in the cointegration
      relation. This can be either "none", "const", or "trend".}
    \item{\code{lag}:}{Object of class \code{"integer"}: The lag order
      for the variables in levels.}
    \item{\code{P}:}{Object of class \code{"integer"}: The count of
      variables.}
    \item{\code{season}:}{Object of class \code{"ANY"}: The frequency of
      the data, if seasonal dummies should be included, otherwise NULL.}
    \item{\code{dumvar}:}{Object of class \code{"ANY"}: A matrix
      containing dummy variables. The row dimension must be equal to
      \code{x}, otherwise NULL.}
    \item{\code{cval}:}{Object of class \code{"ANY"}: The critical
      values of the test at the 1\%, 5\% and 10\% level of significance.}
    \item{\code{teststat}:}{Object of class \code{"ANY"}: The values
      of the test statistics.}
    \item{\code{lambda}:}{Object of class \code{"vector"}: The eigenvalues.}
    \item{\code{Vorg}:}{Object of class \code{"matrix"}: The matrix of
      eigenvectors, such that \eqn{\hat V'S_{kk}\hat V = I}.}
    \item{\code{V}:}{Object of class \code{"matrix"}: The matrix of
      eigenvectors, normalised with respect to the first variable.}
    \item{\code{W}:}{Object of class \code{"matrix"}: The matrix of
      loading weights.}
    \item{\code{PI}:}{Object of class \code{"matrix"}: The coeffcient
      matrix of the lagged variables in levels.}
    \item{\code{DELTA}:}{Object of class \code{"matrix"}: The
      variance/covarinace matrix of \code{V}.}
    \item{\code{GAMMA}:}{Object of class \code{"matrix"}: The
      coeffecient matrix of \code{Z1}.}
    \item{\code{R0}:}{Object of class \code{"matrix"}: The matrix of
      residuals from the regressions in differences.}
    \item{\code{RK}:}{Object of class \code{"matrix"}: The matrix of
      residuals from the regression in lagged levels.}
    \item{\code{bp}:}{Object of class \code{"ANY"}: Potential break
      point, only set if function \code{cajolst} is called, otherwise
      \code{NA}.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e.} `Johansen-Procedure'.}
    \item{\code{spec}:}{Object of class \code{"character"}: The
      specification of the VECM.}
    \item{\code{call}:}{Object of class \code{"call"}: The
      call of function \code{ca.jo}.}
  }
}
\section{Extends}{
Class \code{urca}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="ca.jo")} at the R prompt for a
  complete list of methods which are available for this class.

  Useful methods include
  \describe{
    \item{\code{show}:}{test statistic.}
    \item{\code{summary}:}{like show, but critical values, eigenvectors
      and loading matrix added.}
    \item{\code{plot}:}{The series of the VAR and their potential
      cointegration relations.}
  }
}
\references{

  Johansen, S. (1988), Statistical Analysis of Cointegration Vectors,
  \emph{Journal of Economic Dynamics and Control}, \bold{12}, 231--254.

  Johansen, S. and Juselius, K. (1990), Maximum Likelihood Estimation and
  Inference on Cointegration -- with Applications to the Demand for
  Money, \emph{Oxford Bulletin of Economics and Statistics}, \bold{52,
    2}, 169--210.

  Johansen, S. (1991), Estimation and Hypothesis Testing of
  Cointegration Vectors in Gaussian Vector Autoregressive Models,
  \emph{Econometrica}, \bold{Vol. 59, No. 6}, 1551--1580.

}
\seealso{
  \code{\link{ca.jo}}, \code{\link{plotres}} and \code{\link{urca-class}}.
}
\author{Bernhard Pfaff}
\keyword{classes}
