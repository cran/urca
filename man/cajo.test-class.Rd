\name{cajo.test-class}
\docType{class}
\alias{cajo.test-class}
\title{Representation of class `cajo.test'}
\description{
  This class contains the relevant information by estimating and testing
  a VAR under linear restrictions on \eqn{\bold{\alpha}} and
  \eqn{\bold{\beta}}.
}
\section{Slots}{
  \describe{
    \item{\code{Z0}:}{Object of class \code{"matrix"}: The matrix of the
      differenced series.}
    \item{\code{Z1}:}{Object of class \code{"matrix"}: The regressor
      matrix, except for the lagged variables in levels.}
    \item{\code{ZK}:}{Object of class \code{"matrix"}: The matrix of the
      lagged variables in levels.}
    \item{\code{const}:}{Object of class \code{"logic"}: If TRUE, no
      linear trend is included and \emph{vice versa}.}
    \item{\code{H}:}{Object of class \code{"ANY"}: The matrix
      containing the restrictions placed upon \eqn{\bold{\beta}}.}
    \item{\code{A}:}{Object of class \code{"ANY"}: The matrix
      containing the restrictions placed upon \eqn{\bold{\alpha}}.}
    \item{\code{B}:}{Object of class \code{"ANY"}: The matrix
      orthogonal to matrix \eqn{\bold{A}}.}
    \item{\code{type}:}{Object of class \code{"character"}: The test type.}
    \item{\code{teststat}:}{Object of class \code{"numeric"}: The value
      of the test statistic.}
    \item{\code{pval}:}{Object of class \code{"vector"}: The p-value and
      the degrees of freedom.}
    \item{\code{lambda}:}{Object of class \code{"vector"}: The
      eigenvalues of the restricted model.}
    \item{\code{Vorg}:}{Object of class \code{"matrix"}: The matrix of
      eigenvectors, such that \eqn{\hat V_{\dots}'(H'S_{\dots}H)\hat
	V_{\dots} = I}.}
    \item{\code{V}:}{Object of class \code{"matrix"}: The matrix of the
      restricted eigenvectors, normalised with respect to the first variable.}
    \item{\code{W}:}{Object of class \code{"matrix"}: The matrix of the
      corresponding loading weights.}
    \item{\code{PI}:}{Object of class \code{"matrix"}: The coefficient
      matrix of the lagged variables in levels.}
    \item{\code{DELTA}:}{Object of class \code{"ANY"}: The
      variance/covarinace matrix of \eqn{\bold{V}}.}
    \item{\code{DELTA.bb}:}{Object of class \code{"ANY"}: The
      variance/covarinace matrix of the marginal factor
      \eqn{\bold{B}'\bold{R}_{0t}}.}
    \item{\code{DELTA.ab}:}{Object of class \code{"ANY"}: The
      variance/covarinace matrix of the conditional distribution of
      \eqn{\bold{A}'\bold{R}_{0t}} and \eqn{\bold{R}_{kt}}.}
    \item{\code{DELTA.aa.b}:}{Object of class \code{"ANY"}: The
      variance/covarinace matrix of the restricted loading matrix.}
    \item{\code{GAMMA}:}{Object of class \code{"matrix"}: The
      coefficient matrix of \eqn{\bold{Z1}}.}
    \item{\code{test.name}:}{Object of class \code{"character"}: The
      name of the test, \emph{i.e.} `Johansen-Procedure'.}
  }
}
\section{Extends}{
Class \code{`urca'}, directly.
}
\section{Methods}{
  Type \code{showMethods(classes="cajo.test")} at the R prompt for a
  complete list of methods which are available for this class.

  Useful methods include
  \describe{
    \item{\code{show}:}{test-statistic.}
    \item{\code{summary}:}{like show, but p-value of test statistic,
      restricted eigenvectors, loading matrix and restriction matrices
      \eqn{\bold{H}} and \eqn{\bold{A}}, where applicable, added.}
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
  \code{\link{ablrtest}}, \code{\link{alrtest}}, \code{\link{blrtest}},
  \code{\link{ca.jo}}, \code{\link{ca.jo-class}} and \code{\link{urca-class}}.
}
\author{Bernhard Pfaff}
\keyword{classes}
