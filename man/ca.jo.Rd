\name{ca.jo}
\alias{ca.jo}
\title{Johansen Procedure for VAR}
\description{
  Conducts the Johansen procedure on a given data set. The
  \code{"trace"} or \code{"eigen"} statistics are reported and the
  matrix of eigenvectors as well as the loading matrix.
}
\usage{
ca.jo(x, type = c("eigen", "trace"), constant = FALSE, K = 2, season = NULL)
}
\arguments{
  \item{x}{Data matrix to be investigated for cointegration.}
  \item{type}{The test to be conducted, either \code{`eigen'} or
    \code{`trace'}.}
  \item{constant}{If TRUE, no linear trend is specified.}
  \item{K}{The lag order of the series (levels) in the VAR.}
  \item{season}{If seasonal dummies should be included, the data
    frequency must be set accordingly, \emph{i.e} `4' for quarterly data.}
}
\details{

  The following VAR is estimated by the function \code{ca.jo()}:

  \deqn{\Delta \bold{X}_t = \bold{\Gamma}_1 \Delta \bold{X}_{t-1} +
    \dots + \bold{\Gamma}_{k-1} \Delta \bold{X}_{t-k+1} + \bold{\Pi X}_{t-k} + \bold{\mu} +
    \bold{\Phi D}_t + \bold{\varepsilon}_t}

  where

  \deqn{\bold{\Gamma}_i = - (\bold{I} - \bold{\Pi}_1 - \dots -
    \bold{\Pi}_i), \quad (i = 1, \dots , k-1),}

  and

  \deqn{\bold{\Pi} = -(\bold{I} - \bold{\Pi}_1 - \dots - \bold{\Pi}_k)}

  in a general VAR of the form:

  \deqn{\bold{X}_t = \bold{\Pi}_1 \bold{X}_{t-1} + \dots + \bold{\Pi}_k
    \bold{X}_{t-k} + \bold{\mu} + \bold{\Phi D}_t + \bold{\varepsilon}_t
    , \quad (t = 1, \dots, T).}


  If \code{"season"} is not NULL, centered seasonal dummy variables are
  included. Critical values are only reported for systems with less than
  five variables.
 
}
\value{
  An object of class \code{`ca.jo'}.
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
  \code{\link{plotres}}, \code{\link{blrtest}},
  \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, constant=TRUE, type="eigen", K=2, season=4)
summary(sjd.vecm)
#
data(finland)
sjf <- finland
sjf.vecm <- ca.jo(sjf, constant=FALSE, type="eigen", K=2, season=4)
summary(sjf.vecm)
}
\keyword{regression}
