\name{sumurca-class}
\docType{class}
\alias{sumurca-class}

\title{Representation of class `sumurca'}
\description{
  A class for objects returned by applying method summary() to objects
  from classes: \code{ur.ers}, \code{ca.jo}, \code{cajo.test},
  \code{ur.kpss}, \code{ca.po}, \code{ur.pp}, \code{ur.df}, \code{ur.sp}
  or \code{ur.za}.
 
}
\section{Slots}{
  \describe{
    \item{\code{classname}:}{The class name of the original object to
      which method summary is applied.}
    \item{\code{test.name}:}{The name of the test, \emph{i.e.}
      `Johansen-Procedure'.}
    \item{\code{testreg}:}{The test regression where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{teststat}:}{The test statististic where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{cval}:}{The critical values of the test where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{bpoint}:}{Potential break point where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{signif}:}{Significance level of the test where
      applicable, otherwise set to \code{NULL}.}
    \item{\code{model}:}{Description of the underlying model where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{type}:}{The test type where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{auxstat}:}{The result of an auxiliary regression where
      applicable, otherwise set to \code{NULL}.}
    \item{\code{lag}:}{The number of lags included where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{H}:}{The matrix containing the restrictions placed upon
      \eqn{\bold{\beta}} where applicable, otherwise set to \code{NULL}.}
    \item{\code{A}:}{The matrix containing the restrictions placed upon
      \eqn{\bold{\alpha}} where applicable, otherwise set to \code{NULL}.}
    \item{\code{lambda}:}{The eigenvalues where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{pval}:}{The p-value and the degrees of freedom where
      applicable, otherwise set to \code{NULL}.}
    \item{\code{V}:}{The matrix of eigenvectors, normalised with respect
      to the first variable where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{W}:}{The matrix of loading weights where applicable,
      otherwise set to \code{NULL}.}
    \item{\code{P}:}{The count of variables where applicable,
      otherwise set to \code{NULL}.}
  }
}
\section{Methods}{
  For this class a \code{print} method is available, that nicely prints the
  summary results of objects belonging to either one of the following
  classes: \code{ur.ers}, \code{ca.jo}, \code{cajo.test}, 
  \code{ur.kpss}, \code{ca.po}, \code{ur.pp}, \code{ur.df}, \code{ur.sp} or
  \code{ur.za}.
}
\seealso{
  \code{\link{summary}}, \code{\link{ur.ers-class}},
  \code{\link{ur.kpss-class}}, \code{\link{ca.jo-class}},
  \code{\link{cajo.test-class}}, \code{\link{ca.po-class}},
  \code{\link{ur.pp-class}}, \code{\link{ur.df-class}},
  \code{\link{ur.sp-class}} and \code{\link{ur.za-class}}. 
}
\author{Bernhard Pfaff}
\keyword{classes}
