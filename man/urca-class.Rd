\name{urca-class}
\docType{class}
\alias{urca-class}

\title{Class `urca'. Parent of classes in package `urca'}
\description{
  This class is the parent class of the specific classes
  designed holding the test specific information of the unit
  root/cointegration tests.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("urca", ...)},
but most often the slot \code{test.name} is set by calling one of the
unit root/cointegration functions, \emph{e.g} \code{ur.za}.
}
\section{Slots}{
  \describe{
    \item{\code{test.name}:}{Object of class \code{"character"}. The
      name of the unit root/cointegration test.}
  }
}
\section{Methods}{
No methods defined with class `urca'.
}
\seealso{
  \code{\link{ur.ers-class}}, \code{\link{ur.kpss-class}},
  \code{\link{ca.jo-class}}, \code{\link{ca.po-class}},
  \code{\link{ur.pp-class}}, \code{\link{ur.sp-class}}
  and \code{\link{ur.za-class}}.
}
\author{Bernhard Pfaff}
\keyword{classes}
