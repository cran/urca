\name{show.urca}
\alias{show.urca}
\encoding{latin1}
\title{Function to show objects of classes for unit root tests}
\description{
  The function \code{show.urca} is called within the defined methods
  for classes \code{ca.jo}, \code{cajo.test}, \code{ca.po},
  \code{ur.df}, \code{ur.ers}, \code{ur.kpss}, \code{ur.po}, \code{ur.pp},
  \code{ur.df}, \code{ur.sp} and \code{ur.za}.
 
}
\usage{
show.urca(object)
}
\arguments{
  \item{object}{Object of class contained in \code{urca}.}
}
\details{
  This function is called by method \code{show}.
}
\value{
  The Name and test statistic of a unit root/cointegration test.
}
\author{Bernhard Pfaff}
\keyword{methods}
