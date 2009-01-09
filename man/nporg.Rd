\name{nporg}
\docType{data}
\alias{nporg}
\encoding{latin1}
\title{Nelson \& Plosser original data set}
\description{
  This data set contains the fourteen U.S. economic time series used by
  Nelson \& Plosser in their seminal paper.
}
\usage{data(nporg)}
\format{
  A data frame containing fourteen series.

  \tabular{rl}{
    \code{year} \tab Time index from 1860 until 1970. \cr 
    \code{gnp.r} \tab Real GNP, \cr
    \tab [Billions of 1958 Dollars], [1909 -- 1970] \cr
    \code{gnp.n} \tab Nominal GNP, \cr
    \tab [Millions of Current Dollars], [1909 -- 1970] \cr
    \code{gnp.pc} \tab Real Per Capita GNP, \cr
    \tab [1958 Dollars], [1909 -- 1970] \cr
    \code{ip} \tab Industrial Production Index, \cr
    \tab [1967 = 100], [1860 -- 1970] \cr
    \code{emp} \tab Total Employment, \cr
    \tab [Thousands], [1890 -- 1970] \cr
    \code{ur} \tab Total Unemployment Rate, \cr
    \tab [Percent], [1890 -- 1970] \cr
    \code{gnp.p} \tab GNP Deflator, \cr
    \tab [1958 = 100], [1889 -- 1970] \cr
    \code{cpi} \tab Consumer Price Index, \cr
    \tab [1967 = 100], [1860 -- 1970] \cr
    \code{wg.n} \tab Nominal Wages \cr
    \tab (Average annual earnings per full-time employee in manufacturing), \cr
    \tab [current Dollars], [1900 -- 1970] \cr
    \code{wg.r} \tab Real Wages, \cr
    \tab [Nominal wages/CPI], [1900 -- 1970] \cr
    \code{M} \tab Money Stock (M2), \cr
    \tab [Billions of Dollars, annual averages], [1889 -- 1970] \cr
    \code{vel} \tab Velocity of Money, \cr
    \tab [1869 -- 1970] \cr
    \code{bnd} \tab Bond Yield (Basic Yields of 30-year
        corporate bonds), \cr
    \tab [Percent per annum], [1900 -- 1970] \cr
    \code{sp} \tab Stock Prices, \cr
    \tab [Index; 1941 -- 43 = 100], [1871 -- 1970] \cr
    }
}
\source{
  Nelson, C.R. and Plosser, C.I. (1982), Trends and Random Walks in
  Macroeconomic Time Series, \emph{Journal of Monetary Economics},
  \bold{10}, 139--162.
}
\references{
  \url{http://korora.econ.yale.edu/phillips/index.htm}
}
\author{Bernhard Pfaff}
\keyword{datasets}
\concept{data set Nelson Plosser macroeconomic variables}
