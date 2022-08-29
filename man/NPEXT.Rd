\name{npext}
\docType{data}
\alias{npext}
\encoding{latin1}
\title{Nelson and Plosser extended data set}
\description{
  This data set contains the fourteen U.S. economic time series used by
  Schotman and Dijk. All series are transformed by taking logarithms
  except for the bond yield. The sample period ends in 1988.
}
\usage{data(npext)}
\format{
  A data frame containing fourteen series.
		
  \tabular{rl}{
    \code{year} \tab Time index from 1860 until 1988. \cr 
    \code{realgnp} \tab Real GNP, [Billions of 1958 Dollars], \cr
    \tab [1909 -- 1988] \cr
    \code{nomgnp} \tab Nominal GNP, \cr
    \tab [Millions of Current Dollars], [1909 -- 1988] \cr
    \code{gnpperca} \tab Real Per Capita GNP, \cr
    \tab [1958 Dollars], [1909 -- 1988] \cr
    \code{indprod} \tab Industrial Production Index, \cr
    \tab [1967 = 100], [1860 -- 1988] \cr
    \code{employmt} \tab Total Employment, \cr
    \tab [Thousands], [1890 -- 1988] \cr
    \code{unemploy} \tab Total Unemployment Rate, \cr
    \tab [Percent], [1890 -- 1988] \cr
    \code{gnpdefl} \tab GNP Deflator, \cr
    \tab [1958 = 100], [1889 -- 1988] \cr
    \code{cpi} \tab Consumer Price Index, \cr
    \tab [1967 = 100], [1860 -- 1988] \cr
    \code{wages} \tab Nominal Wages \cr
    \tab (Average annual earnings per full-time employee in manufacturing), \cr
    \tab [current Dollars], [1900 -- 1988] \cr
    \code{realwag} \tab Real Wages, \cr
    \tab [Nominal wages/CPI], [1900 -- 1988] \cr
    \code{M} \tab Money Stock (M2), \cr
    \tab [Billions of Dollars, annual averages], [1889 -- 1988] \cr
    \code{velocity} \tab Velocity of Money, \cr
    \tab [1869 -- 1988] \cr
    \code{interest} \tab Bond Yield (Basic Yields of 30-year
       corporate bonds), \cr
    \tab [Percent per annum], [1900 -- 1988] \cr
    \code{sp500} \tab Stock Prices, \cr
    \tab [Index; 1941 -- 43 = 100], [1871 -- 1988] \cr
    }
}
\source{
  Schotman, P.C. and van Dijk, H.K. (1991), On Bayesian Routes to Unit
  Roots, \emph{Journal of Applied Econometrics}, \bold{6}, 387--401.

  Koop, G. and Steel, M.F.J. (1994), A Decision-Theoretic Analysis of
  the Unit-Root Hypothesis using Mixtures of Elliptical Models,
  \emph{Journal of Business and Economic Statistics}, \bold{12}, 95--107. 
}
\references{\url{https://www.amstat.org}}
\author{Bernhard Pfaff}
\keyword{datasets}
\concept{data set Nelson Plosser macroeconomic variables}
