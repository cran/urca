\name{ecb}
\alias{ecb}
\docType{data}
\encoding{latin1}
\title{Macroeconomic data of the Euro Zone}
\description{
  This data set contains some macroeconomic figures of the Euro Zone in
  order to estimate an exemplary money demand function.
}
\usage{data(ecb)}
\format{
  A data frame containing five series.
  
  \tabular{rl}{
    \code{period} \tab Time index from Q31997 until Q42003.\cr  
    \code{gdp.defl} \tab Gross Domestic Product Deflator,\cr 
    \tab [Index 2000=100, seasonally adjusted] \cr
    \code{gdp.nom} \tab Nominal Gross Domestic Product, \cr 
    \tab [Current prices, EUR billions, seasonally adjusted] \cr 
    \code{m3} \tab Monetary Aggregate M3, \cr
    \tab [outstanding amount at end of quarter, EUR billions, seasonally
    adjusted] \cr 
    \code{rl} \tab Benchmark Government Bond yield with a maturity of 10
    years, \cr
    \tab [percentages per annum, average of last quarter's month].
    }
}
\source{
  European Central Bank, Monthly Bulletins, Frankfurt am Main, Germany.
}
\references{
 \url{http://www.ecb.int}
}
\author{Bernhard Pfaff}
\keyword{datasets}
\concept{data set ECB money demand macroeconomic data Eurozone}
