\name{cajolst}
\alias{cajolst}
\title{Testing Cointegrating Rank with Level Shift at Unknown time}
\description{
  The function \code{cajolst} implements the procedure by Luetkepohl
  \emph{et al.} to test for the cointegration rank of a VAR process with
  a level shift at an unknown time.
}
\usage{
cajolst(x, type = c("eigen", "trace"), constant = FALSE,
        spec = c("longrun", "transitory"), K = 2, season = NULL)
}
\arguments{
  \item{x}{Data matrix to be investigated for cointegration.}
  \item{type}{The test to be conducted, either \code{`eigen'} or
    \code{`trace'}.}
  \item{constant}{If TRUE, no linear trend is specified.}
  \item{K}{The lag order of the series (levels) in the VAR.}
  \item{spec}{Determines the specification of the VECM, see details
    \code{ca.jo()}.}
  \item{season}{If seasonal dummies should be included, the data
    frequency must be set accordingly, \emph{i.e} `4' for quarterly data.}
}
\details{
  Note, that the slot \code{"x"} of the returned object contains the
  adjusted data series, that is, a matrix adjusted for the temptative
  break point, and if applicable, a linear trend and/or seasonal
  effects. The VECM is then estimated and tested for cointegration rank
  subject to the adjusted matrix. The break point is contained in the
  slot \code{"bp"}. 
}
\value{
  Returns an object of class \code{`ca.jo'}.
}
\references{
  
  L\"utkepohl, H., Saikkonen, P. and Trenkler, C. (2004), Testing for the
  Cointegrating Rank of a VAR Process with Level Shift at Unknown Time,
  \emph{Econometrica}, \bold{Vol. 72, No. 2}, 647--662.
  
}
\seealso{
  
  \code{\link{plotres}}, \code{\link{alrtest}}, \code{\link{ablrtest}},
  \code{\link{blrtest}}, \code{\link{ca.jo}}, \code{\link{cajools}},
  \code{\link{lttest}}, \code{\link{ca.jo-class}} and \code{\link{urca-class}}. 
}
\examples{
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.lst <- cajolst(sjd, constant=TRUE, type="eigen", K=2, season=4)
summary(sjd.lst)
}
\keyword{regression}
