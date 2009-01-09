\name{ur.kpss}
\alias{ur.kpss}
\encoding{latin1}
\title{Kwiatkowski et al. Unit Root Test}
\description{
  Performs the KPSS unit root test, where the Null hypothesis is
  stationarity. The test types specify as deterministic component either
  a constant \code{"mu"} or a constant with linear trend \code{"tau"}.
}
\usage{
ur.kpss(y, type = c("mu", "tau"), lags = c("short", "long", "nil"),
        use.lag = NULL)
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{type}{Type of deterministic part.}
  \item{lags}{Maximum number of lags used for error term correction.}
  \item{use.lag}{User specified number of lags.}
}
\details{
  \code{lags="short"} sets the number of lags to \eqn{\root 4 \of {4
      \times (n/100)}},
  whereas \code{lags="long"} sets the number of lags to \eqn{\root 4 \of
    {12 \times (n/100)}}. If \code{lags="nil"} is choosen, then no error
  correction is made. Furthermore, one can specify a different number of
  maximum lags by setting \code{use.lag} accordingly.
}
\value{
  An object of class \code{ur.kpss}.
}
\references{
  Kwiatkowski, D., Phillips, P.C.B., Schmidt, P. and Shin, Y., (1992),
  Testing the Null Hypothesis of Stationarity Against the Alternative of
  a Unit Root: How Sure Are We That Economic Time Series Have a Unit
  Root?, \emph{Journal of Econometrics}, \bold{54}, 159--178.

  Download possible at: \url{http://cowles.econ.yale.edu/}, see rubric
  'Discussion Papers (CFDPs)'.
}
\seealso{\code{\link{ur.kpss-class}}}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
gnp.l <- log(gnp)
kpss.gnp <- ur.kpss(gnp.l, type="tau", lags="short")
summary(kpss.gnp)
}
\author{Bernhard Pfaff}
\keyword{regression}
