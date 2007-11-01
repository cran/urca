\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,ur.ers,missing-method}
\alias{plot,ur.kpss,missing-method}
\alias{plot,ca.jo,missing-method}
\alias{plot,ca.po,missing-method}
\alias{plot,ur.pp,missing-method}
\alias{plot,ur.sp,missing-method}
\alias{plot,ur.za,missing-method}
\alias{plot,ur.df,missing-method}
\title{Methods for Function plot in Package `urca'}
\description{
  Plot methods for objects belonging to classes set in package
  \code{`urca'}. Depending on the unit root/cointegration test a
  suitable graphical presentation is selected.
  }
\section{Methods}{\describe{

\item{x = "ur.ers", y = "missing"}{Diagram of fit of the Elliott,
  Rothenberg \& Stock unit root test of type \code{"DF-GLS"} with
  residual plot and their acfs' and pacfs'.}
\item{x = "ur.kpss", y = "missing"}{Residual plot and their acfs' and
  pacfs' of the KPSS test.}
\item{x = "ca.jo", y = "missing"}{Time series plots and associated
  cointegration relations for the Johansen procedure.}
\item{x = "ca.po", y = "missing"}{Residual plot and their acfs' and
  pacfs' of the cointegration regression(s) for the Phillips \& Ouliaris
  test.}
\item{x = "ur.pp", y = "missing"}{Diagram of fit of the Phillips \&
  Perron unit root test, residual plot and their acfs' and pacfs'.}
\item{x = "ur.sp", y = "missing"}{Diagram of fit of the Schmidt \&
  Phillips unit root test, residual plot and their acfs' and pacfs'.}
\item{x = "ur.za", y = "missing"}{Plot of recursive t-statistics
  as outcome of Zivot \& Andrews unit root test.}
}}
\seealso{
  \code{\link{ur.ers-class}}, \code{\link{ur.kpss-class}},
  \code{\link{ca.jo-class}}, \code{\link{ca.po-class}},
  \code{\link{ur.pp-class}}, \code{\link{ur.sp-class}} and
  \code{\link{ur.za-class}}.}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
gnp.l <- log(gnp)
#
ers.gnp <- ur.ers(gnp, type="DF-GLS", model="trend", lag.max=4)
plot(ers.gnp)
#
kpss.gnp <- ur.kpss(gnp.l, type="tau", lags="short")
plot(kpss.gnp)
#
pp.gnp <- ur.pp(gnp, type="Z-tau", model="trend", lags="short")
plot(pp.gnp)
#
sp.gnp <- ur.sp(gnp, type="tau", pol.deg=1, signif=0.01)
plot(sp.gnp)
#
za.gnp <- ur.za(gnp, model="both", lag=2)
plot(za.gnp)
#
data(denmark)
sjd <- denmark[, c("LRM", "LRY", "IBO", "IDE")]
sjd.vecm <- ca.jo(sjd, ecdet="const", type="eigen", K=2, season=4)
plot(sjd.vecm)
}
\author{Bernhard Pfaff}
\keyword{methods}

