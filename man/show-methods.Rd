\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,ur.ers-method}
\alias{show,ur.kpss-method}
\alias{show,ca.po-method}
\alias{show,ur.pp-method}
\alias{show,ur.sp-method}
\alias{show,ur.za-method}
\title{Methods for Function show in Package `urca'}
\description{
  Displays the outcome of the unit root/cointegration tests.
}
\section{Methods}{\describe{

\item{object = "ca.po"}{Displays the test-statistic of the Phillips \&
  Ouliaris cointegration test.}
\item{object = "ur.ers"}{Displays the test-statistic of the Elliott,
  Rothenberg \& Stock unit root test.}
\item{object = "ur.kpss"}{Displays the test-statistic of the Kwiatkowski
  \emph{et al.} unit root test.}
\item{object = "ur.pp"}{Displays the test-statistic of the Phillips \&
  Perron unit root test.}
\item{object = "ur.sp"}{Displays the test-statistic of the Schmidt \&
  Phillips unit root test.}
\item{object = "ur.za"}{Displays the test-statistic of the Zivot \&
  Andrews unit root test.}
}}
\seealso{
  \code{\link{ca.po-class}}, \code{\link{ur.ers-class}},
  \code{\link{ur.kpss-class}}, \code{\link{ur.pp-class}},
  \code{\link{ur.sp-class}} and \link{ur.za-class}}.
}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
gnp.l <- log(gnp)
#
ers.gnp <- ur.ers(gnp, type="DF-GLS", model="trend", lag.max=4)
show(ers.gnp)
#
kpss.gnp <- ur.kpss(gnp.l, type="tau", lags="short")
show(kpss.gnp)
#
pp.gnp <- ur.pp(gnp, type="Z-tau", model="trend", lags="short")
show(pp.gnp)
#
sp.gnp <- ur.sp(gnp, type="tau", pol.deg=1, signif=0.01)
show(sp.gnp)
#
za.gnp <- ur.za(gnp, model="both", lag=2)
show(za.gnp)
}
\keyword{methods}
