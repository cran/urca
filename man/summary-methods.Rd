\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,ur.ers-method}
\alias{summary,ur.kpss-method}
\alias{summary,ca.jo-method}
\alias{summary,cajo.test-method}
\alias{summary,ca.po-method}
\alias{summary,ur.pp-method}
\alias{summary,ur.df-method}
\alias{summary,ur.sp-method}
\alias{summary,ur.za-method}
\encoding{latin1}
\title{Methods for Function summary in Package `urca'}
\description{
  Summarises the outcome of unit root/cointegration tests by creating a new object of class \code{sumurca}.
}
\section{Methods}{\describe{
    \item{object = "ur.df"}{The test type, its statistic, the test
      regression and the critical values for the Augmented Dickey and
      Fuller test are returned.}  
    \item{object = "ur.ers"}{The test type, its statistic and the
      critical values for the Elliott, Rothenberg \& Stock test are
      returned. In case of test \code{"DF-GLS"} the summary output
      of the test regression is provided, too.}
    \item{object = "ur.kpss"}{The test statistic, the critical value as
      well as the test type and the number of lags used for error
      correction for the Kwiatkowski \emph{et al.} unit root test is
      returned.}
    \item{object = "ca.jo"}{The \code{"trace"} or \code{"eigen"} statistic,
      the critical values as well as the eigenvalues, eigenvectors and
      the loading matrix of the Johansen procedure are reported.}
    \item{object = "cajo.test"}{The test statistic of a restricted VAR
      with respect to \eqn{\bold{\alpha}} and/or \eqn{\bold{\beta}} with
      p-value and degrees of freedom is reported. Furthermore, the
      restriction matrix(ces), the eigenvalues and eigenvectors as well
      as the loading matrix are returned.}  
    \item{object = "ca.po"}{The \code{"Pz"} or \code{"Pu"} statistic,
      the critical values as well as the summary output of the test
      regression for the Phillips \& Ouliaris cointegration test.}
    \item{object = "ur.pp"}{The Z statistic, the critical values as
      well as the summary output of the test regression for the Phillips
      \& Perron test, as well as the test statistics for the
      coefficients of the deterministic part is returned.}
    \item{object = "ur.df"}{The relevant tau statistic, the critical
      values as well as the summary output of the test regression for
      the augmented Dickey-Fuller test is returned.}
    \item{object = "ur.sp"}{The test statistic, the critical value as
      well as the summary output of the test regression for the Schmidt
      \& Phillips test is returned.}
    \item{object = "ur.za"}{The test statistic, the critical values as
      well as the summary output of the test regression for the Zivot \&
      Andrews test is returned.}
}}
\seealso{
  \code{\link{ur.ers-class}}, \code{\link{ur.kpss-class}},
  \code{\link{ca.jo-class}}, \code{\link{cajo.test-class}},
  \code{\link{ca.po-class}}, \code{\link{ur.pp-class}},
  \code{\link{ur.df-class}}, \code{\link{ur.sp-class}},
  \code{\link{ur.za-class}} and \code{\link{sumurca-class}}.}
\examples{
data(nporg)
gnp <- na.omit(nporg[, "gnp.r"])
gnp.l <- log(gnp)
#
ers.gnp <- ur.ers(gnp, type="DF-GLS", model="trend", lag.max=4)
summary(ers.gnp)
#
kpss.gnp <- ur.kpss(gnp.l, type="tau", lags="short")
summary(kpss.gnp)
#
pp.gnp <- ur.pp(gnp, type="Z-tau", model="trend", lags="short")
summary(pp.gnp)
#
df.gnp <- ur.df(gnp, type="trend", lags=4)
summary(df.gnp)
#
sp.gnp <- ur.sp(gnp, type="tau", pol.deg=1, signif=0.01)
summary(sp.gnp)
#
za.gnp <- ur.za(gnp, model="both", lag=2)
summary(za.gnp)
#
data(finland)
sjf <- finland
sjf.vecm <- ca.jo(sjf, ecdet="none", type="eigen", K=2, season=4)
summary(sjf.vecm)
#
HF0 <- matrix(c(-1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), c(4,3))
summary(blrtest(sjf.vecm, H=HF0, r=3))
}
\author{Bernhard Pfaff}
\keyword{methods}
