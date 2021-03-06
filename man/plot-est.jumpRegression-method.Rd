% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plot.R
\docType{methods}
\name{plot,est.jumpRegression-method}
\alias{plot,est.jumpRegression-method}
\title{Plot method for the Bayesian estimation results}
\usage{
\S4method{plot}{est.jumpRegression}(x, par.options, style = c("chains", "acf",
  "density"), par2plot, reduced = FALSE, thinning, burnIn,
  priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...)
}
\arguments{
\item{x}{est.jumpRegression class, created with method \code{\link{estimate,jumpRegression-method}}}

\item{par.options}{list of options for function par()}

\item{style}{one out of "chains", "acf", "density"}

\item{par2plot}{logical vector, which parameters to be plotted, order: \eqn{(\phi, \theta, \gamma^2, \xi, N)}}

\item{reduced}{logical (1), if TRUE, the chains are thinned and burn-in phase is dropped}

\item{thinning}{thinning rate, if missing, the proposed one by the estimation procedure is taken}

\item{burnIn}{burn-in phase, if missing, the proposed one by the estimation procedure is taken}

\item{priorMeans}{logical(1), if TRUE (default), prior means are marked with a line}

\item{col.priorMean}{color of the prior mean line, default 2}

\item{lty.priorMean}{linetype of the prior mean line, default 1}

\item{...}{optional plot parameters}
}
\description{
Plot method for the estimation results of the jump regression model.
}
\examples{
model <- set.to.class("jumpRegression", fun = function(t, N, theta) exp(theta[1]*t) + theta[2]*N,
  parameter = list(theta = c(2, 2), gamma2 = 0.25, xi = c(3, 0.5)),
  Lambda = function(t, xi) (t/xi[2])^xi[1])
data <- simulate(model, t = seq(0, 1, by = 0.01), plot.series = TRUE)
est <- estimate(model, t = seq(0, 1, by = 0.01), data, 1000)  # nMCMC small for example
plot(est)
plot(est, burnIn = 100, thinning = 2, reduced = TRUE)
plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2, 3)), xlab = "iteration")
plot(est, style = "acf", main = "", par2plot = c(TRUE, FALSE, FALSE, TRUE, TRUE))
plot(est, style = "density", lwd = 2, priorMean = FALSE)
plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
plot(est, style = "acf", par.options = list(), par2plot = c(TRUE, rep(FALSE, 4)), main = "")
}

