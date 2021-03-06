% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/estimate.R
\docType{methods}
\name{estimate,jumpRegression-method}
\alias{estimate,jumpRegression-method}
\title{Estimation for regression model dependent on Poisson process}
\usage{
\S4method{estimate}{jumpRegression}(model.class, t, data, nMCMC, propSd,
  adapt = TRUE, proposal = c("normal", "lognormal"), it.xi = 10)
}
\arguments{
\item{model.class}{class of the regression model based on the NHPP including all required information, see \code{\link{jumpRegression-class}}}

\item{t}{vector of time points}

\item{data}{vector of observation variables}

\item{nMCMC}{length of Markov chain}

\item{propSd}{vector of proposal variances for \eqn{(\theta, \xi)}}

\item{adapt}{if TRUE (default), proposal variance is adapted}

\item{proposal}{proposal density for \eqn{\theta}: "normal" (default) or "lognormal" (for positive parameters)}

\item{it.xi}{number of iterations for MH step for \eqn{\xi} inside the Gibbs sampler}
}
\description{
Bayesian estimation of the parameter of the regression model
  \eqn{y_i = f(t_i, N_{t_i}, \theta) + \epsilon_i} with
  \eqn{N_t\sim Pois(\Lambda(t, \xi)), \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t))}.
}
\section{Proposal densities}{

For \eqn{\theta}, there is the possibility to choose "normal" or "lognormal".
The proposal density for \eqn{\xi} depends on the starting value of \eqn{\xi}. If all components are positive, the proposal density is lognormal, and normal otherwise.
}
\examples{
t <- seq(0,1, by = 0.01)
model <- set.to.class("jumpRegression", fun = function(t, N, theta) exp(theta[1]*t) + theta[2]*N,
                   parameter = list(theta = c(2, 2), gamma2 = 0.25, xi = c(3, 0.5)),
                   Lambda = function(t, xi) (t/xi[2])^xi[1])
data <- simulate(model, t = t, plot.series = FALSE)
est <- estimate(model, t, data, 1000)
plot(est)
\dontrun{
# work in progress
est_hid <- estimate(model, t, data$Y, 1000)
plot(est_hid)
}
}
\references{
Heeke, G., S. Hermann, R. Maurer, K. Ickstadt, and C. H. Mueller (2015).
Stochastic Modeling and Statistical Analysis of Fatigue Tests on Prestressed Concrete Beams under Cyclic Loadings.
SFB 823 discussion paper 25/15.
}

