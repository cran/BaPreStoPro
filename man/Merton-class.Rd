% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/class_definitions.R
\docType{class}
\name{Merton-class}
\alias{Merton-class}
\title{S4 class of model informations for a special jump diffusion process, called Merton model}
\description{
Informations of model
\eqn{dY_t = \phi Y_t dt + \gamma^2 Y_t dW_t + \theta Y_tdN_t} with
\eqn{N_t\sim Pois(\Lambda(t, \xi))}. The explicit solution of the SDE is given by
  \eqn{Y_t = y_0 \exp( \phi t - \gamma^2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
}
\section{Slots}{

\describe{
\item{\code{thetaT}}{parameter \eqn{\widetilde{\theta}=\log(1+\theta)}}

\item{\code{phi}}{parameter \eqn{\phi}}

\item{\code{gamma2}}{parameter \eqn{\gamma^2}}

\item{\code{xi}}{parameter \eqn{\xi}}

\item{\code{Lambda}}{function \eqn{\Lambda(t,\xi)}}

\item{\code{prior}}{list of prior parameters for \eqn{\phi, \widetilde{\theta}, \gamma^2}}

\item{\code{priorDensity}}{list of prior density function for \eqn{\xi}}

\item{\code{start}}{list of starting values for the Metropolis within Gibbs sampler}
}}
\examples{
parameter <- list(phi = 0.01, thetaT = 0.1, gamma2 = 0.01, xi = c(2, 0.2))
Lambda <- function(t, xi) (t / xi[2])^xi[1]
# prior density for xi:
priorDensity <- function(xi) dgamma(xi, c(2, 0.2), 1)
# prior parameter for phi (normal), thetaT (normal) and gamma2 (inverse gamma):
prior <- list(m.phi = parameter$phi, v.phi = parameter$phi, m.thetaT = parameter$thetaT,
   v.thetaT = parameter$thetaT, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
start <- parameter
model <- set.to.class("Merton", parameter, prior, start, Lambda = Lambda,
   priorDensity = priorDensity)
summary(class.to.list(model))
# default:
model <- set.to.class("Merton", parameter, Lambda = Lambda)
}

