\name{BaPreStoPro-package}
\alias{BaPreStoPro-package}
\alias{BaPreStoPro}
\docType{package}
\title{
Bayesian Prediction of Stochastic Processes
}
\description{
This package contains simulate, estimate and predict methods for non-homogeneous Poisson processes (NHPP), jump diffusions, (mixed) diffusions, hidden (mixed) diffusion models, regression model including a NHPP, and (mixed) regression models for comparison.
}
\details{
\tabular{ll}{
Package: \tab BaPreStoPro\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-06-07\cr
License: \tab GLP-2, GLP-3\cr
}
Each of the models has its specific class, \code{"jumpDiffusion"}, \code{"Merton"}, \code{"Diffusion"}, \code{"mixedDiffusion"}, \code{"hiddenDiffusion"}, \code{"hiddenmixedDiffusion"}, \code{"jumpRegression"}, \code{"NHPP"}, \code{"Regression"}, \code{"mixedRegression"}, created with the function \code{\link{set.to.class}}. For each of the model classes, a method \code{simulate} and \code{estimate} are provided. The output of method \code{estimate} is a new class object with the prefix "est.". For the estimation classes, methods \code{plot} and \code{predict} are available. An overview of the package can be found in Hermann (2016a) and theoretical details to the prediction procedures in Hermann (2016b).
}
\section{jumpDiffusion}{
The model class  \code{"jumpDiffusion"} contains all information about model defined by the stochastic differential equation (SDE) \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma^2,t,Y_t)dW_t + h(\theta,t,Y_t)dN_t} with
\eqn{N_t\sim Pois(\Lambda(t, \xi))} a non-homogeneous Poisson process and \eqn{W_t} a Brownian motion. The SDE is approximated with the Euler Maruyama approximation, which leads, dependent on the Poisson process variables, to a normal likelihood. 
For more information how to build the class, see examples in \code{\link{jumpDiffusion-class}} and \code{\link{estimate,jumpDiffusion-method}}. 

Estimation is done by a Metropolis-within-Gibbs sampler. For each of the parameters, a Metropolis-Hastings (MH) step is made, where the proposal density can be chosen between normal and lognormal. A proposal standard deviation can be chosen, which, if desired, is adapted after every 50 iterations, see Rosenthal (2011). 

In the case of unobserved variables of the Poisson process, one step of the Gibbs sampler is filtering of the unobserved jump process. Details can be found in Hermann and Ruggeri (2016) or Hermann (2016a).
}
\section{Merton}{
Specific choices of functions \eqn{b, s} and \eqn{h} lead to an explicit solution of the process \eqn{Y_t = y_0 \exp( \phi t - \gamma^2/2 t+\gamma W_t + \log(1+\theta) N_t)}. This model is well-known in the literature as Merton model. There are conjugate prior distributions available for \eqn{\phi, \gamma^2} and \eqn{\log(1+\theta)}: \eqn{\phi} and \eqn{\log(1+\theta)} are assumed to have normal prior distribution with parameters \code{m.phi} (mean) and \code{v.phi} (variance), \code{m.thetaT} (mean) and \code{v.thetaT} (variance) respectively. \eqn{\gamma^2} is assumed to have an inverse gamma distribution with parameters \code{alpha.gamma} and \code{beta.gamma}. An example how to build the model class can be found in \code{\link{Merton-class}} and \code{\link{estimate,Merton-method}}.

Estimation is similar to the jump diffusion process based on the Euler approximation. The difference is that for the parameters with conjugate priors, no MH step is necessary and drawing from the full conditional posterior is possible.
}
\section{Diffusion}{
The special case of \eqn{h(\theta, t, y)=0} leads to a general diffusion process \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma^2,t,Y_t)dW_t}. We here restrict to the special case of \eqn{s(\gamma^2, t, y)=\sqrt{\gamma^2}\widetilde{s}(t, y)}, because a conjugate prior for \eqn{\gamma^2} is available (the inverse gamma with parameters \code{alpha.gamma} and \code{beta.gamma}) in this case. For \eqn{\phi}, a normal prior with parameters \code{m.phi} (mean) and \code{v.phi} (variance) is assumed. An example can be found in \code{\link{Diffusion-class}} and \code{\link{estimate,Diffusion-method}}.

A Gibbs sampler with an MH step for \eqn{\phi} and one step drawing from the full conditional of \eqn{\gamma^2} is implemented. For \eqn{\phi}, the proposal density can be chosen, \code{"normal"} or \code{"lognormal"} and the proposal standard deviation, with the option to adapt, as well.
}
\section{mixedDiffusion}{
The diffusion process is extended to a hierarchical model with \eqn{\phi} as random effect with normal mixture distribution with mean \eqn{\mu} and variance \eqn{\Omega} (diagonal matrix). For \eqn{\mu}, a normal prior (parameters: \code{m.mu} and \code{v.mu}) is conjugate. For each diagonal element of \eqn{\Omega}, an inverse gamma prior (parameters: \code{alpha.omega} and \code{beta.omega}) is conjugate. Further information can be found in \code{\link{mixedDiffusion-class}} and \code{\link{estimate,mixedDiffusion-method}}.

A Gibbs sampler with an MH step for each random effect and each one step drawing from the full conditionals of \eqn{\mu, \Omega} and \eqn{\gamma^2} is implemented.
}
\section{hiddenDiffusion}{
The same model as the diffusion above is taken, but with an added error:
\eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t,
  \epsilon_i\sim N(0,\sigma^2), Y_{t_0}=y_0(\phi, t_0)}.
The inverse gamma prior for \eqn{\sigma^2} is conjugate and is, therefore, implemented with parameters \code{alpha.sigma} and \code{beta.sigma}. 
The diffusion process is a latent variable and has also to be estimated. A conditional sequential Monte Carlo (SMC) approach is implemented. For further details see Andrieu et al. (2010).
Examples can be found in \code{\link{hiddenDiffusion-class}} and \code{\link{estimate,hiddenDiffusion-method}}
}
\section{hiddenmixedDiffusion}{
Here, the hidden diffusion model is extended to a hierarchical model with the random effect \eqn{\phi}, similar to the mixed diffusion model. 
In the Gibbs sampler, one step is filtering the unobserved diffusion process for each observed series. Based on this estimation, the random effects are estimated with an MH step.
With conjugate priors, full conditionals for \eqn{\mu, \Omega, \gamma^2} and \eqn{\sigma^2} are available. 
Examples can be found in \code{\link{hiddenmixedDiffusion-class}} and \code{\link{estimate,hiddenmixedDiffusion-method}}.
}
\section{jumpRegression}{
We here consider the regression model: \eqn{y_i = f(t_i, N_{t_i}, \theta) + \epsilon_i} dependent on the Poisson process
  \eqn{N_t\sim Pois(\Lambda(t, \xi))} with \eqn{\epsilon_i\sim N(0,\gamma^2\widetilde{s}(t))}.

Here, for the case of missing observations for the Poisson process variable, a filtering procedure based on the conditional SMC is implemented. But this is still work in progress.

See examples in \code{\link{jumpRegression-class}} and \code{\link{estimate,jumpRegression-method}}.
}
\section{NHPP}{
Some models base on the non-homogeneous Poisson process. Here, only the NHPP itself is considered. A simple MH algorithm is implemented. A proposal density can be chosen, normal or lognormal, and the proposal standard deviation as well, which is adapted, if desired.

See examples in \code{\link{NHPP-class}} and \code{\link{estimate,NHPP-method}}.
}
\section{Regression}{
For the case of a comparison of regression and diffusion model, as made, for example, in Hermann et al. (2016), estimation and prediction is also made for regression models. All notations are analogously to the diffusion model. 
See examples in \code{\link{Regression-class}} and \code{\link{estimate,Regression-method}}.
}
\section{mixedRegression}{
Analogous hierarchical regression model to the mixed diffusion model above. 
See examples in \code{\link{mixedRegression-class}} and \code{\link{estimate,mixedRegression-method}}.
}
\author{
Simone Hermann
<hermann@statistik.tu-dortmund.de>
}
\references{
Hermann, S. (2016a). BaPreStoPro: an R Package for Bayesian Prediction of Stochastic Processes. SFB 823 discussion paper 28/16.

Hermann, S. (2016b). Bayesian Prediction for Stochastic Processes based on the Euler Approximation Scheme. SFB 823 discussion paper 27/16.

Hermann, S., K. Ickstadt and C. H.  Mueller (2016). Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. \emph{Applied Stochastic Models in Business and Industry, DOI: 10.1002/asmb.2175}.

Hermann, S. and F. Ruggeri (2016). Modelling Wear in Cylinder Liners. SFB 823 discussion paper 06/16. 

Hermann, S., K. Ickstadt, and C. H. Mueller (2015). Bayesian Prediction for a Jump Diffusion
Process with Application to Crack Growth in Fatigue Experiments. SFB 823 discussion paper 30/15. 

Heeke, G., S. Hermann, R. Maurer, K. Ickstadt, and C. H. Mueller (2015). Stochastic Modeling
and Statistical Analysis of Fatigue Tests on Prestressed Concrete Beams under Cyclic Loadings. SFB 823 discussion paper 25/15.

Monte Carlo methods:

Robert, C. P. and G. Casella (2004). Monte Carlo Statistical Methods. Springer, New York.

Adaptive MCMC:

Rosenthal, J. S. (2011). Optimal Proposal Distributions and Adaptive MCMC. In: Handbook of Markov Chain Monte Carlo, pp. 93-112.

particel Gibbs / SMC:

Andrieu, C., A. Doucet and R. Holenstein (2010). Particle Markov Chain Monte Carlo Methods.
Journal of the Royal Statistical Society B 72, pp. 269-342.

}
\keyword{package, (jump) diffusion, mixed (hidden) diffusion, non-homogeneous Poisson process}

\examples{
model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
t <- seq(0, 1, by = 0.1)
data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
est <- estimate(model, t, data, 10) # better: 10000
plot(est)
pred <- predict(est)

}
