

#' S4 class of model informations for the jump diffusion process
#' @description Informations of model 
#' \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma^2,t,Y_t)dW_t + h(\theta,t,Y_t)dN_t} with 
#' \eqn{N_t\sim Pois(\Lambda(t, \xi))}.
#' @slot theta parameter \eqn{\theta}
#' @slot phi parameter \eqn{\phi}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot xi parameter \eqn{\xi}
#' @slot b.fun function \eqn{b(\phi,t,y)}
#' @slot s.fun function \eqn{s(\gamma^2,t,y)}
#' @slot h.fun function \eqn{b(\theta,t,y)}
#' @slot Lambda function \eqn{\Lambda(t,\xi)}
#' @slot priorDensity list of prior density functions, default is a non-informative approach
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(phi = 0.01, theta = 0.1, gamma2 = 0.01, xi = c(2, 0.2))
#' b.fun <- function(phi, t, y) phi * y
#' s.fun <- function(gamma2, t, y) sqrt(gamma2) * y
#' h.fun <- function(theta, t, y) theta * y
#' Lambda <- function(t, xi) (t / xi[2])^xi[1]
#' priorDensity <- list(
#'   phi = function(phi) 1,
#'   theta = function(theta) dnorm(theta, 0.1, 0.001),
#'   gamma2 = function(gamma2) dgamma(1/gamma2, 3, 0.01*2),
#'   xi = function(xi) dgamma(xi, c(2, 0.2), 1)
#' )
#' start <- parameter
#' model <- set.to.class("jumpDiffusion", parameter, start = start, 
#'   b.fun = b.fun, s.fun = s.fun, h.fun = h.fun, Lambda = Lambda, 
#'   priorDensity = priorDensity)

setClass(Class = "jumpDiffusion", representation = representation(theta = "numeric", phi = "numeric", gamma2 = "numeric", xi = "numeric",
                                                                  b.fun = "function", s.fun = "function", h.fun = "function",
                                                                  Lambda = "function", priorDensity = "list", start = "list"))


#' S4 class of model informations for a special jump diffusion process, called Merton model
#' @description Informations of model 
#' \eqn{dY_t = \phi Y_t dt + \gamma^2 Y_t dW_t + \theta Y_tdN_t} with 
#' \eqn{N_t\sim Pois(\Lambda(t, \xi))}. The explicit solution of the SDE is given by 
#'   \eqn{Y_t = y_0 \exp( \phi t - \gamma^2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @slot thetaT parameter \eqn{\widetilde{\theta}=\log(1+\theta)}
#' @slot phi parameter \eqn{\phi}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot xi parameter \eqn{\xi}
#' @slot Lambda function \eqn{\Lambda(t,\xi)}
#' @slot prior list of prior parameters for \eqn{\phi, \widetilde{\theta}, \gamma^2}
#' @slot priorDensity list of prior density function for \eqn{\xi}
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(phi = 0.01, thetaT = 0.1, gamma2 = 0.01, xi = c(2, 0.2))
#' Lambda <- function(t, xi) (t / xi[2])^xi[1]
#' # prior density for xi:
#' priorDensity <- function(xi) dgamma(xi, c(2, 0.2), 1)
#' # prior parameter for phi (normal), thetaT (normal) and gamma2 (inverse gamma):
#' prior <- list(m.phi = parameter$phi, v.phi = parameter$phi, m.thetaT = parameter$thetaT, 
#'    v.thetaT = parameter$thetaT, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
#' start <- parameter
#' model <- set.to.class("Merton", parameter, prior, start, Lambda = Lambda, 
#'    priorDensity = priorDensity)
#' summary(class.to.list(model))
#' # default:
#' model <- set.to.class("Merton", parameter, Lambda = Lambda)
#' 
setClass(Class = "Merton", representation = representation(thetaT = "numeric", phi = "numeric", gamma2 = "numeric", xi = "numeric",
                                                           Lambda = "function", prior = "list", start = "list", priorDensity = "function"))
#' S4 class of model informations for diffusion process
#' @description Informations of model 
#' \eqn{dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t}.
#' @slot phi parameter \eqn{\phi}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot b.fun function \eqn{b(\phi,t,y)}
#' @slot sT.fun function \eqn{\widetilde{s}(t,y)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(phi = 0.1, gamma2 = 0.01)
#' b.fun <- function(phi, t, y) phi * y
#' sT.fun <- function(t, y) y
#' start <- parameter
#' prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, 
#'    alpha.gamma = 3, beta.gamma = 2*parameter$gamma2)
#' model <- set.to.class("Diffusion", parameter, prior, start, 
#'   b.fun = b.fun, sT.fun = sT.fun)

setClass(Class = "Diffusion", representation = representation(phi = "numeric", gamma2 = "numeric",
                                                              b.fun = "function", sT.fun = "function", prior = "list", start = "list"))
#' S4 class of model informations for hierarchical (mixed) diffusion process model
#' @description Informations of model 
#' \eqn{dY_t = b(\phi_j,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, \phi_j\sim N(\mu, \Omega), Y_{t_0}=y_0(\phi, t_0)}.
#' @slot phi parameter \eqn{\phi}
#' @slot mu parameter \eqn{\mu}
#' @slot Omega parameter \eqn{\Omega}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot y0.fun function \eqn{y_0(\phi, t)}
#' @slot b.fun function \eqn{b(\phi,t,y)}
#' @slot sT.fun function \eqn{\widetilde{s}(t,y)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' mu <- c(2, 1); Omega <- c(1, 0.04)
#' phi <- sapply(1:2, function(i) rnorm(21, mu[i], sqrt(Omega[i])))
#' parameter <- list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.01)
#' b.fun <- function(phi, t, y) phi[1] * y
#' sT.fun <- function(t, y) y
#' y0.fun <- function(phi, t) phi[2]
#' start <- parameter
#' prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, 
#'    alpha.omega = rep(3, length(parameter$mu)), beta.omega = parameter$Omega*2, 
#'    alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
#' model <- set.to.class("mixedDiffusion", parameter, prior, start, 
#'   b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun)

setClass(Class = "mixedDiffusion", representation = representation(phi = "matrix", mu = "numeric", Omega = "numeric", gamma2 = "numeric",
                                                                   y0.fun = "function", b.fun = "function", sT.fun = "function", prior = "list", start = "list"))

#' S4 class of model informations for hidden diffusion process
#' @description Informations of model 
#'   \eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, 
#'   \epsilon_i\sim N(0,\sigma^2), Y_{t_0}=y_0(\phi, t_0)}.
#' @slot phi parameter \eqn{\phi}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot sigma2 parameter \eqn{\sigma^2}
#' @slot y0.fun function \eqn{y_0(\phi, t)}
#' @slot b.fun function \eqn{b(\phi,t,y)}
#' @slot sT.fun function \eqn{\widetilde{s}(t,y)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(phi = c(2, 1), gamma2 = 0.1, sigma2 = 0.1)
#' b.fun <- function(phi, t, y) phi[1] * y
#' sT.fun <- function(t, y) y
#' y0.fun <- function(phi, t) phi[2]
#' start <- parameter
#' prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, 
#'    beta.gamma = parameter$gamma2*2, alpha.sigma=3, beta.sigma=parameter$sigma2*2)
#' model <- set.to.class("hiddenDiffusion", parameter, prior, start, 
#'   b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun)

setClass(Class = "hiddenDiffusion", representation = representation(phi = "numeric", gamma2 = "numeric",
                                                                    sigma2 = "numeric",
                                                                    b.fun = "function", sT.fun = "function", y0.fun = "function",
                                                                    prior = "list", start = "list"))
#' S4 class of model informations for hierarchical (mixed) hidden diffusion process
#' @description Informations of model 
#'   \eqn{Z_{ij} = Y_{t_{ij}} + \epsilon_{ij}, dY_t = b(\phi_j,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, \phi_j\sim N(\mu, \Omega), 
#'   Y_{t_0}=y_0(\phi, t_0), \epsilon_{ij}\sim N(0,\sigma^2)}.
#' @slot phi parameter \eqn{\phi}
#' @slot mu parameter \eqn{\mu}
#' @slot Omega parameter \eqn{\Omega}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot sigma2 parameter \eqn{\sigma^2}
#' @slot y0.fun function \eqn{y_0(\phi, t)}
#' @slot b.fun function \eqn{b(\phi,t,y)}
#' @slot sT.fun function \eqn{\widetilde{s}(t,y)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' mu <- c(2, 1); Omega <- c(1, 0.04)
#' phi <- sapply(1:2, function(i) rnorm(21, mu[i], sqrt(Omega[i])))
#' parameter <- list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1, sigma2 = 0.1)
#' b.fun <- function(phi, t, y) phi[1] * y
#' sT.fun <- function(t, y) y
#' y0.fun <- function(phi, t) phi[2]
#' start <- parameter
#' prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, 
#'    alpha.omega = rep(3, length(parameter$mu)), beta.omega = parameter$Omega*2, 
#'    alpha.gamma = 3, beta.gamma = parameter$gamma2*2, 
#'    alpha.sigma = 3, beta.sigma = parameter$sigma2*2)
#' model <- set.to.class("hiddenmixedDiffusion", parameter, prior, start, 
#'   b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun)

setClass(Class = "hiddenmixedDiffusion", representation = representation(phi = "matrix", mu = "numeric", Omega = "numeric", gamma2 = "numeric",
                                                                    sigma2 = "numeric",
                                                                    b.fun = "function", sT.fun = "function", y0.fun = "function",
                                                                    prior = "list", start = "list"))
#' S4 class of model informations for the jump regression model
#' @description Informations of model 
#'   \eqn{y_i = f(t_i, N_{t_i}, \theta) + \epsilon_i} with
#'   \eqn{N_t\sim Pois(\Lambda(t, \xi)), \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t))}.
#' @slot theta parameter \eqn{\theta}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot xi parameter \eqn{\xi}
#' @slot fun function \eqn{f(t, N, \theta)}
#' @slot sT.fun function \eqn{\widetilde{s}(t)}
#' @slot Lambda function \eqn{\Lambda(t,\xi)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(theta = c(3, 1), gamma2 = 0.1, xi = c(2, 0.2))
#' fun <- function(t, N, theta) theta[1]*t + theta[2]*N
#' sT.fun <- function(t) t
#' Lambda <- function(t, xi) (t / xi[2])^xi[1]
#' prior <- list(m.theta = parameter$theta, v.theta = parameter$theta^2, 
#'    alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
#' start <- parameter
#' model <- set.to.class("jumpRegression", parameter, prior, start = start, 
#'   fun = fun, sT.fun = sT.fun, Lambda = Lambda)

setClass(Class = "jumpRegression", representation = representation(theta = "numeric", gamma2 = "numeric", xi = "numeric", fun = "function",
                                                                   Lambda = "function", sT.fun = "function", prior = "list", start = "list"))

#' S4 class of model informations for non-homogeneous Poisson process
#' @description Informations of NHPP with cumulative intensity function 
#' \eqn{\Lambda(t, \xi)}. 
#' @slot xi parameter \eqn{\xi}
#' @slot Lambda function \eqn{\Lambda(t,\xi)}
#' @slot priorDensity prior density function for \eqn{\xi}
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(xi = c(2, 0.2))
#' Lambda <- function(t, xi) (t / xi[2])^xi[1]
#' priorDensity <- function(xi) dgamma(xi, c(2, 0.2), 1)
#' start <- parameter
#' model <- set.to.class("NHPP", parameter, start = start, Lambda = Lambda, 
#'    priorDensity = priorDensity)

setClass(Class = "NHPP", representation = representation(xi = "numeric", Lambda = "function", priorDensity = "function", start = "numeric"))

#' S4 class of model informations for the regression model
#' @description Informations of model 
#'   \eqn{y_i = f(\phi, t_i) + \epsilon_i, \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t_i))}.
#' @slot phi parameter \eqn{\phi}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot fun function \eqn{f(\phi, t)}
#' @slot sT.fun function \eqn{\widetilde{s}(t)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' parameter <- list(phi = c(3, 1), gamma2 = 0.1)
#' fun <- function(phi, t) phi[1] + phi[2]*t
#' sT.fun <- function(t) t
#' prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, 
#'    alpha.gamma = 3, beta.gamma = 2*parameter$gamma2)
#' start <- parameter
#' model <- set.to.class("Regression", parameter, prior, start, fun = fun, sT.fun = sT.fun)

setClass(Class = "Regression", representation = representation(phi = "numeric", gamma2 = "numeric",
                                                              fun = "function", sT.fun = "function", prior = "list", start = "list"))
#' S4 class of model informations for the hierarchical (mixed) regression model
#' @description Informations of model 
#'   \eqn{y_{ij} = f(\phi_j, t_{ij}) + \epsilon_{ij}, \phi_j\sim N(\mu, \Omega),
#'   \epsilon_{ij}\sim N(0,\gamma^2\widetilde{s}(t_{ij}))}.
#' @slot phi parameter \eqn{\phi}
#' @slot mu parameter \eqn{\mu}
#' @slot Omega parameter \eqn{\Omega}
#' @slot gamma2 parameter \eqn{\gamma^2}
#' @slot fun function \eqn{f(\phi, t)}
#' @slot sT.fun function \eqn{\widetilde{s}(t)}
#' @slot prior list of prior parameters
#' @slot start list of starting values for the Metropolis within Gibbs sampler
#' @examples 
#' mu <- c(2, 1); Omega <- c(1, 0.04)
#' phi <- sapply(1:2, function(i) rnorm(21, mu[i], sqrt(Omega[i])))
#' parameter <- list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.01)
#' fun <- function(phi, t) phi[1] + phi[2]*t
#' sT.fun <- function(t) t
#' prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, 
#'    alpha.omega = rep(3, length(parameter$mu)), beta.omega = parameter$Omega*2, 
#'    alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
#' start <- parameter
#' model <- set.to.class("mixedRegression", parameter, prior, start, fun = fun, sT.fun = sT.fun)

setClass(Class = "mixedRegression", representation = representation(phi = "matrix", mu = "numeric", Omega = "numeric", gamma2 = "numeric",
                                                                   fun = "function", sT.fun = "function", prior = "list", start = "list"))

####################
# estimation classes
####################

setClass(Class = "est.jumpDiffusion", representation = representation(theta = "numeric", phi = "numeric", gamma2 = "numeric", xi = "matrix",
                                                                      model = "list", N.est = "matrix", t = "numeric",
                                                                      Y = "numeric", N = "numeric",
                                                                      burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.Merton", representation = representation(thetaT = "numeric", phi = "numeric", gamma2 = "numeric", xi = "matrix",
                                                               model = "list", N.est = "matrix", t = "numeric",
                                                               Y = "numeric", N = "numeric",
                                                               burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.Diffusion", representation = representation(phi = "matrix", gamma2 = "numeric",
                                                                  model = "list", t = "numeric", Y = "numeric",
                                                                  burnIn = "numeric", thinning = "numeric"))
# Y, t as list ?!?
setClass(Class = "est.mixedDiffusion", representation = representation(phi = "list", mu = "matrix", Omega = "matrix", gamma2 = "numeric",
                                                                       model = "list", t = "numeric", Y = "matrix", t.list = "list", Y.list = "list",
                                                                       burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.hiddenDiffusion", representation = representation(phi = "matrix", gamma2 = "numeric", sigma2 = "numeric", Y.est = "matrix",
                                                                        model = "list", t = "numeric", Z = "numeric",
                                                                        burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.hiddenmixedDiffusion", representation = representation(phi = "list", mu = "matrix", Omega = "matrix", gamma2 = "numeric",
                                                                         sigma2 = "numeric", Y.est = "list",
                                                                         model = "list", t = "numeric", Z = "matrix", t.list = "list", Z.list = "list",
                                                                         burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.jumpRegression", representation = representation(theta = "matrix", gamma2 = "numeric", xi = "matrix", N.est = "matrix",
                                                                       model = "list", t = "numeric", Y = "numeric", N = "numeric",
                                                                       burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.NHPP", representation = representation(xi = "matrix", model = "list", N = "numeric", t = "numeric", jumpTimes = "numeric",
                                                             burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.Regression", representation = representation(phi = "matrix", gamma2 = "numeric",
                                                                   model = "list", t = "numeric", Y = "numeric",
                                                                   burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.mixedRegression", representation = representation(phi = "list", mu = "matrix", Omega = "matrix", gamma2 = "numeric",
                                                                        model = "list", t = "numeric", Y = "matrix", t.list = "list", Y.list = "list",
                                                                        burnIn = "numeric", thinning = "numeric"))






