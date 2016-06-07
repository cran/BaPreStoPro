#' Building of model classes
#'
#' @description Definition of the model classes.
#' @param class.name name of model class
#' @param parameter list of parameter values
#' @param prior optional list of prior parameters
#' @param start optional list of starting values
#' @param b.fun drift function b
#' @param s.fun variance function s
#' @param h.fun jump high function h
#' @param sT.fun variance function \eqn{\widetilde{s}}
#' @param y0.fun function for the starting point, if dependent on parameter
#' @param fun regression function
#' @param Lambda intensity rate of Poisson process
#' @param priorDensity list of functions for prior densities, if missing: non-informative estimation
#' @section Details:
#' \code{set.to.class} is the central function to define a S4 model class, where the \code{simulate} and the \code{estimate} methods build up.
#' Main input parameter is \code{class.name}, which is one out of "jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion", "jumpRegression", "NHPP", "Regression" and "mixedRegression", 
#' which is the name of the class object containing all information of the model.
#' If you write \code{set.to.class(class.name)} without any further input parameter, the function tells you which entries the list \code{parameter} has to contain.
#' This is the second central input parameter. If input parameter \code{start} is missing, it is set to \code{parameters}. 
#' If input parameter \code{prior}, which is a list of prior parameters, is missing, they are calculated from \code{parameter} in that way, that prior mean and standard deviation is equal to the entries of \code{parameter}.
#' Functions \code{b.fun, s.fun, h.fun} can be seen in the model definition of the jump diffusion \eqn{dY_t = b(\phi, t, Y_t)dt + s(\gamma^2, t, Y_t)dW_t + h(\theta, t, Y_t)dN_t}.
#' In the case of a continuous diffusion, one out of "Diffusion", "mixedDiffusion", "hiddenDiffusion" or "hiddenmixedDiffusion", variance function \eqn{s(\gamma^2, t, y)} is restricted to the case \eqn{s(\gamma^2, t, y)=\gamma\widetilde{s}(t, y)}. \code{sT.fun} stands for \eqn{\widetilde{s}(t, y)}.
#' In the case of a regression model, "Regression" or "mixedRegression", \code{sT.fun} means the variance function dependent on t of the regression error \eqn{\epsilon_i\sim N(0,\sigma^2\widetilde{s}(t))}.
#' In both cases, default value is \code{sT.fun = function(t, y) 1}.
#' \code{y0.fun} is for the models, where the starting value depends on the parameter phi, "mixedDiffusion", "hiddenDiffusion" or "hiddenmixedDiffusion". Default value is a constant function in 1.
#' \code{fun} is the regression function for the models "Regression", "mixedRegression" and "jumpRegression". In the first two cases, this is \eqn{f(\phi, t)} and in the third \eqn{f(t, N_t, \theta)}.
#' Function \code{Lambda} is the cumulative intensity function in the models including the non-homogeneous Poisson process.
#' Input parameter \code{priorDensity} is for the model class \code{\link{jumpDiffusion-class}} a list of functions for the prior density functions. 
#' For the model classes \code{\link{NHPP-class}} and \code{\link{Merton-class}}, \code{priorDensity} is the density of the intensity rate parameter of the Poisson process.
#' Default is a non-informative approach for all cases.
#' 
#' @examples
#' set.to.class("jumpDiffusion")
#' (names <- set.to.class("jumpDiffusion"))
#' model <- set.to.class("jumpDiffusion", 
#'              parameter = list(theta = 0.1, phi = 0.01, gamma2 = 0.1, xi = 3))
#' summary(class.to.list(model))
#'

#' @export
set.to.class <- function(class.name = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion", "jumpRegression", "NHPP", "Regression", "mixedRegression"),
                         parameter, prior, start, b.fun, s.fun, h.fun, sT.fun, y0.fun, fun, Lambda, priorDensity){
  class.name <- match.arg(class.name)
  df <- defaults(class.name)
  if(missing(parameter)){
    message(paste("parameter has to be list of", toString(df)))
    return(invisible(df))
  }

  if (!is.list(parameter)) 
    stop(
      paste("parameter has to be list of", toString(df)) 
    )
  if (!missing(prior) && !all(getPriorNames(class.name) %in% names(prior)) && !(class.name %in% c("jumpDiffusion", "NHPP"))) 
    stop(
      paste("prior has to be list of", toString(getPriorNames(class.name)))
    )
  if (!missing(start) && !all(df %in% names(start))) 
    stop(
      paste("start has to be a list of", toString(df))
    )
  

  if(missing(prior)) prior <- getPrior(parameter, class.name)
  if(missing(start)) start <- parameter

  if(missing(b.fun)) b.fun <- function(phi, t, y) phi
  if(missing(s.fun)) s.fun <- function(gamma2, t, y) sqrt(gamma2)
  if(missing(h.fun)) h.fun <- function(theta, t, y) theta

  if(missing(sT.fun)) sT.fun <- function(t, y) 1
  if(missing(y0.fun)) y0.fun <- function(phi, t) 1
  if(missing(fun)){
    if(class.name %in% c("Regression", "mixedRegression")) fun1 <- function(phi, t) phi
    if(class.name == "jumpRegression") fun2 <- function(t, N, theta) theta
  }else{
    if(class.name %in% c("Regression", "mixedRegression")) fun1 <- fun
    if(class.name == "jumpRegression") fun2 <- fun
  }

  if(missing(Lambda)) Lambda <- function(t, xi) xi*t


  if(class.name == "jumpDiffusion"){
    if(missing(priorDensity)){
      priorDensity <- list(
        phi = function(phi) 1,
        theta = function(theta) 1,
        gamma2 = function(gamma2) 1,
        xi = function(xi_drawn) 1)
    }
    return(new(Class = class.name, theta = parameter$theta, phi = parameter$phi, gamma2 = parameter$gamma2, xi = parameter$xi,
               b.fun = b.fun, s.fun = s.fun, h.fun = h.fun, Lambda = Lambda, priorDensity = priorDensity,
               start = start))
  }
  if(class.name == "Merton"){
    if(missing(priorDensity)) priorDensity <- function(xi) 1
    return(new(Class = class.name, thetaT = parameter$thetaT, phi = parameter$phi, gamma2 = parameter$gamma2, xi = parameter$xi,
               Lambda = Lambda, prior = prior, start = start, priorDensity = priorDensity))
  }
  if(class.name == "Diffusion"){
    return(new(Class = class.name, phi = parameter$phi, gamma2 = parameter$gamma2,
                 b.fun = b.fun, sT.fun = sT.fun,
                 prior = prior, start = start))
  }
  if(class.name == "mixedDiffusion"){
    return(new(Class = class.name, phi = parameter$phi, mu = parameter$mu, Omega = parameter$Omega, gamma2 = parameter$gamma2,
                 y0.fun = y0.fun, b.fun = b.fun, sT.fun = sT.fun,
                 prior = prior, start = start))
  }
  if(class.name == "hiddenDiffusion"){
    return(new(Class = class.name, phi = parameter$phi, gamma2 = parameter$gamma2, sigma2 = parameter$sigma2,
                 b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun,
                 prior = prior, start = start))
  }
  if(class.name == "hiddenmixedDiffusion"){
    return(new(Class = class.name, phi = parameter$phi, mu = parameter$mu, Omega = parameter$Omega, gamma2 = parameter$gamma2,
                 sigma2 = parameter$sigma2, b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun,
                 prior = prior, start = start))
  }
  if(class.name == "jumpRegression"){
    return(new(Class = class.name, theta = parameter$theta, gamma2 = parameter$gamma2, xi = parameter$xi, fun = fun2, Lambda = Lambda,
               sT.fun = sT.fun, prior = prior, start = start))
  }
  if(class.name == "NHPP"){
    if(missing(priorDensity)) priorDensity <- function(xi) 1
    return(new(Class = class.name, xi = parameter$xi, Lambda = Lambda, start = parameter$xi, priorDensity = priorDensity))
  }
  if(class.name == "Regression"){
    return(new(Class = class.name, phi = parameter$phi, gamma2 = parameter$gamma2,
               fun = fun1, sT.fun = sT.fun,
               prior = prior, start = start))
  }
  if(class.name == "mixedRegression"){
    return(new(Class = class.name, phi = parameter$phi, mu = parameter$mu, Omega = parameter$Omega, gamma2 = parameter$gamma2,
               fun = fun1, sT.fun = sT.fun,
               prior = prior, start = start))
  }

}

defaults <- function(class.name = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion", "jumpRegression", "NHPP", "Regression", "mixedRegression")){

  class.name <- match.arg(class.name)
  if(class.name == "jumpDiffusion"){
    name.vec <- c("theta", "phi", "gamma2", "xi")
  }
  if(class.name == "Merton"){
    name.vec <- c("thetaT", "phi", "gamma2", "xi")
  }
  if(class.name == "Diffusion"){
    name.vec <- c("phi", "gamma2")
  }
  if(class.name == "mixedDiffusion"){
    name.vec <- c("phi", "mu", "Omega", "gamma2")
  }
  if(class.name == "hiddenDiffusion"){
    name.vec <- c("phi", "gamma2", "sigma2")
  }
  if(class.name == "hiddenmixedDiffusion"){
    name.vec <- c("phi", "mu", "Omega", "gamma2", "sigma2")
  }
  if(class.name == "jumpRegression"){
    name.vec <- c("theta", "gamma2", "xi")
  }
  if(class.name == "NHPP"){
    name.vec <- "xi"
  }
  if(class.name == "Regression"){
    name.vec <- c("phi", "gamma2")
  }
  if(class.name == "mixedRegression"){
    name.vec <- c("phi", "mu", "Omega", "gamma2")
  }
  name.vec
}

getPriorNames <- function(model = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion",
                                          "jumpRegression", "NHPP", "Regression", "mixedRegression")){
  model <- match.arg(model)
  
  names <- NULL
  
  if(model == "Merton"){
    names <- c("m.phi", "v.phi", "m.thetaT", "v.thetaT", "alpha.gamma", "beta.gamma")
  }
  if(model=="Diffusion"){
    names <- c("m.phi", "v.phi", "alpha.gamma", "beta.gamma")
  }
  if(model=="mixedDiffusion"){
    names <- c("m.mu", "v.mu", "alpha.omega", "beta.omega", "alpha.gamma", "beta.gamma")
  }
  if(model=="hiddenDiffusion"){
    names <- c("m.phi", "v.phi", "alpha.gamma", "beta.gamma", "alpha.sigma", "beta.sigma")
  }
  if(model=="hiddenmixedDiffusion"){
    names <- c("m.mu", "v.mu", "alpha.omega", "beta.omega", "alpha.gamma", "beta.gamma", "alpha.sigma", "beta.sigma")
  }
  if(model =="jumpRegression"){
    names <- c("m.theta", "v.theta", "alpha.gamma", "beta.gamma")
  }
  if(model == "Regression"){
    names <- c("m.phi", "v.phi", "alpha.gamma", "beta.gamma")
  }
  if(model == "mixedRegression"){
    names <- c("m.mu", "v.mu", "alpha.omega", "beta.omega", "alpha.gamma", "beta.gamma")
  }
  names
}

#' Builds a list from class object
#'
#' @description Class slots are tranfered to list entries.
#' @param cl class object
#' @examples
#' model <- set.to.class("jumpDiffusion", 
#'              parameter = list(theta = 0.1, phi = 0.01, gamma2 = 0.1, xi = 3))
#' summary(class.to.list(model))
#' @export

class.to.list <- function(cl){
  
  sN <- slotNames(cl)
  res <- lapply(sN, function(name) slot(cl, name))
  names(res) <- sN
  res
}

getPrior <- function(parameter, model = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion",
                                          "jumpRegression", "NHPP", "Regression", "mixedRegression")){
  model <- match.arg(model)
  
  if(model=="jumpDiffusion"){
    prior <- list()
  }
  if(model == "Merton"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi, m.thetaT = parameter$thetaT, v.thetaT = parameter$thetaT,
                  alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="Diffusion"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, beta.gamma = 2*parameter$gamma2)
  }
  if(model=="mixedDiffusion"){
    prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                  beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="hiddenDiffusion"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2, alpha.sigma=3, beta.sigma=parameter$sigma2*2)
  }
  if(model=="hiddenmixedDiffusion"){
    prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                  beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2, alpha.sigma = 3, beta.sigma = parameter$sigma2*2)
    
  }
  if(model =="jumpRegression"){
    prior <- list(m.theta = parameter$theta, v.theta = parameter$theta^2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="NHPP"){
    prior <- list()
  }
  if(model == "Regression"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, beta.gamma = 2*parameter$gamma2)
  }
  if(model == "mixedRegression"){
    prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                  beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  prior
}

