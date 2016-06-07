#' Bayesian estimation
#'
#' @description Estimation method for the S4 classes.
#' @param model.class class object with model informations, see \code{\link{set.to.class}}
#' @param t vector or list of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#' @param ... parameters dependent on the model class
#' @return class object \code{est.}\code{model.class} containing Markov chains, data input and model informations
#' @references 
#' Hermann, S. (2016). BaPreStoPro: an R Package for Bayesian Prediction of Stochastic Processes. 
#' SFB 823 discussion paper 28/16.
#'
#' Robert, C. P. and G. Casella (2004). Monte Carlo Statistical Methods. Springer, New York.
#' 
#' Rosenthal, J. S. (2011). Optimal Proposal Distributions and Adaptive MCMC. In: Handbook of Markov Chain Monte Carlo, pp. 93-112.
#'
setGeneric("estimate", function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal"), ...) {
  standardGeneric("estimate")
})


########
#' Estimation for diffusion process
#'
#' @description Bayesian estimation of the parameters \eqn{\phi} and \eqn{\gamma^2} of the stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t}.
#' @param model.class class of the diffusion process model including all required information, see \code{\link{Diffusion-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\phi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#'
#' @examples
#' model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' est_diff <- estimate(model, t, data, 1000)
#' plot(est_diff)
#' 
#' @references 
#' Hermann, S., K. Ickstadt and C. H. Mueller (2016). 
#' Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. 
#' Applied Stochastic Models in Business and Industry, DOI: 10.1002/asmb.2175.
#' 
#' @export
setMethod(f = "estimate", signature = "Diffusion",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal")){
            proposal <- match.arg(proposal)
            
            
            if (!is.vector(t, mode = "numeric")) 
              stop(
                "t has to be a vector"
              )
            if (!is.vector(data, mode = "numeric"))
              stop(
                "data has to be a vector"
              )
            if (length(t) != length(data))
              stop(
                "t and data must have the same length"
              )
            if(!missing(propSd) && length(model.class@start$phi) != length(propSd))
              stop(
                "propSd must have length of phi"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start$phi < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )

            
            
            
    X <- data
    prior <- model.class@prior
    start <- model.class@start
    bSDE <- model.class@b.fun
    sVar <- model.class@sT.fun
    len <- nMCMC
    
    priorDensity <- function(phi) dnorm(phi, prior$m.phi, sqrt(prior$v.phi))
    
    if(missing(propSd)) propSd <- abs(prior$m.phi)/5
    lt <- length(t)
    dt <- t[-1] - t[-lt]
    lphi <- length(propSd)
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    
    postPhi <- function(lastPhi, gamma2, propSd){
      phi_old <- lastPhi
      phi_drawn <- proposals$draw(phi_old, propSd)
      ratio <- prod(priorDensity(phi_drawn) / priorDensity(phi_old))
      ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt)))
      ratio <- ratio* proposals$ratio(phi_drawn, phi_old, propSd)
      if(is.na(ratio)){ratio <- 0}
      if(runif(1) < ratio){
        phi_old <- phi_drawn
      }
      phi_old
    }
    postGamma2 <- function(phi){
      alphaPost <- prior$alpha.gamma + (lt-1)/2
      betaPost <-  prior$beta.gamma + sum( (X[-1] - X[-lt] - bSDE(phi, t[-lt], X[-lt])*dt)^2/(sVar(t[-lt], X[-lt])^2*dt) )/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    phi_out <- matrix(0, len, lphi)
    gamma2_out <- numeric(len)
    
    phi <- start$phi
    gamma2 <- start$gamma2
    
    for(count in 1:len){
      
      phi <- postPhi(phi, gamma2, propSd)
      gamma2 <- postGamma2(phi)
      
      phi_out[count, ] <- phi
      gamma2_out[count] <- gamma2
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:length(phi), function(i){
          ad.propSd(phi_out[(count-50+1):count, i], propSd[i], count/50) })
      }
      
      
    }
    result <- list(phi = phi_out, gamma2 = gamma2_out)
    
    if(nMCMC > 100){
      he <- matrix(0, ncol(result$phi) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$phi)+1)) he[i, ] <- diagnostic(result$phi[,i-1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
      
    } else {
      burnIn <- 0; thinning <- 1
    }

    result <- new(Class = "est.Diffusion", phi = result$phi, gamma2 = result$gamma2,
                  model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    return(result)

})



########
#' Estimation for hierarchical (mixed) diffusion model
#'
#' @description Bayesian estimation of a model
#' \eqn{dY_t = b(\phi_j,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, \phi_j\sim N(\mu, \Omega), Y_{t_0}=y_0(\phi, t_0)}.
#' @param model.class class of the hierarchical diffusion model including all required information, see \code{\link{mixedDiffusion-class}}
#' @param t list or vector of time points
#' @param data list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\phi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' model <- set.to.class("mixedDiffusion", 
#'              parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'              b.fun = function(phi, t, x) phi*x, sT.fun = function(t, x) x)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data[1:20,], 100)  # nMCMC should be much larger
#' plot(est)
#' 
#' # OU
#' b.fun <- function(phi, t, y) phi[1]-phi[2]*y; y0.fun <- function(phi, t) phi[3]
#' mu <- c(10, 5, 0.5); Omega <- c(0.9, 0.01, 0.01)
#' phi <- sapply(1:3, function(i) rnorm(21, mu[i], sqrt(Omega[i])))
#' model <- set.to.class("mixedDiffusion", 
#'                parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'                y0.fun = y0.fun, b.fun = b.fun, sT.fun = function(t, x) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data[1:20,], 100)  # nMCMC should be much larger
#' plot(est)
#' 
#' ##
#' t.list <- list()
#' for(i in 1:20) t.list[[i]] <- t
#' t.list[[21]] <- t[1:50]
#' data.list <- list()
#' for(i in 1:20) data.list[[i]] <- data[i,]
#' data.list[[21]] <- data[21, 1:50]
#' est <- estimate(model, t.list, data.list, 100)
#' pred <- predict(est, t = t[50:101], which.series = "current", ind.pred = 21, 
#'    b.fun.mat = function(phi, t, y) phi[,1]-phi[,2]*y)
#' @references 
#' Hermann, S., K. Ickstadt and C. H. Mueller (2016). 
#' Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. 
#' Applied Stochastic Models in Business and Industry, DOI: 10.1002/asmb.2175.
#' @export
setMethod(f = "estimate", signature = "mixedDiffusion",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal")) {
            proposal <- match.arg(proposal)
            
            if (!(is.vector(t))) 
              stop(
                "t has to be a vector or a list"
              )
            if (!(is.matrix(data) || is.list(data)))
              stop(
                "data has to be a matrix or a list"
              )
            if (is.vector(t, mode = "numeric") && (length(t) != nrow(data) && length(t) != ncol(data)))
              stop(
                "length of t has to be equal to the columns/rows of data"
              )
            if (is.list(t) && ( length(t) != length(data) || !all(sapply(t, length) == sapply(data, length))))
              stop(
                "data must match with t"
              )
            if(!missing(propSd) && ncol(model.class@start$phi) != length(propSd))
              stop(
                "propSd must have length of phi_j"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start$phi < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )

    prior <- model.class@prior
    start <- model.class@start
    y0.fun <- model.class@y0.fun
    bSDE <- model.class@b.fun
    sVar <- model.class@sT.fun
    len <- nMCMC

    if(is.matrix(data)){
      if(nrow(data) == length(t)){
        y <- t(data)
      }
      times <- list()
      y <- list()
      for(i in 1:nrow(data)){
        times[[i]] <- t
        y[[i]] <- data[i,]
      }
    }else{
      y <- data
      times <- t
    }

    postOm <- function(phi, mu){
      postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)
    }
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    
    postPhii_old <- function(lastPhi, mu, Omega, gamma2, X, t, propSd){  # constant y0
      lt <- length(t)
      dt <- diff(t)
      phi_old <- lastPhi
      phi_drawn <- proposals$draw(phi_old, propSd) 
      ratio <- prod(dnorm(phi_drawn, mu, sqrt(Omega)) / dnorm(phi_old, mu, sqrt(Omega)) )
      ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt)))
      ratio <- ratio* proposals$ratio(phi_drawn, phi_old, propSd)
      if(is.na(ratio)) ratio <- 0
      if(runif(1) <= ratio){
        phi_old <- phi_drawn
      }
      phi_old
    }
    postPhii <- function(lastPhi, mu, Omega, gamma2, X, t, propSd){  # X, t vektoren
      lt <- length(t)
      dt <- diff(t)
      phi_old <- lastPhi; lphi <- length(lastPhi)
      phi_drawn <- phi_old + rnorm(length(mu), 0, propSd)
      for(k in 1:lphi){
        
        phitest <- phi_drawn; phitest[k] <- rnorm(1, phitest[k], 0.1)
        if(y0.fun(phi_drawn, t[1]) != y0.fun(phitest, t[1])){
          fun <- function(theta){
            phi <- phi_drawn; phi[k] <- theta
            abs(y0.fun(phi, t[1]) - X[1])
          } 
          phi_drawn[k] <- optimize(f = fun, phi_drawn[k] + c(-1,1)*start$mu[k], maximum = FALSE)$minimum
        }
        
      }
      ratio <- prod(dnorm(phi_drawn, mu, sqrt(Omega)) / dnorm(phi_old, mu, sqrt(Omega)) )
      ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt)))
      if(is.na(ratio)) ratio <- 0
      if(runif(1) <= ratio){
        phi_old <- phi_drawn
      }
      phi_old
    }
    n <- length(y)
    N.all <- sum(sapply(y, length)-1)
    postGamma2 <- function(phi){
      alphaPost <- prior$alpha.gamma + N.all/2
      help <- numeric(n)
      for(i in 1:n){
        ni <- length(times[[i]])
        delta <- diff(times[[i]])
        help[i] <- sum( (y[[i]][-1] - y[[i]][-ni] - bSDE(phi[i,], times[[i]][-ni], y[[i]][-ni])*delta)^2/(sVar(times[[i]][-ni], y[[i]][-ni])^2*delta) )
      }
      betaPost <-  prior$beta.gamma + sum(help)/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    phi_out <- list()
    mu_out <- matrix(0, len, length(start$mu))
    Omega_out <- matrix(0, len, length(start$mu))
    gamma2_out <- numeric(len)
    
    phi <- start$phi
    gamma2 <- start$gamma2
    mu <- start$mu
    Omega <- postOm(phi, mu)
    
    if(missing(propSd)) propSd <- abs(start$mu)/5
    
    for(count in 1:len){
      
      for(i in 1:n){
        phi[i,] <- postPhii(phi[i,], mu, Omega, gamma2, y[[i]], times[[i]], propSd)
      }
      mu <- postmu(phi, prior$m.mu, prior$v.mu, Omega)
      Omega <- postOm(phi, mu)
      gamma2 <- postGamma2(phi)
      
      phi_out[[count]] <- phi
      mu_out[count,] <- mu
      Omega_out[count,] <- Omega
      gamma2_out[count] <- gamma2
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:length(phi[1,]), function(i){
          ad.propSd(sapply(phi_out[(count-50+1):count], function(mat) mat[1,i]), propSd[i], count/50) })
      }
      
    }
    result <- list(phi = phi_out, mu = mu_out, Omega = Omega_out, gamma2 = gamma2_out)
    
    if(nMCMC > 100){
      he <- matrix(0, ncol(result$mu) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$mu)+1)) he[i, ] <- diagnostic(result$mu[,i-1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    
    if(is.matrix(data)){
      result <- new(Class = "est.mixedDiffusion", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                    model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    }else{
      result <- new(Class = "est.mixedDiffusion", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                    model = class.to.list(model.class), t = t[[1]], Y = matrix(data[[1]], 1),
                    t.list = t, Y.list = data, burnIn = burnIn, thinning = thinning)
    }

    return(result)

})


########
#' Estimation for hidden diffusion process
#'
#' @description Bayesian estimation of the model
#'   \eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, 
#'   \epsilon_i\sim N(0,\sigma^2), Y_{t_0}=y_0(\phi, t_0)} with a particle Gibbs sampler.
#' @param model.class class of the hidden diffusion model including all required information, see \code{\link{hiddenDiffusion-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\phi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#' @param Npart number of particles in the particle Gibbs sampler
#'
#' @references 
#' Andrieu, C., A. Doucet and R. Holenstein (2010). Particle Markov Chain Monte Carlo Methods. 
#' Journal of the Royal Statistical Society B 72, pp. 269-342.
#' @examples
#' model <- set.to.class("hiddenDiffusion", y0.fun = function(phi, t) 0.5, 
#'              parameter = list(phi = 5, gamma2 = 1, sigma2 = 0.1))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data$Z, 100)  # nMCMC should be much larger!
#' plot(est)
#' 
#' \dontrun{
#' # OU
#' b.fun <- function(phi, t, y) phi[1]-phi[2]*y
#' model <- set.to.class("hiddenDiffusion", y0.fun = function(phi, t) 0.5, 
#'                parameter = list(phi = c(10, 1), gamma2 = 1, sigma2 = 0.1), 
#'                b.fun = b.fun, sT.fun = function(t, x) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data$Z, 1000)
#' plot(est)
#' }
#' @export
setMethod(f = "estimate", signature = "hiddenDiffusion",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal"), Npart = 100) {
            proposal <- match.arg(proposal)

            if (!is.vector(t, mode = "numeric")) 
              stop(
                "t has to be a vector"
              )
            if (!is.vector(data, mode = "numeric"))
              stop(
                "data has to be a vector"
              )
            if (length(t) != length(data))
              stop(
                "t and data must have the same length"
              )
            if(!missing(propSd) && length(model.class@start$phi) != length(propSd))
              stop(
                "propSd must have length of phi"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start$phi < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )
            
            
    y <- data
    prior <- model.class@prior
    start <- model.class@start
    len <- nMCMC
    sigmaTilde <- model.class@sT.fun
    y0.fun <- model.class@y0.fun
    b.fun <- model.class@b.fun
    maxIt <- 1   # input parameter ?
    
    if(missing(propSd)) propSd <- abs(start$phi)/5
    
    lt <- length(t)
    lphi <- length(start$phi)
    
    phi_out <- matrix(0,len,length(start$phi))
    sigma2_out <- numeric(len)
    gamma2_out <- numeric(len)
    X_out <- matrix(0,len,length(t))
    
    phi <- start$phi
    sigma2 <- start$sigma2
    gamma2 <- start$gamma2
    
    postSigma2 <- function(X){
      alphaPost <- prior$alpha.sigma + lt/2
      betaPost <-  prior$beta.sigma + sum((y-X)^2)/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    postGamma2 <- function(phi, X){
      alphaPost <- prior$alpha.gamma + (lt-1)/2
      help <- sum( (X[-1] - X[-lt] - b.fun(phi,t[-lt],X[-lt])*(t[-1]-t[-lt]))^2/(sigmaTilde(t[-lt],X[-lt])^2*(t[-1]-t[-lt])) )
      betaPost <-  prior$beta.gamma + help/2
      1/rgamma(1, alphaPost, betaPost)
    }
    priorDensity <- function(phi) dnorm(phi, prior$m.phi, sqrt(prior$v.phi))
    
    if(proposal == "lognormal"){
      if(any(start$phi < 0)) warning("Attention: proposal density has positive support")
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    
    estPhi_and_X <- function(X, lastPhi, gamma2, sigma2, B.fixed, propSd){
      likeli <- function(phi,X){
        prod(dnorm(X[-1],X[-lt]+b.fun(phi, t[-lt], X[-lt])*diff(t), sqrt(gamma2*diff(t))*sigmaTilde(t[-lt],X[-lt])))
      }
      # output variables
      phi_out <- lastPhi
      phi <- lastPhi
      for(k in 1:lphi){
        for(count in 1:maxIt){
          phi[k] <- proposals$draw(phi_out[k], propSd[k])
          ratio <- priorDensity(phi)[k]/priorDensity(phi_out)[k]
          ratio <- ratio * proposals$ratio(phi[k], phi_out[k], propSd[k])
          ratio <- ratio*likeli(phi, X)/likeli(phi_out, X)
          ratio <- ratio*dnorm(y[1], y0.fun(phi, t[1]), sqrt(sigma2))/dnorm(y[1], y0.fun(phi_out,t[1]), sqrt(sigma2))
          if(is.na(ratio)) ratio <- 0
          
          if(runif(1) <= ratio){
            phi_out[k] <- phi[k]
            break
          }
        }
      }
      res_SMC <- SMC(phi_out, gamma2, sigma2, Npart, t, y, b.fun, y0.fun, sigmaTilde, Y.fixed = X, B.fixed = B.fixed)
      lign.B <- A.to.B(res_SMC$parents)
      indice <- sample(1:Npart, 1, prob = res_SMC$W[,lt])
      X <- res_SMC$y[indice,]
      
      list(phi = phi_out, X = X, B.fixed = lign.B[indice,])
    }
    res_SMC <- SMC(phi, gamma2, sigma2, Npart, t, y, b.fun, y0.fun, sigmaTilde, conditional = FALSE)
    lign.B <- A.to.B(res_SMC$parents)
    indice <- sample(1:Npart, 1, prob =  res_SMC$W[,lt])
    X <- res_SMC$y[indice,]
    B.fixed <- lign.B[indice,]
    
    for(count in 1:len){
      
      # Particle Gibbs
      result <- estPhi_and_X(X, phi, gamma2, sigma2, B.fixed, propSd)
      phi <- result$phi
      X <- result$X
      B.fixed <- result$B.fixed
      
      gamma2 <- postGamma2(phi, X)
      sigma2 <- postSigma2(X)
      
      phi_out[count,] <- phi
      sigma2_out[count] <- sigma2
      gamma2_out[count] <- gamma2
      X_out[count,] <- X
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:lphi, function(i){
          ad.propSd(phi_out[(count-50+1):count, i], propSd[i], count/50) })
      }
      
      if (count%%1000 == 0) message(paste(count, "iterations done"))
      
    }
    result <- list(phi = phi_out, sigma2 = sigma2_out, gamma2 = gamma2_out, X = X_out)

    if(nMCMC > 100){
      he <- matrix(0, ncol(result$phi) + 2, 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$sigma2)
      for(i in 3:(ncol(result$phi)+2)) he[i, ] <- diagnostic(result$phi[,i-2])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    

    result <- new(Class = "est.hiddenDiffusion", phi = result$phi, gamma2 = result$gamma2, sigma2 = result$sigma2, Y.est = result$X,
                  model = class.to.list(model.class), t = t, Z = data, burnIn = burnIn, thinning = thinning)
    return(result)

})



########
#' Estimation for hierarchical (mixed) hidden diffusion process
#'
#' @description Bayesian estimation of the parameters in the hierarchical model: 
#'   \eqn{Z_{ij} = Y_{t_{ij}} + \epsilon_{ij}, dY_t = b(\phi_j,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, \phi_j\sim N(\mu, \Omega), 
#'   Y_{t_0}=y_0(\phi, t_0), \epsilon_{ij}\sim N(0,\sigma^2)} with the particle Gibbs sampler.
#' @param model.class class of the hierarchical hidden diffusion model including all required information, see \code{\link{hiddenmixedDiffusion-class}}
#' @param t list or vector of time points
#' @param data list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\phi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#' @param Npart number of particles in the particle Gibbs sampler
#' @references 
#' Andrieu, C., A. Doucet and R. Holenstein (2010). Particle Markov Chain Monte Carlo Methods. 
#' Journal of the Royal Statistical Society B 72, pp. 269-342.
#' @examples
#' mu <- c(5, 1); Omega <- c(0.9, 0.04)
#' phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' y0.fun <- function(phi, t) phi[2]
#' model <- set.to.class("hiddenmixedDiffusion", y0.fun = y0.fun, 
#'                  b.fun = function(phi, t, y) phi[1], 
#'                  parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 1, sigma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' 
#' \dontrun{
#' est <- estimate(model, t, data$Z[1:20,], 2000)
#' plot(est)
#' }
#' @export
setMethod(f = "estimate", signature = "hiddenmixedDiffusion",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal"), Npart = 100) {
            proposal <- match.arg(proposal)

            if (!(is.vector(t))) 
              stop(
                "t has to be a vector or a list"
              )
            if (!(is.matrix(data) || is.list(data)))
              stop(
                "data has to be a matrix or a list"
              )
            if (is.vector(t, mode = "numeric") && (length(t) != nrow(data) && length(t) != ncol(data)))
              stop(
                "length of t has to be equal to the columns/rows of data"
              )
            if (is.list(t) && ( length(t) != length(data) || !all(sapply(t, length) == sapply(data, length))))
              stop(
                "data must match with t"
              )
            if(!missing(propSd) && ncol(model.class@start$phi) != length(propSd))
              stop(
                "propSd must have length of phi_j"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start$phi < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )
            

    if(is.matrix(data)){
      if(nrow(data) == length(t)){
        y <- t(data)
      }else{
        if(ncol(data) != length(t)){
          stop("length of t has to be equal to the columns of y")
        }
      }
      times <- list()
      y <- list()
      for(i in 1:nrow(data)){
        times[[i]] <- t
        y[[i]] <- data[i,]
      }
    }else{
      y <- data
      times <- t
    }
    
    prior <- model.class@prior
    start <- model.class@start
    len <- nMCMC
    sigmaTilde <- model.class@sT.fun
    y0.fun <- model.class@y0.fun
    b.fun <- model.class@b.fun
    maxIt <- 1    # input parameter ?

    if(missing(propSd)) propSd <- abs(start$mu)/5
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    
    lphi <- length(start$mu)
    n <- length(y)
    n_all1 <- sum(sapply(y, length) - 1)
    n_all2 <- sum(sapply(y, length))
    
    postSigma2 <- function(X){   
      alphaPost <- prior$alpha.sigma + n_all2/2
      betaPost <-  prior$beta.sigma + sum((unlist(y)-unlist(X))^2)/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    postGamma2 <- function(X, phi){ #  phi matrix, X, t lists, alpha, beta prior parameters
      alphaPost <- prior$alpha.gamma + n_all1/2
      
      help <- numeric(n)
      for(i in 1:n){
        ni <- length(times[[i]])
        dt <- diff(times[[i]])
        help[i] <- sum( (X[[i]][-1] - X[[i]][-ni] - b.fun(phi[i,],times[[i]][-ni],X[[i]][-ni])*dt)^2/(sigmaTilde(times[[i]][-ni],X[[i]][-ni])^2*dt) )
      }
      betaPost <-  prior$beta.gamma + sum(help)/2
      
      1/rgamma(1, alphaPost, betaPost)
    }
    
    postPhii_Xi <- function(y, t, X, lastPhi, mu, Omega, gamma2, sigma2, B.fixed, propSd){
      lt <- length(t)
      likeli <- function(phi,X){
        prod(dnorm(X[-1], X[-lt]+b.fun(phi, t[-lt], X[-lt])*diff(t), sqrt(gamma2*diff(t))*sigmaTilde(t[-lt],X[-lt])))  # true transition density???
      }
      # output variables
      phi_out <- lastPhi
      phi <- lastPhi
      for(k in 1:lphi){
        for(count in 1:maxIt){
          phi[k] <- proposals$draw(phi_out[k], propSd[k])
          ratio <- dnorm(phi[k], mu[k], sqrt(Omega[k]))/dnorm(phi_out[k], mu[k], sqrt(Omega[k]))
          ratio <- ratio * proposals$ratio(phi[k], phi_out[k], propSd[k])
          ratio <- ratio*likeli(phi, X)/likeli(phi_out, X)
          ratio <- ratio*dnorm(y[1], y0.fun(phi, t[1]), sqrt(sigma2))/dnorm(y[1], y0.fun(phi_out, t[1]), sqrt(sigma2))
          if(is.na(ratio)) ratio <- 0
          
          if(runif(1) <= ratio){
            phi_out[k] <- phi[k]
            break
          }
        }
      }
      res_SMC <- SMC(phi_out, gamma2, sigma2, Npart, t, y, b.fun, y0.fun, sigmaTilde, Y.fixed = X, B.fixed = B.fixed)
      lign.B <- A.to.B(res_SMC$parents)
      indice <- sample(1:Npart, 1, prob = res_SMC$W[,lt])
      X <- res_SMC$y[indice,]
      
      list(phi = phi_out, X = X, B.fixed = lign.B[indice,])
      
    }
    
    n <- length(y)
    
    phi_out <- list()
    mu_out <- matrix(0, len, length(start$mu))
    Omega_out <- matrix(0, len, length(start$mu))
    sigma2_out <- numeric(len)
    gamma2_out <- numeric(len)
    X_out <- list()
    
    phi <- start$phi
    sigma2 <- start$sigma2
    gamma2 <- start$gamma2
    mu <- start$mu
    Omega <- postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)
    X <- list()
    B.fixed <- list()
    
    for(i in 1:n){
      
      if(dnorm(y[[i]][1], y0.fun(phi[i,], times[[i]][1]), sqrt(sigma2)) == 0){
        stop("bad starting values")
      }
      
      result <- SMC(phi[i,], gamma2, sigma2, Npart, times[[i]], y[[i]], b.fun, y0.fun, sigmaTilde, conditional = FALSE)
      lign.B <- A.to.B(result$parents)
      indice <- sample(1:Npart, 1, prob = result$W[,length(times[[i]])])
      X[[i]] <- result$y[indice,]
      B.fixed[[i]] <- lign.B[indice,]
    }
    
    for(count in 1:len){
      
      for(i in 1:n){
        help <- postPhii_Xi(y[[i]], times[[i]], X[[i]], phi[i,], mu, Omega, gamma2, sigma2, B.fixed[[i]], propSd)
        X[[i]] <- help$X
        phi[i,] <- help$phi
        B.fixed[[i]] <- help$B.fixed
      }
      
      mu <- postmu(phi, prior$m.mu, prior$v.mu, Omega)
      Omega <- postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)
      
      sigma2 <- postSigma2(X)
      gamma2 <- postGamma2(X, phi)
      
      phi_out[[count]] <- phi
      mu_out[count,] <- mu
      Omega_out[count,] <- Omega
      sigma2_out[count] <- sigma2
      gamma2_out[count] <- gamma2
      X_out[[count]] <- X
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:lphi, function(i){
          ad.propSd(sapply(phi_out[(count-50+1):count], function(mat) mat[1, i]), propSd[i], count/50) })
      }
      
      if (count%%100 == 0) message(paste(count, "iterations done"))
    }
    result <- list(phi = phi_out, mu = mu_out, Omega = Omega_out, sigma2 = sigma2_out, gamma2 = gamma2_out, X = X_out)
    
    
    if(nMCMC > 100){
      he <- matrix(0, ncol(result$mu) + 2, 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$sigma2)
      for(i in 3:(ncol(result$mu)+2)) he[i, ] <- diagnostic(result$mu[,i-2])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    

    if(is.matrix(data)){
      result <- new(Class = "est.hiddenmixedDiffusion", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                    sigma2 = result$sigma2, Y.est = result$X,
                    model = class.to.list(model.class), t = t, Z = data, burnIn = burnIn, thinning = thinning)
    }else{
      result <- new(Class = "est.hiddenmixedDiffusion", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                    sigma2 = result$sigma2, Y.est = result$X,
                    model = class.to.list(model.class), t = t[[1]], Z = matrix(data[[1]], 1),
                    t.list = t, Z.list = data, burnIn = burnIn, thinning = thinning)
    }

    return(result)

})





########
#' Estimation for a non-homogeneous Poisson process
#'
#' @description Bayesian estimation of a non-homogeneous Poisson process (NHPP) with cumulative intensity function \eqn{\Lambda(t, \xi)}.
#' @param model.class class of the NHPP model including all required information, see \code{\link{NHPP-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\xi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#' @references 
#' Hermann, S., K. Ickstadt and C. H. Mueller (2015). 
#' Bayesian Prediction for a Jump Diffusion Process with Application to Crack Growth in Fatigue Experiments.
#' SFB 823 discussion paper 30/15.
#' @examples
#' model <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)), 
#'                    Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data$Times, 10000, proposal = "lognormal")
#' plot(est)
#' 
#' ##
#' model <- set.to.class("NHPP", parameter = list(xi = 5), 
#'                    Lambda = function(t, xi) t*xi)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data$N, 10000)
#' plot(est, par.options = list(mfrow = c(1,1)))
#' @export
setMethod(f = "estimate", signature = "NHPP",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal")) {
            proposal <- match.arg(proposal)

            if (!is.vector(t, mode = "numeric"))
              stop(
                "t has to be a vector"
              )
            if (!is.vector(data, mode = "numeric"))
              stop(
                "data has to be a vector"
              )
            if(!missing(propSd) && ncol(model.class@start) != length(propSd))
              stop(
                "propSd must have length of xi"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )

    if(length(t) != length(data)){
      jumpTimes <- data
    }else{
      jumpTimes <- dNtoTimes(diff(data), t[-1])
    }
    Tend <- max(t)
    start <- model.class@start
    n <- nMCMC
    Lambda <- model.class@Lambda
    priorDensity <- model.class@priorDensity 
    
    lambda <- function(t, xi){
      h <- 1e-05
      (Lambda(t+h,xi)-Lambda(t,xi))/h
    }
    Lik <- function(xi){
      lambda_vec <- lambda(jumpTimes, xi)
      prod(lambda_vec)*exp(-Lambda(Tend, xi))
    }
    
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(xi_old, propSd){
        proposal(xi_old, propSd)
      }
      proposals$ratio <- function(xi_drawn, xi_old, propSd){
        proposalRatio(xi_old, xi_drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(xi_old, propSd){ 
        rnorm(length(xi_old), xi_old, propSd)
      }
      proposals$ratio <- function(xi_drawn, xi_old, propSd) 1
    }
    if(missing(propSd)) propSd <- (abs(start)+0.1)/2
    xi_old <- start
    xi_out <- matrix(0, length(start), n)
    LikOld <- Lik(xi_old)
    
    for(count in 1:n){
      xi_drawn <- proposals$draw(xi_old, propSd)
      
      LikNew <- Lik(xi_drawn)
      ratio <- proposals$ratio(xi_drawn, xi_old, propSd)
      ratio <- ratio*prod(priorDensity(xi_drawn)/priorDensity(xi_old))
      ratio <- ratio*LikNew/LikOld
      if(is.na(ratio)) ratio <- 0
      
      if(runif(1) <= ratio){
        xi_old <- xi_drawn
        LikOld <- LikNew
      }
      xi_out[,count] <- xi_old
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:length(start), function(i){
          ad.propSd(xi_out[i, (count-50+1):count], propSd[i], count/50) })
      }
      
    }
    res <- xi_out
          
    if(nMCMC > 100){
      if(is.vector(res)){
        he <- diagnostic(res); burnIn <- he[1]; thinning <- min( he[2], ceiling((nMCMC-burnIn)/100) )
      }else{
        he <- matrix(0, ncol(res), 2)
        for(i in 1:nrow(res)) he[i,] <- diagnostic(res[i,])
        burnIn <- max(he[, 1])
        thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
      }
    } else {
      burnIn <- 0; thinning <- 1
    }
    
    
            
  if(length(t) != length(data)){
    result <- new(Class = "est.NHPP", xi = t(res), jumpTimes = as.numeric(data), t = t, N = TimestoN(data, t),
                  model = class.to.list(model.class), burnIn = burnIn, thinning = thinning)

  } else {
    result <- new(Class = "est.NHPP", xi = t(res), t = t, N = data, jumpTimes = dNtoTimes(diff(data), t),
                  model = class.to.list(model.class), burnIn = burnIn, thinning = thinning)
  }
    return(result)
 })



########
#' Estimation for jump diffusion process
#'
#' @description Bayesian estimation of a stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma^2,t,Y_t)dW_t + h(\theta,t,Y_t)dN_t}.
#' @param model.class class of the jump diffusion model including all required information, see \code{\link{jumpDiffusion-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{(\phi, \theta, \gamma^2, \xi)}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density for phi, theta: "normal" (default) or "lognormal" (for positive parameters), see description below
#' @param it.xi number of iterations for MH step for \eqn{\xi} inside the Gibbs sampler
#' @section Proposal densities:
#' For \eqn{\gamma^2}, always the lognormal density is taken, since the parameter is always positive.
#' For \eqn{\theta} and \eqn{\phi}, there is the possibility to choose "normal" or "lognormal" (for both together). 
#' The proposal density for \eqn{\xi} depends on the starting value of \eqn{\xi}. If all components are positive, the proposal density is lognormal, and normal otherwise.
#' @examples
#' # non-informative
#' model <- set.to.class("jumpDiffusion", Lambda = function(t, xi) (t/xi[2])^xi[1],
#'                parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(model, t, data, 1000)
#' plot(est)
#' 
#' # informative
#' model <- set.to.class("jumpDiffusion", Lambda = function(t, xi) (t/xi[2])^xi[1],
#'    parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)),
#'    priorDensity = list(phi = function(phi) dnorm(phi, 0.05, 0.01),
#'                        theta = function(theta) dgamma(1/theta, 10, 0.1*9),
#'                        gamma2 = function(gamma2) dgamma(1/gamma2, 10, 0.1*9),
#'                        xi = function(xi) dnorm(xi, c(3, 1/4), c(1,1))))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(model, t, data, 1000)
#' plot(est)
#' 
#' \dontrun{
#' est_hidden <- estimate(model, t, data$Y, 1000)
#' plot(est_hidden)
#' }
#' @export
setMethod(f = "estimate", signature = "jumpDiffusion",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal"), it.xi = 5) {
    proposal <- match.arg(proposal)

    if (!is.vector(t, mode = "numeric")) 
      stop(
        "t has to be a vector"
      )
    if (!is.vector(data))
      stop(
        "data has to be a vector or a list with entries N and Y"
      )
    if (is.vector(data, mode = "numeric") && length(t) != length(data))
      stop(
        "t and data must have the same length"
      )
    if (is.list(data) &&  !(length(t) == length(data$Y) && length(t) == length(data$N)) ) 
      stop(
        "t, data$N and data$Y must have the same length"
      )
    if(!missing(propSd) && length(model.class@start$xi) + 3 != length(propSd))
      stop(
        "propSd must have length of xi + 3"
      )
    if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
      stop(
        "nMCMC has to be a natural number"
      )
    if(any(c(model.class@start$theta, model.class@start$phi) < 0) && proposal == "lognormal")
      stop(
        "lognormal proposal density has positive support and starting value is negative"
      )

    
    if(is.list(data)){  
      X <- data$Y
      N <- data$N
      jumpTimes <- dNtoTimes(diff(N), t[-1])
    }else{
      X <- data
    }  
    Tend <- max(t)
    n <- nMCMC
    start <- model.class@start
    Lambda <- model.class@Lambda
    lambda <- function(t, xi, h = 1e-05) (Lambda(t+h, xi)-Lambda(t, xi))/h
    Lik.N <- function(xi, jumpTimes){
      lambda_vec <- lambda(jumpTimes, xi)
      prod(lambda_vec)*exp(-Lambda(Tend, xi))
    }
    
    b <- model.class@b.fun  
    s <- model.class@s.fun  
    h <- model.class@h.fun  
    priorDensity <- model.class@priorDensity

    dX <- diff(X)
    dt <- diff(t)
    lt <- length(t)
    # starting values
    phi <- start$phi
    theta <- start$theta
    gamma2 <- start$gamma2
    xi <- start$xi
    
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    
    if(all(xi > 0)){
      proposals_xi <- list()
      proposals_xi$draw <- function(xi_old, propSd){
        proposal(xi_old, propSd)
      }
      proposals_xi$ratio <- function(xi_drawn, xi_old, propSd){
        proposalRatio(xi_old, xi_drawn, propSd)
      }
    } else {
      proposals_xi <- list()
      proposals_xi$draw <- function(xi_old, propSd){ 
        rnorm(length(xi_old), xi_old, propSd)
      }
      proposals_xi$ratio <- function(xi_drawn, xi_old, propSd) 1
    }
    
    
    if(is.numeric(data)){
      rangeN <- 2
      post_dN <- function(dN, i, phi, gamma2, theta, xi){
        Di <- dX[i]-b(phi,t[i],X[i])*dt[i]-h(theta,t[i],X[i])*dN
        dLambda <- Lambda(t[i+1],xi)-Lambda(t[i],xi)
        exp(-dLambda)/prod(1:max(dN,1))*exp(-Di^2/(2*s(gamma2,t[i],X[i])^2*dt[i])+dN*log(dLambda))/sqrt(2*pi*s(gamma2,t[i],X[i])^2*dt[i])
      }
      drawN <- function(phi, gamma2, theta, xi, dN_old){
        dN_new <- dN_old
        cands <- lapply(dN_old, function(n) 0:(n*rangeN+5))
        prob <- lapply(1:(lt-1), function(i) sapply(cands[[i]], post_dN, i, phi, gamma2, theta, xi))
        diFu <- lapply(prob, cumsum)
        ind <- sapply(diFu, function(vec) any(is.na(vec)) | any(is.infinite(vec)) )
        u <- numeric(length(diFu))
        if(sum(ind) < length(diFu)){
          u[!ind] <- runif(length(diFu[!ind]), 0, sapply(diFu[!ind], max))
          for(j in (1:(lt-1))[!ind]){
            dN_new[j] <- cands[[j]][which(diFu[[j]] >= u[j])[1]]
          }
        }
        dN_new
      }
      if(is.null(start$N)){
        N <- simN(t, xi, 1, start = c(t[1], 0), Lambda)$N
      }else{
        N <- start$N
      }
      N_out <- matrix(0, lt, n)
      sample.N <- TRUE
    }else{
      sample.N <- FALSE
    }
    dN <- diff(N)
    
    likeli <- function(phi, gamma2, theta, dN){
      dnorm(dX, b(phi, t[-lt], X[-lt])*dt + h(theta, t[-lt], X[-lt])*dN, s(gamma2, t[-lt], X[-lt])*sqrt(dt))
    }
    
    # storage variables
    phi_out <- rep(0, n)
    theta_out <- rep(0, n)
    gamma2_out <- rep(0, n)
    xi_out <- matrix(0, n, length(xi))
    
    if(missing(propSd)){
      propSd_phi <- start$phi
      propSd_theta <- start$theta/5
      propSd_gamma2 <- start$gamma2/5
      propSd_xi <- (abs(start$xi)+0.1)/5
    }else{
      propSd_phi <- propSd[1]
      propSd_theta <- propSd[2]
      propSd_gamma2 <- propSd[3]
      propSd_xi <- propSd[4:length(propSd)]
    }
    
    for(count in 1:n){
      if(sample.N){
        dN <- drawN(phi, gamma2, theta, xi, dN)
        N <- cumsum(c(0, dN))
        jumpTimes <- dNtoTimes(dN, t[-1])
      }
      for(count2 in 1:it.xi){
        xi_drawn <- proposals_xi$draw(xi, propSd_xi)
        ratio <- proposals_xi$ratio(xi_drawn, xi, propSd_xi)
        ratio <- ratio*Lik.N(xi_drawn, jumpTimes)/Lik.N(xi, jumpTimes)
        ratio <- ratio*prod(priorDensity$xi(xi_drawn)/priorDensity$xi(xi))
        if(is.na(ratio)) ratio <- 0
        
        if(runif(1) <= ratio){
          xi <- xi_drawn
        }
      }

      phi_drawn <- proposals$draw(phi, propSd_phi)
      ratio <- prod(likeli(phi_drawn, gamma2, theta, dN)/likeli(phi, gamma2, theta, dN))
      ratio <- ratio*priorDensity$phi(phi_drawn)/priorDensity$phi(phi)
      ratio <- ratio*proposals$ratio(phi_drawn, phi, propSd_phi)
      phi[runif(1) <= ratio] <- phi_drawn
      
      theta_drawn <- proposals$draw(theta, propSd_theta)
      ratio <- prod(likeli(phi, gamma2, theta_drawn, dN)/likeli(phi, gamma2, theta, dN))
      ratio <- ratio*priorDensity$theta(theta_drawn)/priorDensity$theta(theta)
      ratio <- ratio*proposals$ratio(theta_drawn, theta, propSd_theta)
      theta[runif(1) <= ratio] <- theta_drawn
      
      gamma2_drawn <- proposal(gamma2, propSd_gamma2)
      ratio <- prod(likeli(phi, gamma2_drawn, theta, dN)/likeli(phi, gamma2, theta, dN))
      ratio <- ratio*proposalRatio(gamma2, gamma2_drawn, propSd_gamma2)
      ratio <- ratio*priorDensity$gamma2(gamma2_drawn)/priorDensity$gamma2(gamma2)
      gamma2[runif(1) <= ratio] <- gamma2_drawn
      
      # storage
      phi_out[count] <- phi
      theta_out[count] <- theta
      gamma2_out[count] <- gamma2
      xi_out[count, ] <- xi
      
      if(sample.N){
        N_out[, count] <- N
        if(count %% 1000 == 0) message(paste(count, "iterations are calculated"))
      }
      
      if (adapt && count%%50 == 0){
        propSd_phi <- ad.propSd(phi_out[(count-50+1):count], propSd_phi, count/50)
        propSd_theta <- ad.propSd(theta_out[(count-50+1):count], propSd_theta, count/50)
        propSd_gamma2 <- ad.propSd(gamma2_out[(count-50+1):count], propSd_gamma2, count/50)
        propSd_xi <- sapply(1:length(xi), function(i){
          ad.propSd(xi_out[(count-50+1):count, i], propSd_xi[i], count/50) })
      }
    }

    result <- list(phi = phi_out, gamma2 = gamma2_out, theta = theta_out, xi = xi_out)
    
    if(nMCMC > 100){
      he <- matrix(0, 3 + ncol(result$xi), 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$phi); he[3,] <- diagnostic(result$theta)
      for(i in 4:(3+ncol(result$xi))) he[i,] <- diagnostic(result$xi[, i-3])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    

    if(is.list(data)){
      result <- new(Class = "est.jumpDiffusion", theta = result$theta, phi = result$phi, gamma2 = result$gamma2, xi = result$xi,
                    model = class.to.list(model.class), t = t, Y = data$Y, N = data$N, burnIn = burnIn, thinning = thinning)
    } else {
      result <- new(Class = "est.jumpDiffusion", theta = result$theta, phi = result$phi, gamma2 = result$gamma2, xi = result$xi,
                    N.est = N_out, model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    }

    return(result)
})


########
#' Estimation for jump diffusion process
#'
#' @description Bayesian estimation of a stochastic process
#'   \eqn{Y_t = y_0 \exp( \phi t - \gamma^2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @param model.class class of the jump diffusion model including all required information, see \code{\link{Merton-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\xi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density for xi: "normal" (default) or "lognormal"
#' @param it.xi number of iterations for MH step for \eqn{\xi} inside the Gibbs sampler
#' @references 
#' Hermann, S. and F. Ruggeri (2016). Modelling Wear Degradation in Cylinder Liners. 
#' SFB 823 discussion paper 06/16.
#' 
#' Hermann, S., K. Ickstadt and C. H. Mueller (2015). 
#' Bayesian Prediction for a Jump Diffusion Process with Application to Crack Growth in Fatigue Experiments.
#' SFB 823 discussion paper 30/15.
#' @examples
#' model <- set.to.class("Merton", parameter = list(thetaT = 0.1, phi = 0.05, gamma2 = 0.1, xi = 10))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(model, t, data, 1000)
#' plot(est)
#' \dontrun{
#' est_hidden <- estimate(model, t, data$Y, 1000)
#' plot(est_hidden)
#' }
#' @export
setMethod(f = "estimate", signature = "Merton",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal"), it.xi = 10) {
            proposal <- match.arg(proposal)
            
            if (!is.vector(t, mode = "numeric")) 
              stop(
                "t has to be a vector"
              )
            if (!is.vector(data))
              stop(
                "data has to be a vector or a list with entries N and Y"
              )
            if (is.vector(data, mode = "numeric") && length(t) != length(data))
              stop(
                "t and data must have the same length"
              )
            if (is.list(data) &&  !(length(t) == length(data$Y) && length(t) == length(data$N)) ) 
              stop(
                "t, data$N and data$Y must have the same length"
              )
            if(!missing(propSd) && length(model.class@start$xi) != length(propSd))
              stop(
                "propSd must have length of phi"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )

    if(is.list(data)){  
      X <- data$Y
      N <- data$N
      jumpTimes <- dNtoTimes(diff(N), t[-1])
    }else{
      X <- data
    } 
    Tend <- max(t)
            
    n <- nMCMC
    start <- model.class@start
    prior <- model.class@prior        
    Lambda <- model.class@Lambda
    lambda <- function(t, xi, h = 1e-05) (Lambda(t+h, xi)-Lambda(t, xi))/h
    Lik.N <- function(xi, jumpTimes){
      lambda_vec <- lambda(jumpTimes, xi)
      prod(lambda_vec)*exp(-Lambda(Tend, xi))
    }
    priorDensity <- model.class@priorDensity        

    Delta <- diff(t)
    t.l <- t
    
    dlX <- diff(log(X))
    x0 <- X[1]
    X <- X[-1]
    t <- t[-1]
    if(is.list(data)){
      dN <- diff(N)
      N <- N[-1]
    }
    logX <- log(X)
    logx0 <- log(x0)
    l <- length(t)
    if(any(t == 0)){
      help <- matrix(rep(t, l-1), nrow = l-1, ncol = l-1)
      Th <- matrix(0, l-1, l-1)
      Th[lower.tri(Th, diag = TRUE)] <- t(help)[lower.tri(Th, diag = TRUE)]
      Th[upper.tri(Th)] <- help[upper.tri(help)]
      T_1 <- solve(Th)
      T_1 <- cbind(rep(0, l), rbind(rep(0, l-1), T_1))
    }else{
      help <- matrix(rep(t, l), nrow = l, ncol = l)
      Th <- matrix(0, l, l)
      Th[lower.tri(Th, diag = TRUE)] <- t(help)[lower.tri(Th, diag = TRUE)]
      Th[upper.tri(Th)] <- help[upper.tri(help)]
      T_1 <- solve(Th)
    }
    postPhi <- function(gamma2, thetaT, N){
      Vpost <- 1/( t%*%T_1%*%t/gamma2 + 1/prior$v.phi )
      mpost <- Vpost*( prior$m.phi/prior$v.phi + 1/gamma2*t%*%T_1%*%(logX - logx0 + gamma2*t/2 - thetaT*N) )
      
      rnorm(1, mpost, sqrt(Vpost))
    }
    postThetaT <- function(gamma2, phi, N){
      Vpost <- 1/( N%*%T_1%*%N/gamma2 + 1/prior$v.thetaT )
      mpost <- Vpost*( prior$m.thetaT/prior$v.thetaT + 1/gamma2*N%*%T_1%*%(logX - logx0 - phi*t + gamma2*t/2) )
      
      rnorm(1, mpost, sqrt(Vpost))
    }
    proposalGamma <- function(gamma2, phi, thetaT, N){
      lt <- length(t)
      aPost <- prior$alpha.gamma+lt/2
      
      u_X <- logx0 + (phi-gamma2/2)*t + thetaT*N
      bPost <- prior$beta.gamma+(logX-u_X)%*%T_1%*%(logX-u_X)/2
      1/rgamma(1, aPost, bPost)
    }
    RatioGamma <- function(gamma2, gamma2_drawn, phi, thetaT, N){
      h1 <- t%*%T_1%*%t
      h2 <- (logX - logx0 - phi*t - thetaT*N)%*%T_1%*%t
      exp((1/gamma2_drawn - 1/gamma2)*( (gamma2-gamma2_drawn)*h2/2 + (gamma2-gamma2_drawn)*h1/8))
    }

    # starting values
    phi <- start$phi
    thetaT <- start$thetaT
    gamma2 <- start$gamma2
    xi <- start$xi
    if(missing(propSd)) propSd <- (abs(start$xi)+0.1)/10
    
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(xi_old, propSd){
        proposal(xi_old, propSd)
      }
      proposals$ratio <- function(xi_drawn, xi_old, propSd){
        proposalRatio(xi_old, xi_drawn, propSd)
      }
    } else {
      proposals <- list()
      proposals$draw <- function(xi_old, propSd){ 
        rnorm(length(xi_old), xi_old, propSd)
      }
      proposals$ratio <- function(xi_drawn, xi_old, propSd) 1
    }
    
    if(is.numeric(data)){
      rangeN <- 2
      
      post_dN <- function(dN, i, phi, gamma2, thetaT, xi){
        Di <- dlX[i]-(phi-gamma2/2)*Delta[i]-thetaT*dN
        dLambda <- Lambda(t.l[i], xi)-Lambda(t.l[i-1], xi)
        exp(-dLambda)/prod(1:max(dN,1))*exp(-Di^2/(2*gamma2*Delta[i])+dN*log(dLambda))/sqrt(2*pi*gamma2*Delta[i])
      }
      drawN <- function(phi, gamma2, thetaT, xi, dN_old){
        dN_new <- dN_old
        cands <- lapply(dN_old, function(n) 0:(n*rangeN+5))
        prob <- lapply(1:l, function(i) sapply(cands[[i]], post_dN, i, phi, gamma2, thetaT, xi))
        diFu <- lapply(prob, cumsum)
        ind <- sapply(diFu, function(vec) any(is.na(vec)) | any(is.infinite(vec)) )
        u <- numeric(length(diFu))
        if(sum(ind) < length(diFu)){
          u[!ind] <- runif(length(diFu[!ind]), 0, sapply(diFu[!ind], max))
          for(j in (1:l)[!ind]){
            dN_new[j] <- cands[[j]][which(diFu[[j]]>=u[j])[1]]
          }
        }
        dN_new
      }
      if(is.null(start$dN)){
        dN <- diff(simN(t.l, xi, 1, start = c(t[1], 0), Lambda = Lambda)$N)
      }else{
        dN <- start$dN
      }
      N <- cumsum(dN)
      sample.N <- TRUE
      
      N_out <- matrix(0, l, n)
      
    }else{
      sample.N <- FALSE
    }
    
    # storage variables
    phi_out <- rep(0, n)
    thetaT_out <- rep(0, n)
    gamma2_out <- rep(0, n)
    xi_out <- matrix(0, n, length(xi))
    
    for(count in 1:n){
      
      if(sample.N){
        dN <- drawN(phi, gamma2, thetaT, xi, dN)
        N <- cumsum(dN)
        jumpTimes <- dNtoTimes(dN, t)
        if(count %% 1000 == 0) message(paste(count, "iterations are calculated"))
      }
      for(count2 in 1:it.xi){
        xi_drawn <- proposals$draw(xi, propSd)
        ratio <- proposals$ratio(xi_drawn, xi, propSd)
        ratio <- ratio*Lik.N(xi_drawn, jumpTimes)/Lik.N(xi, jumpTimes)
        ratio <- ratio*prod(priorDensity(xi_drawn)/priorDensity(xi))
        if(is.na(ratio)) ratio <- 0
        
        if(runif(1) <= ratio){
          xi <- xi_drawn
        }
      }
      
      phi <- postPhi(gamma2, thetaT, N)
      
      thetaT <- postThetaT(gamma2, phi, N)
      
      gamma2_drawn <- proposalGamma(gamma2, phi, thetaT, N)
      gamma2[runif(1) <= RatioGamma(gamma2, gamma2_drawn, phi, thetaT, N)] <- gamma2_drawn
      
      # storage
      phi_out[count] <- phi
      thetaT_out[count] <- thetaT
      gamma2_out[count] <- gamma2
      xi_out[count, ] <- xi
      if(sample.N){
        N_out[, count] <- N
      }
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:length(xi), function(i){
          ad.propSd(xi_out[(count-50+1):count, i], propSd[i], count/50) })
      }
    }

    result <- list(phi = phi_out, gamma2 = gamma2_out, thetaT = thetaT_out, xi = xi_out)
    
    if(nMCMC > 100){
      he <- matrix(0, 3 + ncol(result$xi), 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$phi); he[3,] <- diagnostic(result$thetaT)
      for(i in 4:(3+ncol(result$xi))) he[i,] <- diagnostic(result$xi[, i-3])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    
            
            
    if(is.list(data)){
      result <- new(Class = "est.Merton", thetaT = result$thetaT, phi = result$phi, gamma2 = result$gamma2, xi = result$xi,                    
                    model = class.to.list(model.class), t = t.l, Y = data$Y, N = data$N, burnIn = burnIn, thinning = thinning)
    } else {
      result <- new(Class = "est.Merton", thetaT = result$thetaT, phi = result$phi, gamma2 = result$gamma2,
                    xi = result$xi, N.est = N_out,
                    model = class.to.list(model.class), t = t.l, Y = data, burnIn = burnIn, thinning = thinning)
    }
     return(result)

})


########
#' Estimation for regression model dependent on Poisson process
#'
#' @description Bayesian estimation of the parameter of the regression model
#'   \eqn{y_i = f(t_i, N_{t_i}, \theta) + \epsilon_i} with
#'   \eqn{N_t\sim Pois(\Lambda(t, \xi)), \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t))}.
#' @param model.class class of the regression model based on the NHPP including all required information, see \code{\link{jumpRegression-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{(\theta, \xi)}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density for \eqn{\theta}: "normal" (default) or "lognormal" (for positive parameters)
#' @param it.xi number of iterations for MH step for \eqn{\xi} inside the Gibbs sampler
#' @section Proposal densities:
#' For \eqn{\theta}, there is the possibility to choose "normal" or "lognormal". 
#' The proposal density for \eqn{\xi} depends on the starting value of \eqn{\xi}. If all components are positive, the proposal density is lognormal, and normal otherwise.
#'
#' @references 
#' Heeke, G., S. Hermann, R. Maurer, K. Ickstadt, and C. H. Mueller (2015). 
#' Stochastic Modeling and Statistical Analysis of Fatigue Tests on Prestressed Concrete Beams under Cyclic Loadings.
#' SFB 823 discussion paper 25/15.
#' @examples
#' t <- seq(0,1, by = 0.01)
#' model <- set.to.class("jumpRegression", fun = function(t, N, theta) exp(theta[1]*t) + theta[2]*N, 
#'                    parameter = list(theta = c(2, 2), gamma2 = 0.25, xi = c(3, 0.5)),
#'                    Lambda = function(t, xi) (t/xi[2])^xi[1])
#' data <- simulate(model, t = t, plot.series = FALSE)
#' est <- estimate(model, t, data, 1000) 
#' plot(est)
#' \dontrun{
#' # work in progress
#' est_hid <- estimate(model, t, data$Y, 1000)
#' plot(est_hid)
#' }
#' @export
setMethod(f = "estimate", signature = "jumpRegression",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal"), it.xi = 10) {

    proposal <- match.arg(proposal)

    if (!is.vector(t, mode = "numeric")) 
      stop(
        "t has to be a vector"
      )
    if (!is.vector(data))
      stop(
        "data has to be a vector or a list with entries N and Y"
      )
    if (is.vector(data, mode = "numeric") && length(t) != length(data))
      stop(
        "t and data must have the same length"
      )
    if (is.list(data) &&  !(length(t) == length(data$Y) && length(t) == length(data$N)) ) 
      stop(
        "t, data$N and data$Y must have the same length"
      )
    if(!missing(propSd) && length(model.class@start$theta) != length(propSd))
      stop(
        "propSd must have length of phi"
      )
    if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
      stop(
        "nMCMC has to be a natural number"
      )
    if(any(model.class@start$theta < 0) && proposal == "lognormal")
      stop(
        "lognormal proposal density has positive support and starting value is negative"
      )

    
    if(is.list(data)){
      Y <- data$Y
      N <- data$N
      jumpTimes <- dNtoTimes(diff(N), t[-1])
    }else{
      Y <- data
    }
    Tend <- max(t)
            
    fun <- model.class@fun
    n <- nMCMC
    start <- model.class@start
    prior <- model.class@prior
    Lambda <- model.class@Lambda
    lambda <- function(t, xi, h = 1e-05) (Lambda(t+h, xi)-Lambda(t, xi))/h
    Lik.N <- function(xi, jumpTimes){
      lambda_vec <- lambda(jumpTimes, xi)
      prod(lambda_vec)*exp(-Lambda(Tend, xi))
    }
    sVar <- model.class@sT.fun
    
    lt <- length(t)
    # starting values
    theta <- start$theta
    gamma2 <- start$gamma2
    xi <- start$xi
    
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    if(all(xi > 0)){
      proposals_xi <- list()
      proposals_xi$draw <- function(xi_old, propSd){
        proposal(xi_old, propSd)
      }
      proposals_xi$ratio <- function(xi_drawn, xi_old, propSd){
        proposalRatio(xi_old, xi_drawn, propSd)
      }
    } else {
      proposals <- list()
      proposals_xi$draw <- function(xi_old, propSd){ 
        rnorm(length(xi_old), xi_old, propSd)
      }
      proposals_xi$ratio <- function(xi_drawn, xi_old, propSd) 1
    }
    
    if(is.numeric(data)){
      
      Npart <- 10
      rangeN <- 2
      A.to.B = function(A){
        B <- matrix(0, Npart, lt)
        B[,lt] <- 1:Npart;
        for (t in seq(lt, 2, by = -1)){
          B[,t-1] <- A[B[,t], t-1]}
        return(B)
      }
      
      CSMC = function(theta, gamma2, xi, N.cond, B.fixed, conditional = TRUE){# conditional SMC
        # N.cond = the old samples
        x <- matrix(0, Npart, lt)
        w <- matrix(1, Npart, lt)
        W <- matrix(1, Npart, lt)
        parents <-  matrix(1, Npart, lt-1)
        
        # initialisation
        x[,1] <- 0  # poisson always starts in 0
        if(conditional){
          x[B.fixed[1], 1] <- N.cond[1]
        }else{
          N.cond <- numeric(lt)
        }
        w[,1] <- dnorm(Y[1], mean = fun(t[1], x[,1], theta), sd = sqrt(gamma2))
        W[,1] <- w[,1]/sum(w[,1])
        
        for (n in 2:lt){
          if(conditional){
            set.parents  <- (1:Npart)[-B.fixed[n]]
            
            On_1 <- rmultinom(1, Npart-1, W[,n-1])
            O <- On_1[B.fixed[n-1]] + 1
            he <- sample(set.parents, O-1)
            
            parents[B.fixed[n], n-1] <- B.fixed[n-1]
            parents[he, n-1] <- B.fixed[n-1]
            parents[-c(B.fixed[n],he), n-1] <- sample(set.parents, Npart-O, replace = TRUE, prob =  W[-B.fixed[n], n-1])
            
          }else{
            set.parents <- 1:Npart
            parents[, n-1] <- sample(1:Npart, Npart, replace = TRUE, prob = W[,n-1])
          }
          
          x[,1:(n-1)] <- x[parents[,n-1], 1:(n-1)]
          x.past <- x[,n-1]
          
          dr <- function(Ni, dNi_old){
            cands <- 0:(dNi_old*rangeN + 5)
            prob <- dnorm(Y[n], fun(t[n], Ni+cands, theta), sqrt(2*gamma2))*
              dpois(Ni+cands, Lambda(t[n], xi))
            diFu <- cumsum(prob)
            u <- runif(1, 0, max(diFu))
            Ni + cands[which(diFu >= u)[1]]
          }
          x.new <- sapply(x.past, dr, diff(N.cond)[n-1])
          
          x[set.parents, n] <- x.new[set.parents]
          x[-set.parents, n] <- N.cond[n]
          wn <- dpois(x[, n], Lambda(t[n], xi))*
            dnorm(Y[n], fun(t[n], x[, n], theta), sqrt(gamma2))/
            (dnorm(Y[n], fun(t[n], x[, n], theta), sqrt(2*gamma2))*
               dpois(x[, n], Lambda(t[n], xi)))
          if(sum(wn) == 0) wn <- rep(1, length(wn))
          w[, n] <- wn
          W[, n] <- w[, n]/sum(w[, n])
        }
        lign.B <- A.to.B(parents)
        indice <- sample(1:Npart, 1, prob = W[, lt])
        X <- x[indice, ]
        B.fixed <- lign.B[indice,]
        return(list(N = X, B.fixed = B.fixed) )
      }
      
      if(is.null(start$N)){
        #      N <- simN(t, xi, 1, start = c(t[1], 0), Lambda = Lambda)$N
        he <- CSMC(theta, gamma2, xi, conditional = FALSE)
        N <- he$N
        B.fixed <- he$B.fixed
      }else{
        N <- start$N
      }
      N_out <- matrix(0, lt, n)
      sample.N <- TRUE
    }else{
      sample.N <- FALSE
    }
    ltheta <- length(start$theta)

    if(missing(propSd)){
      propSd_theta <- abs(prior$m.theta)/20
      propSd_xi <- (abs(start$xi)+0.1)/10
    }else{
      propSd_theta <- propSd[1:ltheta]
      propSd_xi <- propSd[(ltheta+1):length(propSd)]
    }
    
    
    postTheta <- function(N, gamma2, lastPhi, propSd){  
      phi_old <- lastPhi
      phi_drawn <- proposals$draw(phi_old, propSd)
      ratio <- prod(dnorm(phi_drawn, prior$m.theta, sqrt(prior$v.theta)))/prod(dnorm(phi_old, prior$m.theta, sqrt(prior$v.theta)))
      ratio <- ratio* prod(dnorm(Y, fun(t, N, phi_drawn), sqrt(gamma2*sVar(t)))/dnorm(Y, fun(t, N, phi_old), sqrt(gamma2*sVar(t))))
      ratio <- ratio * proposals$ratio(phi_drawn, phi_old, propSd)
      if(is.na(ratio)) ratio <- 0
      if(runif(1) < ratio){
        phi_old <- phi_drawn
      }
      phi_old
    }
    
    postGamma2 <- function(theta, N){
      alphaPost <- prior$alpha.gamma + lt/2
      betaPost <-  prior$beta.gamma + sum((Y-fun(t, N, theta))^2/sVar(t))/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    theta_out <- matrix(0, n, ltheta)
    gamma2_out <- numeric(n)
    xi_out <- matrix(0, n, length(xi))
    
    for(count in 1:n){
      if(sample.N){
        he <- CSMC(theta, gamma2, xi, N, B.fixed)
        N <- he$N
        jumpTimes <- dNtoTimes(diff(N), t[-1])
        B.fixed <- he$B.fixed
      }
      for(count2 in 1:it.xi){
        xi_drawn <- proposals_xi$draw(xi, propSd_xi)
        ratio <- proposals_xi$ratio(xi_drawn, xi, propSd_xi)
        ratio <- ratio*Lik.N(xi_drawn, jumpTimes)/Lik.N(xi, jumpTimes)
        if(is.na(ratio)) ratio <- 0
        
        if(runif(1) <= ratio){
          xi <- xi_drawn
        }
      }
      
      theta <- postTheta(N, gamma2, theta, propSd_theta)
      
      gamma2 <- postGamma2(theta, N)
      
      # storage
      theta_out[count, ] <- theta
      gamma2_out[count] <- gamma2
      xi_out[count, ] <- xi
      if(sample.N){
        N_out[, count] <- N
      }
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:ltheta, function(i){
          ad.propSd(theta_out[(count-50+1):count, i], propSd_theta[i], count/50) })
        propSd_xi <- sapply(1:length(xi), function(i){
          ad.propSd(xi_out[(count-50+1):count, i], propSd_xi[i], count/50) })
      }
      
    }

    
    result <- list(gamma2 = gamma2_out, theta = theta_out, xi = xi_out)
    
    if(nMCMC > 100){
      he <- matrix(0, ncol(result$theta) + ncol(result$xi) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$theta)+1)) he[i, ] <- diagnostic(result$theta[,i-1])
      for(i in (ncol(result$theta)+2):(ncol(result$theta) + ncol(result$xi) + 1)) he[i, ] <- diagnostic(result$xi[,i - ncol(result$theta) - 1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    
    

    
    if(is.list(data)){
      result <- new(Class = "est.jumpRegression", theta = result$theta, gamma2 = result$gamma2, xi = result$xi,
                    model = class.to.list(model.class), t = t, Y = data$Y, N = data$N, burnIn = burnIn, thinning = thinning)

    }else{
      result <- new(Class = "est.jumpRegression", theta = result$theta, gamma2 = result$gamma2,
                    xi = result$xi, N.est = N_out,
                    model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    }
    return(result)
})


########
#' Estimation for regression model
#'
#' @description Bayesian estimation of the parameter of the regression model
#'   \eqn{y_i = f(\phi, t_i) + \epsilon_i, \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t_i))}.
#' @param model.class class of the regression model including all required information, see \code{\link{Regression-class}}
#' @param t vector of time points
#' @param data vector of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\phi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#'
#' @references 
#' Hermann, S., K. Ickstadt, and C. H. Mueller (2016). 
#' Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. 
#' Applied Stochastic Models in Business and Industry, DOI: 10.1002/asmb.2175.
#' @examples
#' t <- seq(0,1, by = 0.01)
#' model <- set.to.class("Regression", fun = function(phi, t) phi[1]*t + phi[2], 
#'                    parameter = list(phi = c(1,2), gamma2 = 0.1))
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est <- estimate(model, t, data, 1000)
#' plot(est)
#' @export
setMethod(f = "estimate", signature = "Regression",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal")) {
            proposal <- match.arg(proposal)

            if (!is.vector(t, mode = "numeric")) 
              stop(
                "t has to be a vector"
              )
            if (!is.vector(data, mode = "numeric"))
              stop(
                "data has to be a vector"
              )
            if (length(t) != length(data))
              stop(
                "t and data must have the same length"
              )
            if(!missing(propSd) && length(model.class@start$phi) != length(propSd))
              stop(
                "propSd must have length of phi"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start$phi < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )
            

    y <- data
    prior <- model.class@prior
    start <- model.class@start
    fODE <- model.class@fun
    sVar <- model.class@sT.fun
    len <- nMCMC
    
    
    if(missing(propSd)) propSd <- abs(prior$m.phi)/50
    lt <- length(t)
    lphi <- length(start$phi)
    
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    postPhi <- function(lastPhi, gamma2, propSd){
      phi_old <- lastPhi
      phi_drawn <- proposals$draw(phi_old, propSd)
      ratio <- prod(dnorm(phi_drawn, prior$m.phi, sqrt(prior$v.phi)) / dnorm(phi_old, prior$m.phi, sqrt(prior$v.phi)))
      ratio <- ratio* prod(dnorm(y, fODE(phi_drawn, t), sqrt(gamma2*sVar(t)))/dnorm(y, fODE(phi_old, t), sqrt(gamma2*sVar(t))))
      ratio <- ratio * proposals$ratio(phi_drawn, phi_old, propSd)
      if(is.na(ratio)) ratio <- 0
      if(runif(1) < ratio){
        phi_old <- phi_drawn
      }
      phi_old
    }
    
    postGamma2 <- function(phi){
      alphaPost <- prior$alpha.gamma + lt/2
      betaPost <-  prior$beta.gamma + sum((y-fODE(phi, t))^2/sVar(t))/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    phi_out <- matrix(0, len, length(prior$m.phi))
    gamma2_out <- numeric(len)
    
    phi <- start$phi
    gamma2 <- start$gamma2
    
    for(count in 1:len){
      
      phi <- postPhi(phi, gamma2, propSd)
      gamma2 <- postGamma2(phi)
      
      phi_out[count,] <- phi
      gamma2_out[count] <- gamma2
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:length(phi), function(i){
          ad.propSd(phi_out[(count-50+1):count, i], propSd[i], count/50) })
      }
      
      
    }
    result <- list(phi = phi_out, gamma2 = gamma2_out)
    
    if(nMCMC > 100){
      he <- matrix(0, ncol(result$phi) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$phi)+1)) he[i, ] <- diagnostic(result$phi[,i-1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    

    result <- new(Class = "est.Regression", phi = result$phi, gamma2 = result$gamma2,
                  model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    return(result)

 })


########
#' Estimation for the hierarchical (mixed) regression model
#'
#' @description Bayesian estimation of the parameter of the hierarchical regression model
#'   \eqn{y_{ij} = f(\phi_j, t_{ij}) + \epsilon_{ij}, \phi_j\sim N(\mu, \Omega),
#'   \epsilon_{ij}\sim N(0,\gamma^2\widetilde{s}(t_{ij}))}.
#' @param model.class class of the hierarchical regression model including all required information, see \code{\link{mixedRegression-class}}
#' @param t list or vector of time points
#' @param data list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @param propSd vector of proposal variances for \eqn{\phi}
#' @param adapt if TRUE (default), proposal variance is adapted
#' @param proposal proposal density: "normal" (default) or "lognormal" (for positive parameters)
#' @references 
#' Hermann, S., K. Ickstadt, and C. H. Mueller (2016). 
#' Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. 
#' Applied Stochastic Models in Business and Industry, DOI: 10.1002/asmb.2175.
#' @examples
#' mu <- c(10, 5); Omega <- c(0.9, 0.01)
#' phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' model <- set.to.class("mixedRegression", 
#'                  parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'                  fun = function(phi, t) phi[1]*t + phi[2], sT.fun = function(t) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = FALSE)
#' est <- estimate(model, t, data[1:20,], 1000)
#' plot(est)
#'
#' @export
setMethod(f = "estimate", signature = "mixedRegression",
          definition = function(model.class, t, data, nMCMC, propSd, adapt = TRUE, proposal = c("normal", "lognormal")) {
            proposal <- match.arg(proposal)

            if (!(is.vector(t))) 
              stop(
                "t has to be a vector or a list"
              )
            if (!(is.matrix(data) || is.list(data)))
              stop(
                "data has to be a matrix or a list"
              )
            if (is.vector(t, mode = "numeric") && (length(t) != nrow(data) && length(t) != ncol(data)))
              stop(
                "length of t has to be equal to the columns/rows of data"
              )
            if (is.list(t) && ( length(t) != length(data) || !all(sapply(t, length) == sapply(data, length))))
              stop(
                "data must match with t"
              )
            if(!missing(propSd) && ncol(model.class@start$phi) != length(propSd))
              stop(
                "propSd must have length of phi_j"
              )
            if(!is.numeric(nMCMC) || length(nMCMC) > 1 || nMCMC < 1)
              stop(
                "nMCMC has to be a natural number"
              )
            if(any(model.class@start$phi < 0) && proposal == "lognormal")
              stop(
                "lognormal proposal density has positive support and starting value is negative"
              )

    prior <- model.class@prior
    start <- model.class@start
    fODE <- model.class@fun
    sVar <- model.class@sT.fun
    len <- nMCMC
    if(missing(propSd)) propSd <- abs(start$mu)/5
    
    if(is.matrix(data)){
      if(nrow(data) == length(t)){
        y <- t(data)
      }else{
        if(ncol(data) != length(t)){
          stop("length of t has to be equal to the columns of y")
        }
      }
      times <- list()
      y <- list()
      for(i in 1:nrow(data)){
        times[[i]] <- t
        y[[i]] <- data[i,]
      }
    }else{
      y <- data
      times <- t
    }
    
    postOm <- function(phi,mu){
      postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)
    }
    if(proposal == "lognormal"){
      proposals <- list()
      proposals$draw <- function(old, propSd){
        proposal(old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd){
        proposalRatio(old, drawn, propSd)
      }
    }else{
      proposals <- list()
      proposals$draw <- function(old, propSd){ 
        rnorm(length(old), old, propSd)
      }
      proposals$ratio <- function(drawn, old, propSd) 1
    }
    postPhii <- function(lastPhi, mu, Omega, gamma2, y, t, propSd){
      lt <- length(t)
      phi_old <- lastPhi
      
      phi_drawn <- proposals$draw(phi_old, propSd)
      ratio <- prod(dnorm(phi_drawn, mu, sqrt(Omega))/dnorm(phi_old, mu, sqrt(Omega)))
      ratio <- ratio* prod(dnorm(y, fODE(phi_drawn,t), sqrt(gamma2*sVar(t)))/dnorm(y, fODE(phi_old,t), sqrt(gamma2*sVar(t))))
      ratio <- ratio * proposals$ratio(phi_drawn, phi_old, propSd)
      if(is.na(ratio)){ratio <- 0}
      if(runif(1) < ratio){
        phi_old <- phi_drawn
      }
      phi_old
    }
    N_all <- sum(sapply(y,length))
    n <- length(y)
    postGamma2 <- function(phi){
      alphaPost <- prior$alpha.gamma + N_all/2
      help <- numeric(n)
      for(i in 1:n){
        help[i] <- sum((y[[i]]-fODE(phi[i,], times[[i]]))^2/sVar(times[[i]]))
      }
      betaPost <-  prior$beta.gamma + sum(help)/2
      1/rgamma(1, alphaPost, betaPost)
    }
    
    phi_out <- list()
    mu_out <- matrix(0, len, length(start$mu))
    Omega_out <- matrix(0, len, length(start$mu))
    gamma2_out <- numeric(len)
    
    phi <- start$phi
    gamma2 <- start$gamma2
    mu <- start$mu
    Omega <- postOm(phi, mu)
    
    for(count in 1:len){
      
      for(i in 1:n){
        phi[i,] <- postPhii(phi[i,], mu, Omega, gamma2, y[[i]], times[[i]], propSd)
      }
      mu <- postmu(phi, prior$m.mu, prior$v.mu, Omega)
      Omega <- postOm(phi, mu)
      gamma2 <- postGamma2(phi)
      
      phi_out[[count]] <- phi
      mu_out[count,] <- mu
      Omega_out[count,] <- Omega
      gamma2_out[count] <- gamma2
      
      if (adapt && count%%50 == 0){
        propSd <- sapply(1:length(phi[1,]), function(i){
          ad.propSd(sapply(phi_out[(count-50+1):count], function(mat) mat[1,i]), propSd[i], count/50) })
      }
      
    }
    result <- list(phi = phi_out, mu = mu_out, Omega = Omega_out, gamma2 = gamma2_out)

    if(nMCMC > 100){
      he <- matrix(0, ncol(result$mu) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$mu)+1)) he[i, ] <- diagnostic(result$mu[,i-1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    } else {
      burnIn <- 0; thinning <- 1
    }
    

    if(is.matrix(data)){
      result <- new(Class = "est.mixedRegression", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                    model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    }else{
      result <- new(Class = "est.mixedRegression", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                    model = class.to.list(model.class), t = t[[1]], Y = matrix(data[[1]], 1),
                    t.list = t, Y.list = data, burnIn = burnIn, thinning = thinning)
    }
    
    
    return(result)

})

