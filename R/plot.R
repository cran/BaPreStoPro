#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the jump diffusion model.
#' @param x est.jumpDiffusion class, created with method \code{\link{estimate,jumpDiffusion-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \theta, \gamma^2, \xi, N)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("jumpDiffusion", Lambda = function(t, xi) (t/xi[2])^xi[1],
#' parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)))
#' data <- simulate(model, t = seq(0, 1, by = 0.01), y0 = 0.5, plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data, 1000)  # nMCMC small for example
#' plot(est)
#' plot(est, burnIn = 100, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2, 3)), xlab = "iteration")
#' # plot only for phi and xi ...
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, FALSE, FALSE, TRUE, TRUE))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), par2plot = c(TRUE, rep(FALSE, 4)), main = "")
#' @export
#' @import stats
#' @import graphics
#' @import methods
#' 
setMethod(f = "plot", signature = "est.jumpDiffusion", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...) {
  old.settings <- par(no.readonly = TRUE)
  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, length(x@phi), by = thinning)
  }else{
    ind <- 1:length(x@phi)
  }

  p <- ncol(x@xi)
  if(style == "chains"){
    p1 <- c(1, 0)[c(sum(dim(x@N.est)) > 0, sum(dim(x@N.est)) == 0)]
  } else {
    p1 <- 0
  }
  if(missing(par2plot)) par2plot <- rep(TRUE, 3+p+p1)
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
  }
  
  par(par.options)
  
  if(style == "chains"){
      
    if(par2plot[1]){
      plot(x@phi[ind], type = "l", ylab = expression(phi), ...)
      if(priorMeans) abline(h = x@model$phi, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[2]){
      plot(x@theta[ind], type = "l", ylab = expression(theta), ...)
      if(priorMeans)  abline(h = x@model$theta, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[3]){
      plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
      if(priorMeans)  abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(any(par2plot[3 + seq_len(p)])){
      for(i in 1:p){
        if(par2plot[3 + i]){
          plot(x@xi[ind, i], type = "l", ylab = bquote(xi[.(i)]), ...)
          if(priorMeans)  abline(h = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
    }
    
    if(sum(dim(x@N.est) > 0)){
      if(par2plot[3 + p + 1]){
        plot(x@t, x@N.est[,length(x@phi)], type = "l", xlab = "t", ylab = "N", ylim = range(x@N.est[,ind]), ...)
        for(i in ind) lines(x@t, x@N.est[,i])
      }  
    }
  }
  if(style == "acf"){

    if(par2plot[1]) acf(x@phi[ind], xlab = expression(phi), ...)
    if(par2plot[2]) acf(x@theta[ind], xlab = expression(theta), ...)
    if(par2plot[3]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
    for(i in 1:p){
      if(par2plot[3 + i]) acf(x@xi[ind, i], xlab = bquote(xi[.(i)]), ...)
    }
  }
  if(style == "density"){

    if(par2plot[1]){
      plot(density(x@phi[ind]), xlab = expression(phi), ...)
      if(priorMeans) abline(v = x@model$phi, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[2]){
      plot(density(x@theta[ind]), xlab = expression(theta), ...)
      if(priorMeans)  abline(v = x@model$theta, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[3]){
      plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
      if(priorMeans)  abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    for(i in 1:p){
      if(par2plot[3 + i]){
        plot(density(x@xi[ind, i]), xlab = bquote(xi[.(i)]), ...)
        if(priorMeans)  abline(v = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
  }

  par(old.settings[which(names(old.settings) %in% names(par.options))])
})


#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the Merton model.
#' @param x est.Merton class, created with method \code{\link{estimate,Merton-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \widetilde{\theta}, \gamma^2, \xi, N)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("Merton", Lambda = function(t, xi) (t/xi[2])^xi[1],
#' parameter = list(thetaT = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)))
#' data <- simulate(model, t = seq(0, 1, by = 0.01), y0 = 0.5, plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data, 1000)  # nMCMC small for example
#' plot(est)
#' plot(est, burnIn = 100, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2, 3)), xlab = "iteration")
#' # plot only for phi and xi ...
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, FALSE, FALSE, TRUE, TRUE))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), par2plot = c(TRUE, rep(FALSE, 4)), main = "")
#' @export
setMethod(f = "plot", signature = "est.Merton", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...) {
    old.settings <- par(no.readonly = TRUE)
    style <- match.arg(style) 
    if(reduced){
      if(missing(thinning)) thinning <- x@thinning
      if(missing(burnIn)) burnIn <- x@burnIn
      ind <- seq(burnIn + 1, length(x@phi), by = thinning)
    }else{
      ind <- 1:length(x@phi)
    }
    
    p <- ncol(x@xi)
    if(style == "chains"){
      p1 <- c(1, 0)[c(sum(dim(x@N.est)) > 0, sum(dim(x@N.est)) == 0)]
    } else {
      p1 <- 0
    }
    if(missing(par2plot)) par2plot <- rep(TRUE, 3+p+p1)
    if(missing(par.options)){
      if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
      if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
    }
    
    par(par.options)
    
    if(style == "chains"){
      
      if(par2plot[1]){
        plot(x@phi[ind], type = "l", ylab = expression(phi), ...)
        if(priorMeans) abline(h = x@model$phi, col = col.priorMean, lty = lty.priorMean)
      }
      if(par2plot[2]){
        plot(x@thetaT[ind], type = "l", ylab = expression(tilde(theta)), ...)
        if(priorMeans)  abline(h = x@model$thetaT, col = col.priorMean, lty = lty.priorMean)
      }
      if(par2plot[3]){
        plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
        if(priorMeans)  abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
      }
      if(any(par2plot[3 + seq_len(p)])){
        for(i in 1:p){
          if(par2plot[3 + i]){
            plot(x@xi[ind, i], type = "l", ylab = bquote(xi[.(i)]), ...)
            if(priorMeans)  abline(h = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
          }
        }
      }
      
      if(sum(dim(x@N.est) > 0)){
        if(par2plot[3 + p + 1]){
          plot(x@t, x@N.est[,length(x@phi)], type = "l", xlab = "t", ylab = "N", ylim = range(x@N.est[,ind]), ...)
          for(i in ind) lines(x@t, x@N.est[,i])
        }  
      }
    }
    if(style == "acf"){
      
      if(par2plot[1]) acf(x@phi[ind], xlab = expression(phi), ...)
      if(par2plot[2]) acf(x@thetaT[ind], xlab = expression(tilde(theta)), ...)
      if(par2plot[3]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
      for(i in 1:p){
        if(par2plot[3 + i]) acf(x@xi[ind, i], xlab = bquote(xi[.(i)]), ...)
      }
    }
    if(style == "density"){
      
      if(par2plot[1]){
        plot(density(x@phi[ind]), xlab = expression(phi), ...)
        if(priorMeans) abline(v = x@model$phi, col = col.priorMean, lty = lty.priorMean)
      }
      if(par2plot[2]){
        plot(density(x@thetaT[ind]), xlab = expression(tilde(theta)), ...)
        if(priorMeans)  abline(v = x@model$theta, col = col.priorMean, lty = lty.priorMean)
      }
      if(par2plot[3]){
        plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
        if(priorMeans)  abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
      }
      for(i in 1:p){
        if(par2plot[3 + i]){
          plot(density(x@xi[ind, i]), xlab = bquote(xi[.(i)]), ...)
          if(priorMeans)  abline(v = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
    }
    
    par(old.settings[which(names(old.settings) %in% names(par.options))])
})


#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the diffusion model.
#' @param x est.Diffusion class, created with method \code{\link{estimate,Diffusion-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \gamma^2)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("Diffusion", b.fun = function(phi, t, y) phi[1]-phi[2]*y,
#'     parameter = list(phi = c(10, 1), gamma2 = 0.1))
#' data <- simulate(model, t = seq(0, 1, by = 0.01), y0 = 0.5, plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data, 1000)  # nMCMC small for example
#' plot(est)
#' plot(est, burnIn = 100, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(3,1)), xlab = "iteration")
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, TRUE, FALSE))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), main = "", par2plot = c(FALSE, FALSE, TRUE))
#' @export
setMethod(f = "plot", signature = "est.Diffusion", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1,...) {
  old.settings <- par(no.readonly = TRUE)
  
  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
  }else{
    ind <- 1:length(x@gamma2)
  }
  
  p <- ncol(x@phi)

  if(missing(par2plot)) par2plot <- rep(TRUE, 1+p)
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
  }

  par(par.options)
  
  if(style == "chains"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(x@phi[ind, i], type = "l", ylab = bquote(phi[.(i)]), ...)
        if(priorMeans) abline(h = x@model$phi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[p+1]){
      plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
      if(priorMeans) abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
  }
  if(style == "acf"){
    for(i in 1:p) {
      if(par2plot[i]){
        acf(x@phi[ind, i], xlab = bquote(phi[.(i)]), ...)
      }
    }
    if(par2plot[p+1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
  }
  if(style == "density"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(density(x@phi[ind, i]), xlab = bquote(phi[.(i)]), ...)
        if(priorMeans) abline(v = x@model$phi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[p+1]){
      plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
      if(priorMeans) abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
  }

  par(old.settings[which(names(old.settings) %in% names(par.options))])
})

#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the hierarchical (mixed) diffusion model.
#' @param x est.mixedDiffusion class, created with method \code{\link{estimate,mixedDiffusion-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density", "int.phi"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\mu, \Omega, \gamma^2)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param level level for style = "int.phi"
#' @param phi in the case of simulation study: known values for phi
#' @param ... optional plot parameters
#' @examples 
#' mu <- c(10, 3, 1); Omega = c(1, 0.4, 0.01)
#' phi <- sapply(1:3, function(i) rnorm(20, mu[i], sqrt(Omega[i])))
#' model <- set.to.class("mixedDiffusion", b.fun = function(phi, t, y) phi[1]-phi[2]*y,
#'     parameter = list(mu = mu, Omega = Omega, phi = phi, gamma2 = 0.1),
#'     y0 = function(phi, t) phi[3], sT.fun = function(t, x) sqrt(abs(x)))
#' data <- simulate(model, t = seq(0, 1, by = 0.02), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.02), data, 100)  # nMCMC small for example
#' plot(est, burnIn = 10, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2,1)), xlab = "iteration")
#' plot(est, style = "acf", main = "")
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), main = "", par2plot = c(rep(FALSE, 6), TRUE))
#' plot(est, style = "int.phi", phi = phi, par2plot = c(TRUE, FALSE, FALSE))
#' @export
setMethod(f = "plot", signature = "est.mixedDiffusion", 
          definition = function(x, par.options, style = c("chains", "acf", "density", "int.phi"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, level = 0.05, phi, ...) {
  old.settings <- par(no.readonly = TRUE)

  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
  }else{
    ind <- 1:length(x@gamma2)
  }
  
  p <- ncol(x@mu)
  k <- nrow(x@phi[[1]])
  
  if(missing(par2plot)) par2plot <- rep(TRUE, 1+2*p)
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
    if(style == "int.phi") par.options <- list(mfrow = c(sum(par2plot[1:p]), 1))
  }

  par(par.options)
  
  if(style == "chains"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(x@mu[ind, i], type = "l", ylab = bquote(mu[.(i)]), ...)
        if(priorMeans) abline(h = x@model$mu[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    for(i in 1:p) {
      if(par2plot[p+i]){
        plot(x@Omega[ind, i], type = "l", ylab = bquote(omega[.(i)]), ...)
        if(priorMeans) abline(h = x@model$Omega[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[2*p+1]){
      plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
      if(priorMeans) abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
  }
  if(style == "acf"){
    for(i in 1:p) {
      if(par2plot[i]) acf(x@mu[ind, i], xlab = bquote(mu[.(i)]), ...)
    }
    for(i in 1:p){
      if(par2plot[p+i]) acf(x@Omega[ind, i], xlab = bquote(omega[.(i)]), ...)
    }
    if(par2plot[2*p+1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
  }
  if(style == "density"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(density(x@mu[ind, i]), xlab = bquote(mu[.(i)]), ...)
        if(priorMeans) abline(v = x@model$mu[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    for(i in 1:p) {
      if(par2plot[p+i]){
        plot(density(x@Omega[ind, i]), xlab = bquote(omega[.(i)]), ...)
        if(priorMeans) abline(v = x@model$Omega[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[2*p+1]){
      plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
      if(priorMeans) abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
  }
  if(style == "int.phi"){
    
    qu <- lapply(1:p, function(i) {
      apply( sapply(x@phi[ind], function(m) m[, i]), 1, function(ve) quantile(ve, c(level/2, 0.5, 1-level/2)) )
    })
    for(i in 1:p){
      if(par2plot[i]){
        plot(qu[[i]][2, ], pch = 20, ylim = range(qu[[i]]), ylab = bquote(phi[.(i)]), ...)
        segments(1:k, qu[[i]][1, ], 1:k, qu[[i]][3, ], ...)
        segments(1:k-0.1, qu[[i]][1, ], 1:k+0.1, qu[[i]][1, ], ...)
        segments(1:k-0.1, qu[[i]][3, ], 1:k+0.1, qu[[i]][3, ], ...)
        if(!missing(phi)) points(phi[, i], pch = 20, col = col.priorMean)
      }
    }
    par(old.settings[which(names(old.settings) %in% names(par.options))])
    return(invisible(qu))
  }
  
  par(old.settings[which(names(old.settings) %in% names(par.options))])
})

#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the hidden diffusion model.
#' @param x est.hiddenDiffusion class, created with method \code{\link{estimate,hiddenDiffusion-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \gamma^2, \sigma^2, Y)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("hiddenDiffusion", b.fun = function(phi, t, y) phi[1]-phi[2]*y,
#'     parameter = list(phi = c(10, 1), gamma2 = 1, sigma2 = 0.1),
#'     y0 = function(phi, t) 0.5)
#' data <- simulate(model, t = seq(0, 1, by = 0.01), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data$Y, 100)  # nMCMC small for example
#' plot(est)
#' plot(est, par2plot = c(rep(FALSE, 3), TRUE, FALSE), ylim = c(0.001, 0.1), par.options = list())
#' plot(est, burnIn = 10, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(3,1)), xlab = "iteration")
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, TRUE, FALSE, FALSE))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), main = "", par2plot = c(FALSE, FALSE, TRUE, TRUE))
#' @export
setMethod(f = "plot", signature = "est.hiddenDiffusion", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...) {
  old.settings <- par(no.readonly = TRUE)
  
  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
  }else{
    ind <- 1:length(x@gamma2)
  }
  
  p <- ncol(x@phi)
  
  if(missing(par2plot)){
    if(style == "chains")  par2plot <- rep(TRUE, 3+p)
    if(style %in% c("acf", "density"))  par2plot <- rep(TRUE, 2+p)
    
  }
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
  }
  
  par(par.options)
  
  if(style == "chains"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(x@phi[ind, i], type = "l", ylab = bquote(phi[.(i)]), ...)
        if(priorMeans) abline(h = x@model$phi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[p+1]){
      plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
      if(priorMeans) abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[p+2]){
      plot(x@sigma2[ind], type = "l", ylab = expression(sigma^2), ...)
      if(priorMeans) abline(h = x@model$sigma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[p+3]){
      plot(x@t, x@Z, pch = 20, xlab = "t", ylab = "data", col = col.priorMean, ylim = range(c(x@Z, range(x@Y.est[ind, ]))))
      for(i in ind) lines(x@t, x@Y.est[i,])
      points(x@t, x@Z, pch = 20, col = col.priorMean)
    }  
  }
  if(style == "acf"){
    for(i in 1:p) {
      if(par2plot[i]){
        acf(x@phi[ind, i], xlab = bquote(phi[.(i)]), ...)
      }
    }
    if(par2plot[p+1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
    if(par2plot[p+2]) acf(x@sigma2[ind], xlab = expression(sigma^2), ...)
  }
  if(style == "density"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(density(x@phi[ind, i]), xlab = bquote(phi[.(i)]), ...)
        if(priorMeans) abline(v = x@model$phi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[p+1]){
      plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
      if(priorMeans) abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[p+2]){
      plot(density(x@sigma2[ind]), xlab = expression(sigma^2), ...)
      if(priorMeans) abline(v = x@model$sigma2, col = col.priorMean, lty = lty.priorMean)
    }
  }

  par(old.settings[which(names(old.settings) %in% names(par.options))])
})

#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the hidden hierarchical diffusion model.
#' @param x est.hiddenmixedDiffusion class, created with method \code{\link{estimate,hiddenmixedDiffusion-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density", "int.phi"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\mu, \Omega, \gamma^2, \sigma^2, Y)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param level level for style = "int.phi"
#' @param phi in the case of simulation study: known values for phi
#' @param ... optional plot parameters
#' @examples 
#' \dontrun{
#' mu <- c(10, 3, 1); Omega = c(1, 0.4, 0.01)
#' phi <- sapply(1:3, function(i) rnorm(20, mu[i], sqrt(Omega[i])))
#' model <- set.to.class("hiddenmixedDiffusion", b.fun = function(phi, t, y) phi[1]-phi[2]*y,
#'     parameter = list(mu = mu, Omega = Omega, phi = phi, gamma2 = 1, sigma2 = 0.1),
#'     y0 = function(phi, t) phi[3])
#' data <- simulate(model, t = seq(0, 1, by = 0.02), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.02), data$Z, 1000) 
#' plot(est, burnIn = 10, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2,1)), xlab = "iteration")
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, TRUE, rep(FALSE, 7)))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE, 
#'    par2plot = c(rep(FALSE, 6), TRUE, TRUE, FALSE))
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), main = "", par2plot = c(rep(FALSE, 6), TRUE, TRUE))
#' plot(est, style = "int.phi", phi = phi, par2plot = c(TRUE, FALSE, FALSE))
#' }
#' @export
setMethod(f = "plot", signature = "est.hiddenmixedDiffusion", 
          definition = function(x, par.options, style = c("chains", "acf", "density", "int.phi"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, level = 0.05, phi, ...) {
  old.settings <- par(no.readonly = TRUE)

  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
  }else{
    ind <- 1:length(x@gamma2)
  }
  
  p <- ncol(x@mu)
  k <- nrow(x@phi[[1]])
  
  if(missing(par2plot)){
    if(style == "chains")  par2plot <- rep(TRUE, 3+2*p)
    if(style %in% c("acf", "density"))  par2plot <- rep(TRUE, 2+2*p)
    if(style == "int.phi") par2plot <- rep(TRUE, p)
  }
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
    
    if(style == "int.phi") par.options <- list(mfrow = c(sum(par2plot[1:p]), 1))
  }
  
  par(par.options)
  
  if(style == "chains"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(x@mu[ind, i], type = "l", ylab = bquote(mu[.(i)]), ...)
        if(priorMeans) abline(h = x@model$mu[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    for(i in 1:p) {
      if(par2plot[p+i]){
        plot(x@Omega[ind, i], type = "l", ylab = bquote(omega[.(i)]), ...)
        if(priorMeans) abline(h = x@model$Omega[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[2*p+1]){
      plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
      if(priorMeans) abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[2*p+2]){
      plot(x@sigma2[ind], type = "l", ylab = expression(sigma^2), ...)
      if(priorMeans) abline(h = x@model$sigma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[2*p+3]){
      plot(x@t, x@Z[1,], pch = 20, xlab = "t", ylab = "first data series", col = col.priorMean,
           ylim = range(c(x@Z[1,], range(sapply(ind, function(i) x@Y.est[[i]][[1]])))))
      for(i in ind) lines(x@t, x@Y.est[[i]][[1]])
      points(x@t, x@Z[1,], pch = 20, col = col.priorMean)
    }  
  }
  if(style == "acf"){
    for(i in 1:p) {
      if(par2plot[i]) acf(x@mu[ind, i], xlab = bquote(mu[.(i)]), ...)
    }
    for(i in 1:p){
      if(par2plot[p+i]) acf(x@Omega[ind, i], xlab = bquote(omega[.(i)]), ...)
    }
    if(par2plot[2*p+1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
    if(par2plot[2*p+2]) acf(x@sigma2[ind], xlab = expression(sigma^2), ...)
  }
  if(style == "density"){
    for(i in 1:p) {
      if(par2plot[i]){
        plot(density(x@mu[ind, i]), xlab = bquote(mu[.(i)]), ...)
        if(priorMeans) abline(v = x@model$mu[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    for(i in 1:p) {
      if(par2plot[p+i]){
        plot(density(x@Omega[ind, i]), xlab = bquote(omega[.(i)]), ...)
        if(priorMeans) abline(v = x@model$Omega[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[2*p+1]){
      plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
      if(priorMeans) abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    if(par2plot[2*p+2]){
      plot(density(x@sigma2[ind]), xlab = expression(sigma^2), ...)
      if(priorMeans) abline(v = x@model$sigma2, col = col.priorMean, lty = lty.priorMean)
    }
  }
  if(style == "int.phi"){
    
    qu <- lapply(1:p, function(i) {
      apply( sapply(x@phi[ind], function(m) m[, i]), 1, function(ve) quantile(ve, c(level/2, 0.5, 1-level/2)) )
    })
    for(i in 1:p){
      if(par2plot[i]){
        plot(qu[[i]][2, ], pch = 20, ylim = range(qu[[i]]), ylab = bquote(phi[.(i)]), ...)
        segments(1:k, qu[[i]][1, ], 1:k, qu[[i]][3, ], ...)
        segments(1:k-0.1, qu[[i]][1, ], 1:k+0.1, qu[[i]][1, ], ...)
        segments(1:k-0.1, qu[[i]][3, ], 1:k+0.1, qu[[i]][3, ], ...)
        if(!missing(phi)) points(phi[, i], pch = 20, col = col.priorMean)
      }
    }
    par(old.settings[which(names(old.settings) %in% names(par.options))])
    return(invisible(qu))
  }
  

  par(old.settings[which(names(old.settings) %in% names(par.options))])
})

#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the jump regression model.
#' @param x est.jumpRegression class, created with method \code{\link{estimate,jumpRegression-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \theta, \gamma^2, \xi, N)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("jumpRegression", fun = function(t, N, theta) exp(theta[1]*t) + theta[2]*N,
#'   parameter = list(theta = c(2, 2), gamma2 = 0.25, xi = c(3, 0.5)),
#'   Lambda = function(t, xi) (t/xi[2])^xi[1])
#' data <- simulate(model, t = seq(0, 1, by = 0.01), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data, 1000)  # nMCMC small for example
#' plot(est)
#' plot(est, burnIn = 100, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2, 3)), xlab = "iteration")
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, FALSE, FALSE, TRUE, TRUE))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), par2plot = c(TRUE, rep(FALSE, 4)), main = "")
#' @export
setMethod(f = "plot", signature = "est.jumpRegression", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...) {
  old.settings <- par(no.readonly = TRUE)
  
  
  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
  }else{
    ind <- 1:length(x@gamma2)
  }
  
  p1 <- ncol(x@theta)
  p2 <- ncol(x@xi)
  if(style == "chains"){
    p3 <- c(1, 0)[c(sum(dim(x@N.est)) > 0, sum(dim(x@N.est)) == 0)]
  } else {
    p3 <- 0
  }
  if(missing(par2plot)) par2plot <- rep(TRUE, p1+p2+p3+1)
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
  }
  
  par(par.options)
  
  if(style == "chains"){
    
    for(i in 1:p1) {
      if(par2plot[i]){
        plot(x@theta[ind, i], type = "l", ylab = bquote(theta[.(i)]), ...)
        if(priorMeans) abline(h = x@model$theta[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[p1 + 1]){
      plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
      if(priorMeans)  abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    for(i in 1:p2){
      if(par2plot[p1 + 1 + i]){
        plot(x@xi[ind, i], type = "l", ylab = bquote(xi[.(i)]), ...)
        if(priorMeans)  abline(h = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    
    if(sum(dim(x@N.est) > 0)){
      if(par2plot[p1 +p2 + 1]){
        plot(x@t, x@N.est[,length(x@gamma2)], type = "l", xlab = "t", ylab = "N", ylim = range(x@N.est[,ind]), ...)
        for(i in ind) lines(x@t, x@N.est[,i])
      }  
    }
  }
  if(style == "acf"){
    
    for(i in 1:p1) {
      if(par2plot[i]) acf(x@theta[ind, i], xlab = bquote(theta[.(i)]), ...)
    }
    if(par2plot[p1 + 1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
    
    for(i in 1:p2){
      if(par2plot[p1 + 1 + + i]) acf(x@xi[ind, i], xlab = bquote(xi[.(i)]), ...)
    }
  }
  if(style == "density"){
    
    for(i in 1:p1) {
      if(par2plot[i]){
        plot(density(x@theta[ind, i]), xlab = bquote(theta[.(i)]), ...)
        if(priorMeans) abline(v = x@model$theta[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(par2plot[p1 + 1]){
      plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
      if(priorMeans)  abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
    }
    for(i in 1:p2){
      if(par2plot[p1 + 1 + i]){
        plot(density(x@xi[ind, i]), xlab = bquote(xi[.(i)]), ...)
        if(priorMeans)  abline(v = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }

  }

  par(old.settings[which(names(old.settings) %in% names(par.options))])
})

#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the NHPP.
#' @param x est.NHPP class, created with method \code{\link{estimate,NHPP-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \theta, \gamma^2, \xi, N)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)),
#'   Lambda = function(t, xi) (t/xi[2])^xi[1])
#' data <- simulate(model, t = seq(0, 1, by = 0.01), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data$Times, 10000)  # nMCMC small for example
#' plot(est)
#' plot(est, burnIn = 1000, thinning = 20, reduced = TRUE)
#' plot(est, xlab = "iteration")
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, FALSE), par.options = list(mfrow = c(1, 1)))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), par2plot = c(FALSE, TRUE), main = "")
#' @export
setMethod(f = "plot", signature = "est.NHPP", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...) {
  old.settings <- par(no.readonly = TRUE)
  
  
  style <- match.arg(style) 
  if(reduced){
    if(missing(thinning)) thinning <- x@thinning
    if(missing(burnIn)) burnIn <- x@burnIn
    ind <- seq(burnIn + 1, nrow(x@xi), by = thinning)
  }else{
    ind <- 1:nrow(x@xi)
  }
  
  p <- ncol(x@xi)

  if(missing(par2plot)) par2plot <- rep(TRUE, p)
  if(missing(par.options)){
    if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
    if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
  }
  
  par(par.options)
  
  if(style == "chains"){
    
    for(i in 1:p){
      if(par2plot[i]){
        plot(x@xi[ind, i], type = "l", ylab = bquote(xi[.(i)]), ...)
        if(priorMeans)  abline(h = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
    
  }
  if(style == "acf"){

    for(i in 1:p){
      if(par2plot[i]) acf(x@xi[ind, i], xlab = bquote(xi[.(i)]), ...)
    }
  }
  if(style == "density"){

    for(i in 1:p){
      if(par2plot[i]){
        plot(density(x@xi[ind, i]), xlab = bquote(xi[.(i)]), ...)
        if(priorMeans)  abline(v = x@model$xi[i], col = col.priorMean, lty = lty.priorMean)
      }
    }
  }
  
  par(old.settings[which(names(old.settings) %in% names(par.options))])
})

#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the regression model.
#' @param x est.Regression class, created with method \code{\link{estimate,Regression-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\phi, \gamma^2)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param ... optional plot parameters
#' @examples 
#' model <- set.to.class("Regression", fun = function(phi, t) phi[1]*t + phi[2],
#'     parameter = list(phi = c(1, 2), gamma2 = 0.1))
#' data <- simulate(model, t = seq(0, 1, by = 0.01), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.01), data, 1000)  # nMCMC small for example
#' plot(est)
#' plot(est, burnIn = 100, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(3,1)), xlab = "iteration")
#' plot(est, style = "acf", main = "", par2plot = c(TRUE, TRUE, FALSE))
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), main = "", par2plot = c(FALSE, FALSE, TRUE))
#' @export
setMethod(f = "plot", signature = "est.Regression", 
          definition = function(x, par.options, style = c("chains", "acf", "density"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, ...) {
    old.settings <- par(no.readonly = TRUE)
    
    style <- match.arg(style) 
    if(reduced){
      if(missing(thinning)) thinning <- x@thinning
      if(missing(burnIn)) burnIn <- x@burnIn
      ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
    }else{
      ind <- 1:length(x@gamma2)
    }
    
    p <- ncol(x@phi)
    
    if(missing(par2plot)) par2plot <- rep(TRUE, 1+p)
    if(missing(par.options)){
      if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
      if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
    }
    
    par(par.options)
    
    if(style == "chains"){
      for(i in 1:p) {
        if(par2plot[i]){
          plot(x@phi[ind, i], type = "l", ylab = bquote(phi[.(i)]), ...)
          if(priorMeans) abline(h = x@model$phi[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
      if(par2plot[p+1]){
        plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
        if(priorMeans) abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(style == "acf"){
      for(i in 1:p) {
        if(par2plot[i]){
          acf(x@phi[ind, i], xlab = bquote(phi[.(i)]), ...)
        }
      }
      if(par2plot[p+1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
    }
    if(style == "density"){
      for(i in 1:p) {
        if(par2plot[i]){
          plot(density(x@phi[ind, i]), xlab = bquote(phi[.(i)]), ...)
          if(priorMeans) abline(v = x@model$phi[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
      if(par2plot[p+1]){
        plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
        if(priorMeans) abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
      }
    }
    
    par(old.settings[which(names(old.settings) %in% names(par.options))])
})



#' Plot method for the Bayesian estimation results
#' 
#' @description Plot method for the estimation results of the hierarchical (mixed) regression model.
#' @param x est.mixedRegression class, created with method \code{\link{estimate,mixedRegression-method}}
#' @param par.options list of options for function par()
#' @param style one out of "chains", "acf", "density", "int.phi"
#' @param par2plot logical vector, which parameters to be plotted, order: \eqn{(\mu, \Omega, \gamma^2)}
#' @param reduced logical (1), if TRUE, the chains are thinned and burn-in phase is dropped
#' @param thinning thinning rate, if missing, the proposed one by the estimation procedure is taken
#' @param burnIn burn-in phase, if missing, the proposed one by the estimation procedure is taken
#' @param priorMeans logical(1), if TRUE (default), prior means are marked with a line 
#' @param col.priorMean color of the prior mean line, default 2
#' @param lty.priorMean linetype of the prior mean line, default 1
#' @param level level for style = "int.phi"
#' @param phi in the case of simulation study: known values for phi
#' @param ... optional plot parameters
#' @examples 
#' mu <- c(1, 3); Omega = c(0.4, 0.01)
#' phi <- sapply(1:2, function(i) rnorm(20, mu[i], sqrt(Omega[i])))
#' model <- set.to.class("mixedRegression", fun = function(phi, t) phi[1]*t + phi[2],
#'     parameter = list(mu = mu, Omega = Omega, phi = phi, gamma2 = 0.1))
#' data <- simulate(model, t = seq(0, 1, by = 0.02), plot.series = TRUE)
#' est <- estimate(model, t = seq(0, 1, by = 0.02), data, 100)  # nMCMC small for example
#' plot(est, burnIn = 10, thinning = 2, reduced = TRUE)
#' plot(est, par.options = list(mar = c(5, 4.5, 4, 2) + 0.1, mfrow = c(2,1)), xlab = "iteration")
#' plot(est, style = "acf", main = "")
#' plot(est, style = "density", lwd = 2, priorMean = FALSE)
#' plot(est, style = "density", col.priorMean = 1, lty.priorMean = 2, main = "posterior")
#' plot(est, style = "acf", par.options = list(), main = "", par2plot = c(rep(FALSE, 4), TRUE))
#' plot(est, style = "int.phi", phi = phi, par2plot = c(TRUE, FALSE))
#' @export
setMethod(f = "plot", signature = "est.mixedRegression", 
          definition = function(x, par.options, style = c("chains", "acf", "density", "int.phi"), par2plot, reduced = FALSE, 
                                thinning, burnIn, priorMeans = TRUE, col.priorMean = 2, lty.priorMean = 1, level = 0.05, phi, ...) {
    old.settings <- par(no.readonly = TRUE)
    
    style <- match.arg(style) 
    if(reduced){
      if(missing(thinning)) thinning <- x@thinning
      if(missing(burnIn)) burnIn <- x@burnIn
      ind <- seq(burnIn + 1, length(x@gamma2), by = thinning)
    }else{
      ind <- 1:length(x@gamma2)
    }
    
    p <- ncol(x@mu)
    k <- nrow(x@phi[[1]])
    
    if(missing(par2plot)) par2plot <- rep(TRUE, 1+2*p)
    if(missing(par.options)){
      if(sum(par2plot) == 1) par.options <- list(mfrow = c(1, 1))
      if(sum(par2plot) > 1) par.options <- list(mfrow = c(ceiling(sum(par2plot)/2), 2))
      if(style == "int.phi") par.options <- list(mfrow = c(sum(par2plot[1:p]), 1))
    }
    
    par(par.options)
    
    if(style == "chains"){
      for(i in 1:p) {
        if(par2plot[i]){
          plot(x@mu[ind, i], type = "l", ylab = bquote(mu[.(i)]), ...)
          if(priorMeans) abline(h = x@model$mu[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
      for(i in 1:p) {
        if(par2plot[p+i]){
          plot(x@Omega[ind, i], type = "l", ylab = bquote(omega[.(i)]), ...)
          if(priorMeans) abline(h = x@model$Omega[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
      if(par2plot[2*p+1]){
        plot(x@gamma2[ind], type = "l", ylab = expression(gamma^2), ...)
        if(priorMeans) abline(h = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(style == "acf"){
      for(i in 1:p) {
        if(par2plot[i]) acf(x@mu[ind, i], xlab = bquote(mu[.(i)]), ...)
      }
      for(i in 1:p){
        if(par2plot[p+i]) acf(x@Omega[ind, i], xlab = bquote(omega[.(i)]), ...)
      }
      if(par2plot[2*p+1]) acf(x@gamma2[ind], xlab = expression(gamma^2), ...)
    }
    if(style == "density"){
      for(i in 1:p) {
        if(par2plot[i]){
          plot(density(x@mu[ind, i]), xlab = bquote(mu[.(i)]), ...)
          if(priorMeans) abline(v = x@model$mu[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
      for(i in 1:p) {
        if(par2plot[p+i]){
          plot(density(x@Omega[ind, i]), xlab = bquote(omega[.(i)]), ...)
          if(priorMeans) abline(v = x@model$Omega[i], col = col.priorMean, lty = lty.priorMean)
        }
      }
      if(par2plot[2*p+1]){
        plot(density(x@gamma2[ind]), xlab = expression(gamma^2), ...)
        if(priorMeans) abline(v = x@model$gamma2, col = col.priorMean, lty = lty.priorMean)
      }
    }
    if(style == "int.phi"){
      
      qu <- lapply(1:p, function(i) {
        apply( sapply(x@phi[ind], function(m) m[, i]), 1, function(ve) quantile(ve, c(level/2, 0.5, 1-level/2)) )
      })
      for(i in 1:p){
        if(par2plot[i]){
          plot(qu[[i]][2, ], pch = 20, ylim = range(qu[[i]]), ylab = bquote(phi[.(i)]), ...)
          segments(1:k, qu[[i]][1, ], 1:k, qu[[i]][3, ], ...)
          segments(1:k-0.1, qu[[i]][1, ], 1:k+0.1, qu[[i]][1, ], ...)
          segments(1:k-0.1, qu[[i]][3, ], 1:k+0.1, qu[[i]][3, ], ...)
          if(!missing(phi)) points(phi[, i], pch = 20, col = col.priorMean)
        }
      }
      par(old.settings[which(names(old.settings) %in% names(par.options))])
      return(invisible(qu))
    }
    
    par(old.settings[which(names(old.settings) %in% names(par.options))])
})


