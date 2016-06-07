

########
#' Simulation of diffusion process
#'
#' @description Simulation of a stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t}.
#' @param object class object of parameters: "Diffusion"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param y0 starting point of the process
#' @param mw mesh width for finer Euler approximation to simulate time-continuity
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "Diffusion",
          definition = function(object, nsim = 1, seed = NULL, t, y0, mw = 1, plot.series = TRUE) {
            set.seed(seed)
            if(nsim > 1){
              result <- matrix(0, nsim, length(t))
              for(i in 1:nsim){
                result[i,] <- drawSDE(object@phi, object@gamma2, t, object@b.fun, y0.fun = function(phi, t) y0,
                                  sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
              }
              if(plot.series){
                plot(t, result[1,], type = "l", ylab = "Y", ylim = range(result))
                for(i in 2:nsim) lines(t, result[i,], col = i)
              }
            }else{
              result <- drawSDE(object@phi, object@gamma2, t, object@b.fun, y0.fun = function(phi, t) y0,
                                sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
              if(plot.series){
                plot(t, result, type = "l", ylab = "Y")
              }
            }

    return(result)
})


########
#' Simulation of hierarchical (mixed) diffusion model
#'
#' @description Simulation of the stochastic process model
#'   \eqn{dY_t = b(\phi_j,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param object class object of parameters: "mixedDiffusion"
#' @param nsim number of data sets to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param mw mesh width for finer Euler approximation to simulate time-continuity
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' model <- set.to.class("mixedDiffusion", y0.fun = function(phi, t) 0.5, 
#'   parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'   b.fun = function(phi, t, x) phi*x, sT.fun = function(t, x) x)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "mixedDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, mw = 1, plot.series = TRUE) {
      set.seed(seed)
      if(nsim > 1){
        result <- list()
        for(j in 1:nsim){
          Y <- matrix(0, nrow(object@phi), length(t))
          for(i in 1:nrow(object@phi)){
            Y[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, y0.fun = object@y0.fun,
                                  sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
          }
          result[[j]] <- Y
        }
        if(plot.series){
          plot(t, result[[1]][1,], type = "l", ylab = "Y", ylim = range(result[[1]]))
          for(i in 2:nrow(object@phi)) lines(t, result[[1]][i,])
        }
      }else{
        result <- matrix(0, nrow(object@phi), length(t))
        for(i in 1:nrow(object@phi)){
          result[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, y0.fun = object@y0.fun,
                                sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
        }
        if(plot.series){
          plot(t, result[1,], ylim = range(result), type = "l", ylab="Y")
          for(i in 2:nrow(object@phi)) lines(t, result[i,])
        }
      }

      return(result)
})



########
#' Simulation of regression model
#'
#' @description Simulation of the regression model
#'   \eqn{y_i = f(\phi, t_i) + \epsilon_i, \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t_i))}.
#' @param object class object of parameters: "Diffusion"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("Regression", parameter = list(phi = 5, gamma2 = 0.1), 
#'    fun = function(phi, t) phi*t)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "Regression",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            if(nsim > 1){
              result <- matrix(0, nsim, length(t))
              for(i in 1:nsim){
                result[i,] <- object@fun(object@phi, t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
              }
              if(plot.series){
                plot(t, result[1,], ylab = "Y", ylim = range(result))
                for(i in 2:nsim) points(t, result[i,], col = i)
              }
            }else{
              result <- object@fun(object@phi, t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
              if(plot.series){
                plot(t, result, ylab = "Y")
              }

            }
            return(result)
          })


########
#' Simulation of hierarchical (mixed) regression model
#'
#' @description Simulation of regression model
#'   \eqn{y_{ij} = f(\phi_j, t_{ij}) + \epsilon_{ij}, \phi_j\sim N(\mu, \Omega),
#'   \epsilon_{ij}\sim N(0,\gamma^2\widetilde{s}(t_{ij}))}.
#' @param object class object of parameters: "mixedRegression"
#' @param nsim number of data sets to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' model <- set.to.class("mixedRegression", 
#'    parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'    fun = function(phi, t) phi*t, sT.fun = function(t) t)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "mixedRegression",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            if(nsim > 1){
              result <- list()
              for(j in 1:nsim){
                Y <- matrix(0, nrow(object@phi), length(t))
                for(i in 1:nrow(object@phi)){
                  Y[i,] <- object@fun(object@phi[i,], t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
                }
                result[[j]] <- Y
              }
              if(plot.series){
                plot(t, result[[1]][1,], type = "l", ylab = "Y", ylim = range(result[[1]]))
                for(i in 2:nrow(object@phi)) lines(t, result[[1]][i,])
              }
            }else{
              result <- matrix(0, nrow(object@phi), length(t))
              for(i in 1:nrow(object@phi)){
                result[i,] <- object@fun(object@phi[i,], t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
              }
              if(plot.series){
                plot(t, result[1,], ylim = range(result), type = "l", ylab="Y")
                for(i in 2:nrow(object@phi)) lines(t, result[i,])
              }
            }

            return(result)
          })



########
#' Simulation of hidden diffusion process
#'
#' @description Simulation of a hidden stochastic process model
#'   \eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, 
#'   \epsilon_i\sim N(0,\sigma^2), Y_{t_0}=y_0(\phi, t_0)}.
#' @param object class object of parameters: "hiddenDiffusion"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param mw mesh width for finer Euler approximation to simulate time-continuity
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("hiddenDiffusion", parameter = list(phi = 0.5, gamma2 = 0.01, sigma2 = 0.1))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "hiddenDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, mw = 10, plot.series = TRUE) {
            set.seed(seed)
            if(nsim > 1){
              result <- matrix(0, nsim, length(t))
              for(i in 1:nsim){
                result[i,] <- drawSDE(object@phi, object@gamma2, t, object@b.fun, y0.fun = object@y0.fun,
                                      sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
              }
              obs <- apply(result, 2, function(vec) rnorm(length(vec), vec, sqrt(object@sigma2)))
              if(plot.series){
                plot(t, result[1,], type = "l", ylab = "Y", ylim = range(obs))
                for(i in 2:nsim) lines(t, result[i,], col = i)
                for(i in 1:nsim) points(t, obs[i,], col = i)
              }
            }else{
              result <- drawSDE(object@phi, object@gamma2, t, object@b.fun, y0.fun = object@y0.fun,
                                sT = object@sT.fun, mw = mw, strictly.positive = FALSE)

              obs <- rnorm(length(t), result, sqrt(object@sigma2))
              if(plot.series){
                plot(t, obs, ylab="simulations"); lines(t, result)
                legend("topleft", "hidden process", col = 1, lty = 1, box.lty = 0, inset = 0.01)
              }
              result <- list(Z = obs, Y = result)

            }

            return(result)
          })


########
#' Simulation of hierarchical (mixed) hidden diffusion model
#'
#' @description Simulation of a stochastic process
#'   \eqn{Z_{ij} = Y_{t_{ij}} + \epsilon_{ij}, dY_t = b(\phi_j,t,Y_t)dt + \gamma \widetilde{s}(t,Y_t)dW_t, \phi_j\sim N(\mu, \Omega), 
#'   Y_{t_0}=y_0(\phi, t_0), \epsilon_{ij}\sim N(0,\sigma^2)}.
#' @param object class object of parameters: "hiddenmixedDiffusion"
#' @param nsim number of data sets to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param mw mesh width for finer Euler approximation to simulate time-continuity
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' mu <- c(5, 1); Omega <- c(0.9, 0.04)
#' phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' y0.fun <- function(phi, t) phi[2]
#' model <- set.to.class("hiddenmixedDiffusion", y0.fun = y0.fun, 
#'    b.fun = function(phi, t, y) phi[1], 
#'    parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 1, sigma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t)
#' @export
setMethod(f = "simulate", signature = "hiddenmixedDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, mw = 10, plot.series = TRUE) {
            set.seed(seed)
            if(nsim > 1){
              result1 <- list()
              for(j in 1:nsim){
                Y <- matrix(0, nrow(object@phi), length(t))
                for(i in 1:nrow(object@phi)){
                  Y[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, y0.fun = object@y0.fun,
                                   sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
                }
                result1[[j]] <- Y
              }
              result2 <- lapply(result1, function(Y) apply(Y, 2, function(vec) rnorm(length(vec), vec, sqrt(object@sigma2))))
              if(plot.series){
                old.settings <- par(no.readonly = TRUE)
                
                par(mfrow = c(2,1))
                plot(t, result1[[1]][1,], type = "l", ylab = "Y", ylim = range(result2[[1]]))
                for(i in 2:nrow(object@phi)) lines(t, result1[[1]][i,])
                plot(t, result2[[1]][1,], type = "l", ylab = "Z", ylim = range(result2[[1]]))
                for(i in 2:nrow(object@phi)) lines(t, result2[[1]][i,])
                par(old.settings)
                
              }
            }else{
              result1 <- matrix(0, nrow(object@phi), length(t))
              for(i in 1:nrow(object@phi)){
                result1[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, y0.fun = object@y0.fun,
                                      sT = object@sT.fun, mw = mw, strictly.positive = FALSE)
              }

              result2 <- apply(result1, 2, function(vec) rnorm(length(vec), vec, sqrt(object@sigma2)))
              if(plot.series){
                old.settings <- par(no.readonly = TRUE)
                
                par(mfrow = c(2,1))
                plot(t, result1[1,], ylim = range(result2), type = "l", ylab = "Y")
                for(i in 2:nrow(object@phi)) lines(t, result1[i,])
                plot(t, result2[1,], ylim = range(result2), type = "l", ylab = "Z")
                for(i in 2:nrow(object@phi)) lines(t, result2[i,])
                par(old.settings)
                
              }
              result <- list(Z = result2, Y = result1)

           }

            return(result)
          })




########
#' Simulation of Poisson process
#'
#' @description Simulation of non-homogeneous Poisson process with cumulative intensity function \eqn{\Lambda(t, \xi)}.
#' @param object class object of parameters: "NHPP"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)), 
#'                       Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t)
#' @export
setMethod(f = "simulate", signature = "NHPP",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            result <- simN(t, object@xi, len = nsim, Lambda = object@Lambda)
            if(plot.series){
              if(nsim > 1){
                plot(t, result$N[1,], ylab="N", type="l", ylim = range(result$N))
                for(i in 2:nsim) lines(t, result$N[i,])
              }else{
                plot(t, result$N, ylab="N", type="l")
              }
            }

            return(result)
          })



########
#' Simulation of jump diffusion process
#'
#' @description Simulation of jump diffusion process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t + h(\eta,t,Y_t)dN_t}.
#' @param object class object of parameters: "jumpDiffusion"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param y0 starting point of process
#' @param start vector: start[1] starting point time, start[2] starting point for Poisson process
#' @param mw mesh width for finer Euler approximation to simulate time-continuity
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("jumpDiffusion", 
#'    parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)), 
#'    Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5)
#' @export
setMethod(f = "simulate", signature = "jumpDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, y0, start = c(0,0), mw = 1, plot.series = TRUE) {
    set.seed(seed)
    N <- simN(t, object@xi, len = nsim, start = start, Lambda = object@Lambda)$N
    
    phi <- object@phi
    theta <- object@theta
    gamma2 <- object@gamma2
    b.fun <- object@b.fun
    s.fun <- object@s.fun
    h.fun <- object@h.fun
    start <- y0
    
    lt <- length(t)
    lt2 <- (lt-1)*mw+1
    t2 <- seq(min(t), max(t), length = lt2)
    dt2 <- t2[2] - t2[1]
    
    if(nsim > 1){
      dN <- t(apply(N, 1, diff))
      number <- nrow(N)
      dN2 <- matrix(0, number, lt2-1)
      dN2[, seq(1, lt2-1, by = mw)] <- dN
      Y <- matrix(y0, number, lt2)
      for(i in 2:lt2){
        W <- rnorm(number, 0, sqrt(dt2))
        Y[,i] <- Y[,i-1] + b.fun(phi, t2[i-1], Y[,i-1])*dt2 + s.fun(gamma2, t2[i-1], Y[,i-1])*W +
          h.fun(theta, t2[i-1], Y[,i-1])*dN2[,i-1]
      }
      Y <- Y[, seq(1, lt2, by = mw)]
    }
    if(nsim == 1){
      Y <- rep(y0, lt2)
      dN2 <- numeric(lt2-1)
      dN2[seq(1, lt2-1, by = mw)] <- diff(N)
      for(i in 2:lt2){
        W <- rnorm(1, 0, sqrt(dt2))
        Y[i] <- Y[i-1] + b.fun(phi, t2[i-1], Y[i-1])*dt2 + s.fun(gamma2, t2[i-1], Y[i-1])*W +
          h.fun(theta, t2[i-1], Y[i-1])*dN2[i-1]
      }
      Y <- Y[seq(1, lt2, by = mw)]
    }
   
    result <- list(N = N, Y = Y)
    if(plot.series){
      old.settings <- par(no.readonly = TRUE)
      
      if(nsim > 1){
        
        par(mfrow = c(2,1))
        plot(t, N[1,], type = "l", ylab = "N", ylim = range(N))
        for(i in 2:nsim) lines(t, N[i,])
        plot(t, result$Y[1,], type = "l", ylab = "Y", ylim = range(result$Y))
        for(i in 1:nsim) lines(t, result$Y[i,])
        
      }else{
        
        par(mfrow = c(2,1))
        plot(t, N, type = "l")
        plot(t, result$Y, type = "l", ylab = "Y")
        
      }
      par(old.settings)
      
    }
    return(result)
})





########
#' Simulation of jump diffusion process
#'
#' @description Simulation of jump diffusion process
#'   \eqn{Y_t = y_0 \exp( \phi t - \gamma^2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @param object class object of parameters: "Merton"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param y0 starting point of process
#' @param start vector: start[1] starting point time, start[2] starting point for Poisson process
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("Merton", parameter = list(thetaT = 0.1, phi = 0.05, gamma2 = 0.1, xi = 10))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5)
#' @export
setMethod(f = "simulate", signature = "Merton",
          definition = function(object, nsim = 1, seed = NULL, t, y0, start = c(0,0), plot.series = TRUE) {
    set.seed(seed)

    N <- simN(t, object@xi, len = nsim, start = start, Lambda = object@Lambda)$N
    
    phi <- object@phi
    thetaT <- object@thetaT
    gamma2 <- object@gamma2
    
    
    l <- length(t)
    if(nsim > 1){
      Y <- t(sapply(1:nsim, function(i){
        y0*exp( (phi-gamma2/2)*(t-start[1])+sqrt(gamma2)*cumsum(c(0, rnorm(l-1, 0, sqrt(diff(t[-l])))))+thetaT*(N[i,]-start[2]) )
      }) )
    }
    if(nsim == 1){
      Y <- y0*exp( (phi-gamma2/2)*(t-start[1])+sqrt(gamma2)*cumsum(c(0, rnorm(l-1, 0, sqrt(diff(t[-l])))))+thetaT*(N-start[2]) )
    }

    result <- list(N = N, Y = Y)

    if(plot.series){
      old.settings <- par(no.readonly = TRUE)
      
      if(nsim > 1){
        par(mfrow = c(2,1))
        plot(t, N[1,], type = "l", ylab = "N", ylim = range(N))
        for(i in 2:nsim) lines(t, N[i,])
        plot(t, result$Y[1,], type = "l", ylab = "Y", ylim = range(result$Y))
        for(i in 1:nsim) lines(t, result$Y[i,])
      }else{
        par(mfrow = c(2,1))
        plot(t, N, type = "l")
        plot(t, result$Y, type = "l", ylab = "Y")
      }
      par(old.settings)
      
    }
    return(result)
})



########
#' Simulation of regression model dependent on Poisson process
#'
#' @description Simulation of of the regression model
#'   \eqn{y_i = f(t_i, N_{t_i}, \theta) + \epsilon_i} with
#'   \eqn{N_t\sim Pois(\Lambda(t, \xi)), \epsilon_i\sim N(0,\gamma^2\widetilde{s}(t))}.
#' @param object class object of parameters: "jumpRegression"
#' @param nsim number of trajectories to simulate. Default is 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("jumpRegression", fun = function(t, N, theta) theta[1]*t + theta[2]*N, 
#'    parameter = list(theta = c(1,2), gamma2 = 0.1, xi = 10))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t)
#' @export
setMethod(f = "simulate", signature = "jumpRegression",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

    N <- simN(t, object@xi, len = nsim, Lambda = object@Lambda)$N
    if(nsim > 1){
      result <- apply(N, 2, function(Nt) object@fun(t, Nt, object@theta) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t))))
    }else{
      result <- object@fun(t, N, object@theta) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
    }

    result <- list(N = N, Y = result)

    if(plot.series){
      old.settings <- par(no.readonly = TRUE)
      
      if(nsim > 1){
        par(mfrow = c(2,1))
        plot(t, N[1,], type = "l", ylab = "N", ylim = range(N))
        for(i in 2:nsim) lines(t, N[i,])
        plot(t, result$Y[1,], type = "l", ylab = "Y", ylim = range(result$Y))
        for(i in 1:nsim) lines(t, result$Y[i,])
      }else{
        par(mfrow = c(2,1))
        plot(t, N, type = "l")
        plot(t, result$Y, ylab = "Y")
      }
      par(old.settings)
      
    }
    return(result)
})


drawSDE <- function(phi, gamma2, t, b, y0.fun, sT, mw = 10, strictly.positive = TRUE){
  lt <- length(t)
  lt2 <- (lt-1)*mw+1
  t2 <- seq(min(t), max(t), length = lt2)
  dt2 <- t2[2] - t2[1]
  X <- numeric(lt2)
  X[1] <- -1

  if(strictly.positive){
    while(any(X < 0)){
      err <- rnorm(lt2-1, 0, sqrt(dt2))
      X[1] <- y0.fun(phi, t2[1])
      for(j in 2:lt2){
        X[j] <- X[j-1] + b(phi, t2[j-1], X[j-1])*dt2 + sqrt(gamma2)*sT(t2[j-1], X[j-1])*err[j-1]
      }
    }
  } else {
    err <- rnorm(lt2-1, 0, sqrt(dt2))
    X[1] <- y0.fun(phi, t2[1])
    for(j in 2:lt2){
      X[j] <- X[j-1] + b(phi, t2[j-1], X[j-1])*dt2 + sqrt(gamma2)*sT(t2[j-1], X[j-1])*err[j-1]
    }
  }

  X[seq(1, lt2, by = mw)]
}

simN <- function(t, xi, len, start = c(0,0), Lambda, int = c("Weibull","Exp")){  # start=(T_n,n)
  if(missing(Lambda)){
    int <- match.arg(int)
    if(int == "Weibull"){
      Lambda <- function(t){
        (t/xi[2])^xi[1]
      }
      F_1 <- function(u, Tn_1){
        xi[2]*(Lambda(Tn_1)-log(1-u))^(1/xi[1])
      }
    }else{
      F_1 <- function(u, Tn_1){
        (log(exp(xi[1]*Tn_1+xi[2])-log(1-u))-xi[2])/xi[1]
      }
    }
    drawTn <- function(Tn_1){
      help <- F_1(runif(1), Tn_1)
      if(help < t[1]){
        return(t[1])
      }else{
        return(t[which(min(abs(help-t)) == abs(help-t))[1]])
      }
    }
  } else {

    drawTn <- function(Tn_1){
      prob <- 1-exp(-(Lambda(t, xi)-Lambda(Tn_1, xi)))
      u <- runif(1)
      t[which(abs(u-prob) == min(abs(u-prob)))]
    }
  }
  lt <- length(t)
  drawTn_vec <- function(Tn_1){
    sapply(Tn_1, drawTn)
  }

  Tn <- drawTn_vec(rep(start[1], len))
  times <- Tn

  while(min(Tn, na.rm = T) < t[lt]){
    ind <- which(Tn < t[lt])
    Tn[ind] <- drawTn_vec(Tn[ind])
    Tn[-ind] <- NA
    if(all(is.na(Tn))) break
    times <- cbind(times, Tn)
  }
  if(len > 1){
    N_out <- t(apply(times, 1, TimestoN, t)) + start[2]
    Times <- times
  }else{
    N_out <- TimestoN(times, t) + start[2]
    Times <- as.numeric(times)
  }
  list(Times = Times, N = N_out)
}


