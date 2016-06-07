#'
#'
#' Transformation of event times to NHPP
#'
#' @description Transformation of vector of event times to the corresponding counting process variables.
#' @param times vector of event times
#' @param t times of counting process
#'
#' @export
TimestoN <- function(times, t){
  lt <- length(t)
  dN <- numeric(lt)
  ind <- unlist(sapply(times[times <= max(t)], function(a){which(abs(t-a) == min(abs(t-a)))}))
  for(a in 1:length(ind)) dN[ind[a]] <- dN[ind[a]] + 1
  cumsum(dN)
}


#' Transformation of NHPP variables to event times
#'
#' @description Vector of Poisson process differences are translated to a vector of event times.
#' @param dN vector of differences of counting process
#' @param t times of counting process
#' @export
dNtoTimes <- function(dN, t){
  if(any(dN > 1)){
    m <- sum(dN > 0)
    res <- NULL
    he1 <- dN[dN > 0]
    he2 <- t[dN > 0]
    for(mi in 1:m){
      res <- c(res, rep(he2[mi], he1[mi]))
    }
  }else{
    res <- t[dN > 0]
  }
  res
}



#' Sampling from lognormal proposal density
#'
#' @description Drawing one sample from the lognormal distribution with mean \code{parOld} and standard deviation \code{propSd}. Used in Metropolis Hastings algorithms.
#' @param parOld the parameter from the last iteration step
#' @param propSd proposal standard deviation
#'
#' @examples
#' plot(replicate(100, proposal(1, 0.1)), type = "l")
#' @export
proposal <- function(parOld, propSd){
  if(any(parOld < 1e-150)) parOld[parOld < 1e-150] <- 1e-150  # 1e-320 equal to zero ...
  mu <- log(parOld) - log( propSd^2/(parOld^2) + 1)/2
  sigma2 <- log( propSd^2/(parOld^2)+1)
  rlnorm(length(parOld), mu, sqrt(sigma2))
}

#' Proposal ratio of lognormal proposal density
#'
#' @description Calculation of proposal ratio, see also \code{\link{proposal}}.
#' @param parOld the parameter from the last iteration step
#' @param parNew drawn candidate
#' @param propSd proposal standard deviation
#'
#' @examples
#' cand <- proposal(1, 0.01)
#' proposalRatio(1, cand, 0.01)
#' @export
proposalRatio <- function(parOld, parNew, propSd){
  muOld <- log(parOld) - log( propSd^2/exp(2*log(parOld)) + 1)/2
  sigma2Old <- log( propSd^2/exp(2*log(parOld))+1)
  muNew <- log(parNew) - log( propSd^2/exp(2*log(parNew)) + 1)/2
  sigma2New <- log( propSd^2/exp(2*log(parNew))+1)

  prod(dlnorm(parOld, muNew, sqrt(sigma2New))/dlnorm(parNew, muOld, sqrt(sigma2Old)))
}


#' Inversion Method
#'
#' @description Algorithm to sample from cumulative distribution function, if no inverse function is analytically available.
#' @param Fun cumulative distribution function
#' @param len number of samples
#' @param candArea candidate area
#' @param grid fineness degree
#' @param method vectorial ("vector") or not ("free")
#' @examples
#' test <- InvMethod(function(x) pnorm(x, 5, 1), 1000, candArea = c(0, 10), method = "free")
#' plot(density(test))
#' curve(dnorm(x, 5, 1), col = 2, add = TRUE)
#' @references 
#' Devroye, L. (1986). Non-Uniform Random Variate Generation. New York: Springer.
#' @export
InvMethod <- function(Fun, len, candArea, grid = 1e-05, method = c("vector", "free")){
  method <- match.arg(method)
  if(missing(candArea)){
    candArea <- findCandidateArea(Fun)
    if(method == "vector") candArea <- seq(candArea[1], candArea[2], by = grid)
  }else{
    if(method == "vector" & length(candArea) == 2) candArea <- seq(candArea[1], candArea[2], by = grid)
    if(method == "free" & length(candArea) > 2) candArea <- c(min(candArea), max(candArea))
  }

  if(method == "vector"){
    diFu <- sapply(candArea, Fun)
    U <- runif(len, 0, max(diFu))
    res <- sapply(U, function(u) candArea[which(diFu >= u)[1]])
  }
  if(method == "free"){
    res <- numeric(len)
    U <- runif(len)
    for(i in 1:len){
      lower <- candArea[1]
      upper <- candArea[2]

      diff <- upper - lower
      while(diff >= grid){
        if(Fun(lower+diff/2) < U[i]){
          lower <- lower+diff/2
        }else{
          upper <- lower+diff/2
        }
        diff <- upper - lower
      }
      res[i] <- (lower+upper)/2
    }
  }
  res
}

#' Rejection Sampling Algorithm
#'
#' @description Algorithm to sample from an arbitrary density function.
#' @param Fun cumulative distribution function
#' @param dens density
#' @param len number of samples
#' @param cand candidate area
#' @param grid fineness degree
#' @param method vectorial ("vector") or not ("free")
#' @references 
#' Devroye, L. (1986). Non-Uniform Random Variate Generation. New York: Springer.
#' @examples
#' plot(density(RejSampling(dens = function(x) dnorm(x, 5, 1), 
#'    len = 500, cand = seq(2, 9, by = 0.001), method = "free")))
#' lines(density(RejSampling(dens = function(x) dnorm(x, 5, 1), len = 500, 
#'       cand = seq(2, 9, by = 0.001), method = "vector")), col=2)
#' curve(dnorm(x, 5, 1), from = 2, to = 8, add = TRUE, col = 3)
#' @export
RejSampling <- function(Fun, dens, len, cand, grid = 1e-03, method = c("vector", "free")){  # for negative support?!?
  method <- match.arg(method)
  if(method == "free"){
    res <- numeric(len)
    for(i in 1:len){
      if(missing(cand)){
        ca <- findCandidateArea(function(t) Fun(t))
      }else{
        ca <- range(cand)
      }
      mp <- optimize(f = function(t) dens(t), ca, maximum = TRUE)$objective
      resi <- NULL
      while(is.null(resi)){
        u <- runif(1,0,mp)
        candi <- runif(1, ca[1], ca[2])
        prob <- dens(candi)
        if(u <= prob){
          resi <- candi
        }
      }
      res[i] <- resi
    }
  }
  if(method == "vector"){
    res <- numeric(len)
    if(missing(cand)){
      ca <- findCandidateArea(Fun)
      cand <- seq(ca[1], ca[2], by = grid)
    }
    prob <- vapply(cand, function(v) dens(v), FUN.VALUE = numeric(1))
    cand <- cand[prob != 0]
    prob <- prob[prob != 0]
    mp <- max(prob)
    count <- 1
    while(count <= len){
      u <- runif(1, 0, mp)
      ind <- sample(1:length(cand), 1)
      if(u <= prob[ind]){
        res[count] <- cand[ind]
        count <- count + 1
      }
    }
  }
  return(res)
}


#' Adaptation of proposal standard deviation
#'
#' @description Adaptive MCMC: if acceptance rate of the chain is smaller than \code{lower} or larger than \code{upper}, 
#' the proposal standard deviation \code{propSd}=exp(l) is adapted with respect to function \code{delta.n}, 
#' that means, the new proposal standard deviation
#' is equal to exp(l-\code{delta.n(batch)}), respectively exp(l+\code{delta.n(batch)}).
#'  
#' @param chain Markov chain
#' @param propSd current proposal standard deviation
#' @param batch number of batch (of chain)
#' @param lower lower bound
#' @param upper upper bound
#' @param delta.n function of batch number
#'
#' @return adapted proposal standard deviation
#' @references Rosenthal, J. S. (2011). Optimal Proposal Distributions and Adaptive MCMC. In: Handbook of Markov Chain Monte Carlo, pp. 93-112.
#' @export
ad.propSd <- function(chain, propSd, batch, lower = 0.3, upper = 0.6, delta.n = function(n) min(0.05, 1/sqrt(n))){
  ar <- length(unique(chain))/length(chain)
  lsi <- log(propSd)

  lsi[ar < lower] <- lsi - delta.n(batch)
  lsi[ar > upper] <- lsi + delta.n(batch)
  exp(lsi)
}

findCandidateArea <- function(VFun, start = 1, pos.support = TRUE, quasi.null = 1e-05){
  if(pos.support){
    cand <- start
    while(1 - VFun(cand) > quasi.null){
      cand <- cand*2
    }
    upper <- cand

    cand <- start
    while(VFun(cand) > 0){
      cand <- cand/2
    }
    lower <- cand
  } else{
    cand <- start
    while(1 - VFun(cand) > quasi.null){
      cand <- cand*2
    }
    upper <- cand

    cand <- -start
    while(VFun(cand) > quasi.null){
      cand <- cand*2
    }
    lower <- cand
  }

  c(lower, upper)
}


#' Calculation of a proposal for burn-in phase and thinning rate
#'
#' @description The proposed burn-in is calculated by dividing the Markov chains into \code{m} blocks and calculate the 95\% credibility intervals and the respective mean. 
#' Starting in the first one, the block is taken as burn-in as long as the mean of the current block is not in the credibility interval of the following block or vice versa. 
#' The thinning rate is proposed by the first lag which leads to a chain autocorrelation less than \code{dependence}. 
#' It is not easy to automate these choices, so it is highly recommended to verify the chains manually.
#' @param chain vector of Markov chain samples
#' @param dependence allowed dependence for the chain
#' @param m number of blocks
#' @return vector of burn-in and thinning
#' @export
diagnostic <- function(chain, dependence = 0.8, m = 10) {
  lc <- length(chain)
  K <- floor(lc/m)
  thinning <- min(which(acf(chain[-(1:floor(lc/5))], plot=F)$acf <= dependence)[1], floor(K/10), na.rm = TRUE)
  he1 <- sapply(1:m, function(i) chain[((i - 1) * K + 1):(i * K)][seq(1, K, by = thinning)])

  he2 <- apply(he1, 2, quantile, c(0.025, 0.975))
  he.mean <- apply(he1, 2, mean)
  is.in <- (he.mean[-m] >= he2[1, -1] & he.mean[-m] <= he2[2, -1]) | (he.mean[-1] >= he2[1, -m] & he.mean[-1] <= he2[2, -m])
  #    burnIn <- 0
  burnIn <- K
  for (i in 1:(m-1)) {
    if (sum(is.in) < length(is.in)) {
      is.in <- is.in[-1]
      burnIn <- burnIn + K
    }
  }
  burnIn <- min((m-1)*K, burnIn)
  return(c(burnIn = burnIn, thinning = thinning))
}

