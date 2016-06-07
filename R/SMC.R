

####################
# helping functions...

A.to.B <- function(A){
  lt <- dim(A)[2] + 1;
  Npart <- dim(A)[1];
  B <- matrix(0, Npart, lt)
  B[,lt] <- 1:Npart;
  for (t in seq(lt, 2, by = -1)){
    B[,t-1] <- A[B[,t],t-1]}
  return(B)
}

SMC <- function(phi, gamma2, sigma2, Npart, times, Z, b.fun, y0.fun, sigmaTilde, conditional = TRUE, Y.fixed, B.fixed){

  lt <- length(times)
  dt <-  diff(times)
  y <- matrix(y0.fun(phi, times[1]), Npart, lt)
  w <- Weights <- matrix(1,Npart,lt)
  parents <-  matrix(1, Npart, lt-1)

  w[,1] <- dnorm(Z[1], mean = y[,1], sd = sqrt(sigma2))
  Weights[,1] <- w[,1]/sum(w[,1]) # or only 1/Npart...
  if(conditional){
    y[B.fixed[1],1] <- Y.fixed[1]
  }else{
    if(missing(Y.fixed)) Y.fixed <- numeric(lt)
  }

  for(n in 2:lt){
    if(conditional){
      ##############################
      # sampling of A_{n-1}^{-B_n^K}
      ##############################
      set.parents  <- (1:Npart)[-B.fixed[n]]

      On_1 <- rmultinom(1, Npart-1, Weights[,n-1])
      O <- On_1[B.fixed[n-1]] + 1
      he <- sample(set.parents, O-1)

      parents[B.fixed[n], n-1] <- B.fixed[n-1]
      parents[he, n-1] <- B.fixed[n-1]
      parents[-c(B.fixed[n],he), n-1] <- sample(set.parents, Npart-O, replace = TRUE, prob =  Weights[-B.fixed[n],n-1])
    }else{
      set.parents <- 1:Npart
      parents[, n-1] <- sample(1:Npart, Npart, replace = TRUE, prob = Weights[,n-1])
    }

    y[,1:(n-1)] <- y[parents[,n-1], 1:(n-1)]
    y.past <- y[,n-1]

    ##########################################
    # simulation of particles
    ##########################################

    Vpost <- 1/(1/(gamma2*sigmaTilde(times[n-1], y.past)^2*dt[n-1]) + 1/sigma2)
    mpost <- Vpost*((y.past+b.fun(phi, times[n-1], y.past)*dt[n-1])/(gamma2*sigmaTilde(times[n-1], y.past)^2*dt[n-1])
                    + Z[n]/sigma2)
    y.new <- rnorm(Npart, mpost, sqrt(Vpost))
    y[set.parents,n] <- y.new[set.parents]
    y[-set.parents,n] <- Y.fixed[n]
    wn <- dnorm(y[,n], y.past+b.fun(phi, times[n-1], y.past)*dt[n-1], sqrt(gamma2*sigmaTilde(times[n-1], y.past)^2*dt[n-1]))*
      dnorm(Z[n], y[,n], sqrt(sigma2))/dnorm(y[,n], mpost, sqrt(Vpost) )

    w[,n] <- wn
    Weights[,n] <- w[,n]/sum(w[,n])
  }
  return(list(W = Weights, y = y , parents = parents) )
}

