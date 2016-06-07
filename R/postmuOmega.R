

# postmu <- function(phi, m, v, Omega){  # phi matrix
#   n <- nrow(phi)
#   V_ <- diag(1/v)
#   Vpost <- solve(V_ + n*solve(Omega))
#   mpost <- Vpost%*%( V_%*%m + apply((solve(Omega)%*%t(phi)),1,sum) )
#
#   rmvnorm(1,mpost,Vpost)
# }
# in the case of Omega = vector... possibly faster...
postmu <- function(phi, m, v, Omega){  # phi nxp-matrix, m mean of mu, v diagonal of variance matrix
  n <- nrow(phi)
  V_ <- 1/v
  Vpost <- 1/(V_ + n/Omega)
  mpost <- Vpost*( m/v + apply(t(phi)/Omega, 1, sum) )

  rnorm(length(mpost), mpost, sqrt(Vpost))
}

postOmega <- function(alpha, beta, phi, mu){  # length(alpha)=length(beta)=length(mu)
  p <- length(mu)
  Dia <- numeric(p)
  for(i in 1:p){
    Dia[i] <- 1/rgamma(1, alpha[i] + nrow(phi)/2, beta[i] + sum((phi[,i]-mu[i])^2)/2)
  }
  Dia
}


# postOmega_matrix <- function(R, phi, mu){
#   Rpost <- solve(R + (t(phi)-as.vector(mu))%*%t((t(phi)-as.vector(mu))))
#   solve( rWishart(1,nrow(phi)+length(mu)+1,Rpost)[,,1])
# }
