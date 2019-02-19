# Title:       Replicate Hoff 2009 Chapter 7
# Description: GLM Mixed effects with inverse-Wishart distribution for the var-cov matrix

library(mvtnorm)

# Data (tumor example)
  XY.tumor <- dget("http://www.stat.washington.edu/~hoff/Book/Data/data/XY.tumor")
  Y <- XY.tumor$Y
  X <- XY.tumor$X
  m <- dim(Y)[1]
  p <- dim(X)[2]

# Priors
  # Beta multivariate normal prior
  BETA <- NULL
  for(j in 1:m) {
    BETA <- rbind(BETA, lm(log(Y[j, ]+1/20) ~ -1+X[,,j])$coef)
  }
  # 
  mu0 <- apply(BETA, 2, mean) # prior guess based on OLS estiamted values
  # Sigma inverse-Wishart
  S0  <- cov(BETA) # sample covairance matrix as prior guess
    iL0 <- iSigma <- solve(S0)
  eta0 <- p + 2 # usual rule for the first invwish parameter

## MCMC
set.seed(1)
THETA.post <- NULL
for(s in 1:50000) {
  ## update theta
  Lm <- solve( iL0 + m*iSigma )
  mum <- Lm %*% ( iL0%*%mu0 + iSigma%*%apply(BETA, 2, sum) ) # maybe multiply second part by m?
  theta <- t(rmvnorm(1, mum, Lm))
  
  ##update Sigma
  mtheta <- matrix(theta, m, p, 
                   byrow=TRUE)
  iSigma <- rwish(eta0 + m,
                  solve( S0 + t(BETA - mtheta)%*%(BETA - mtheta) ))
  
  ##update beta
  Sigma <- solve(iSigma) ; dSigma <- det(Sigma)
  for(j in 1:m) {
    # Sample of beta.j from proposal distirbution
    beta.p <- t( rmvnorm(1,
                         BETA[j ,],
                         .5*Sigma) )
    # Acceptance Ratio
    lr <- sum( dpois(Y[j,], exp(X[,,j]%*%beta.p),
                     log = TRUE ) -
                 dpois(Y[j,], exp(X[,,j]%*%BETA[j,]),
                       log = TRUE ) ) +
      dmvnorm( t(beta.p),
               theta,
               Sigma,
               log = TRUE) -
      dmvnorm( t(BETA[j,]),
               theta,
               Sigma,
               log = TRUE )
    
    if( log(runif(1)) < lr ) { BETA[j,] <- beta.p }
  }
  ##store some output (only 1 every 10 samples)
  if(s %% 10 == 0){THETA.post <- rbind(THETA.post, t( theta ))}
}
