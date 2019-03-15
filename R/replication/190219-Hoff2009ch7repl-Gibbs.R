# Title:       Replicate Hoff 2009 Chapter 7
# Description: Multivariate Normal Model to learn how to use the inverse-Wishart distribution

library(MCMCpack)
library(mvtnorm)

# Models ------------------------------------------------------------------
#Data
  #Manually copied from plot on page 114
  y1 <- c(27, 32, 34, 34, 34, 34, 35, 38, 41, 42, 43, 44, 50, 51, 53, 55, 59, 60, 63, 70, 72)
  y2 <- c(58, 27, 38, 47, 48, 55, 34, 43, 55, 38, 39, 60, 55, 57, 60, 69, 77, 75, 86, 55, 59)
  n     <- length(y1)         # number of observations
  Y <- cbind(y1, y2)
  Sigma <- cov(cbind(y1, y2))            # covariance between the variables at hand
  ybar  <- apply(Y, 2, mean) # sample group means

#Priors
  mu0 <-  c(50, 50) # guess for the prior menas on the two (p) variables
  L0  <- matrix(c(625, 312.5,
                  312.5, 625),
                nrow = 2, ncol = 2)
  nu0 <- 4 # based on p + 2 logic
  S0  <- matrix(c(625 , 312.5,
                  312.5 ,625),
                nrow = 2, ncol = 2)
# Simulation
  THETA <- SIGMA <- NULL     # containers for simulation
  set.seed(1) 
  for(s in 1:5000) {
    ###update theta
    s <- 1
    Ln    <- solve ( solve(L0) + n*solve(Sigma) ) # An
    mun   <- Ln %*% ( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
    theta <- rmvnorm(1, mun, Ln)
  
    ###update Sigma
    Sn    <- S0 + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta) )
    Sigma <- solve( rwish(v = nu0 + n,
                          S = solve(Sn)) )
  
    ### save results
    THETA <- rbind(THETA, theta)
    SIGMA <- rbind(SIGMA, c(Sigma))
  }
  
# Credibility Intervals
quantile( THETA[,2]-THETA[ ,1] , prob=c(.025 ,.5 ,.975) )

# Posterior probability of average score in second exam higher than in first
mean( THETA[,2]>THETA[,1])
  