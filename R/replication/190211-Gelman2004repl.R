# Title:       Replicate findings in Gelman et al 2004 adn Gelman 2006, appendix C
# Description: Trying to replicate to understand relationship between OpenBUGS and R Gibbs samplers
#              Also trying to lay the basis for replicating Gelman 2006

# Set up
  library(R2OpenBUGS)
  # Set the Wine working directory and the directory to OpenBUGS
  WINE = "/usr/local/bin/wine"
  WINEPATH = "/usr/local/bin/winepath"
  OpenBUGS.pgm = "/Applications/OpenBUGS323/OpenBUGS.exe"

# C.2 Fitting a hierarchical model in Bugs --------------------------------
  # PDF page: 631
  # This is the same that happens in Gelman 2006 in
  # section 5.1: Noninformative prior disitbrution for the 8-schools problem (p.8)
#> Prep ####
  # DATA
    data("schools")
    J <- nrow(schools)
    y <- schools$estimate
    sigma.y <- schools$sd
  # Parameters definition
    parameters <- c ("theta", "mu.theta", "sigma.theta")
  data <- list ("J", "y", "sigma.y")
  
#> Model1: Uniform Prior ####
  # Model file
    model_Gelman2004 <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/BugsModels/Gelman2004-schools.txt"
  # Initial values (function)
    inits <- function(){
      list (theta       = rnorm(J, 0, 100), 
            mu.theta    = rnorm(1, 0, 100),
            sigma.theta = runif(1,0,100))
    }
  # Simulation
    schools.sim <- bugs (data, inits, parameters,
                         model.file = model_Gelman2004,
                         n.chains = 3, n.iter = 1000,
                         OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T) # very imp on mac
  # Results
    print(schools.sim)
    plot(schools.sim)
  # Access posterior simulations in R
    attributes(schools.sim$sims.list)
  # Plot posterior
    par(mfcol=c(3,1))
    hist(schools.sim$sims.list$sigma.theta,
         yaxt = "n", ylab = ".",
         xlim = c(0, 30),
         main = "8 Schools: posterior on sigma_a given uniform prior on sigma_a",
         breaks = 20)
  # Posterior predictive simulation
    theta <- schools.sim$sims.list$theta
    y.rep <- array (NA, c(1000, J))
    for (sim in 1:1000)
      y.rep[sim,] <- rnorm(J, theta[sim, ], sigma.y)
    
#> Model2: IG(1,1) Prior ####
  model_Gelman2004 <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/BugsModels/Gelman2004-schools-gammaPrior.txt"
  inits <- function (){
    list (theta     = rnorm(J, 0, 100), 
          mu.theta  = rnorm(1, 0, 100),
          tau.theta = runif(1, 0, 100)) # tau is the precision (which is 1/sigma)
                                        # this is different beacuse of the different prior
                                        # (and its parametrization)
  }
  schools.sim.GP <- bugs (data, inits, parameters,
                          model.file = model_Gelman2004,
                          n.chains = 3, n.iter = 1000,
                          OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T)
  print(schools.sim.GP)
  plot(schools.sim.GP)
  # Plot posterior
  hist(schools.sim.GP$sims.list$sigma.theta,
       yaxt = "n", ylab = ".",
       xlim = c(0, 30),
       main = "8 Schools: posterior on sigma_a given IG(1,1) prior on sigma_a^2",
       breaks = 20)
  
#> Model2: IG(.001, .001) Prior ####
  model_Gelman2004 <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/BugsModels/Gelman2004-schools-gammaPrior2.txt"
  inits <- function (){
    list (theta     = rnorm(J, 0, 100), 
          mu.theta  = rnorm(1, 0, 100),
          tau.theta = runif(1, 0, 100)) # tau is the precision (which is 1/sigma)
                                        # this is different beacuse of the different prior
                                        # (and its parametrization)
  }
  schools.sim.GP2 <- bugs (data, inits, parameters,
                          model.file = model_Gelman2004,
                          n.chains = 3, n.iter = 1000,
                          OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T)
  print(schools.sim.GP2)
  plot(schools.sim.GP2)
  # Plot posterior
  hist(schools.sim.GP2$sims.list$sigma.theta,
       yaxt = "n", ylab = ".",
       xlim = c(0, 30),
       main = "8J: post sigma_a w/ IG(.001, .001) prior on sigma_a^2",
       breaks = 50)
  
  
# C.4 Fitting a hierarchical model in R -----------------------------------
# PDF Pages: 640, and 163 (for Section 5.4)
  data("schools")
  J <- nrow(schools) # number of schools
  y <- schools$estimate # vecotr of data values
  sigma.y <- schools$sd # vecotr of standard deviations
  
#> Marginal and Conditional draws ###
  mu.hat <- function (tau, y, sigma.y) {
    sum(y/(sigma.y^2 + tau^2))/sum(1/(sigma.y^2 + tau^2))
  }
  V.mu <- function (tau, y, sigma.y){
    1/sum(1/(sigma.y^2 + tau^2))
  }
  
  # 1. Posteriror distribution of tau (marginal)
    # set up a grid for tau (tua = population btw-groups statndard deviation)
  n.grid <- 2000
  tau.grid <- seq(.01, 40, length = n.grid) # points equally spread from 0 to 40
  log.p.tau <- rep (NA, n.grid)
  for (i in 1:n.grid){
    mu <- mu.hat(tau.grid[i], y, sigma.y)
    V <- V.mu(tau.grid[i], y, sigma.y)
    log.p.tau[i] <- .5*log(V) - .5*sum(log(sigma.y^2 + tau.grid[i]^2)) -
      .5*sum((y-mu)^2/(sigma.y^2 + tau.grid[i]^2))
                                           # log posterior normal disitbrution
  }
  # Notes: compute the posterior density for tau on the log scale and rescale 
  # it to eliminate the possibility of computational overflow or underflow 
  # that can occur when multiplying many factors.
  log.p.tau <- log.p.tau - max(log.p.tau)
  p_tau <- exp(log.p.tau)
  p_tau <- p_tau/sum(p_tau)
  n_sims <- 1000
  tau <- sample (tau.grid, n_sims, replace=TRUE, prob=p_tau)
  
  # 2. Posterior distribution of mu given tau (conditional)
  mu <- rep (NA, n_sims)
  theta <- array (NA, c(n_sims,J))
  for (i in 1:n_sims){
    mu[i] <- rnorm (1, mu.hat(tau[i],y,sigma.y), sqrt(V.mu(tau[i],y,sigma.y)))
    theta_mean <- (mu[i]/tau[i]^2 + y/sigma.y^2)/(1/tau[i]^2 + 1/sigma.y^2)
    theta_sd <- sqrt(1/(1/tau[i]^2 + 1/sigma.y^2))
    theta[i,] <- rnorm (J, theta_mean, theta_sd)
  }
  
#> Gibbs Sampler ####
  theta_update <- function (){
    theta_hat <- (mu/tau^2 + y/sigma.y^2)/(1/tau^2 + 1/sigma.y^2)
    V_theta <- 1/(1/tau^2 + 1/sigma.y^2)
    rnorm (J, theta_hat, sqrt(V_theta))
  }
  mu_update <- function (){
    rnorm (1, mean(theta), tau/sqrt(J))
  }
  tau_update <- function (){
    sqrt(sum((theta-mu)^2)/rchisq(1,J-1))
  }
  
  chains <- 5 #generate five independent Gibbs sampling sequences of length 1000
  iter <- 1000
  sims <- array (NA, c(iter, chains, J+2))
  dimnames (sims) <- list (NULL, NULL,
                           c (paste ("theta[", 1:8, "]", sep=""), "mu", "tau"))
  for (m in 1:chains){
    mu <- rnorm(1, mean(y), sd(y))
    tau <- runif(1, 0, sd(y))
    for (t in 1:iter){
      theta <- theta_update()
      mu <- mu_update()
      tau <- tau_update()
      sims[t,m,] <- c(theta, mu, tau)
    }
  }
  # check the mixing of the sequences
  monitor(sims)
  rstan::monitor(sims)