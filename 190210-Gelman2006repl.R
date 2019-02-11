# Replicate findings in Gelman 2006
# Trying to replicate what was done in Gelman2006

# Sources: Gelman2006
# Set up
  library(R2OpenBUGS)
  # Set the Wine working directory and the directory to OpenBUGS
  WINE = "/usr/local/bin/wine"
  WINEPATH = "/usr/local/bin/winepath"
  OpenBUGS.pgm = "/Applications/OpenBUGS323/OpenBUGS.exe"

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
       main = "8 Schools: posterior on sigma_a given IG(1,1) prior on sigma_a**2",
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
       main = "8 Schools: posterior on sigma_a given IG(.001,.001) prior on sigma_a**2",
       breaks = 50)
  
# BUGS model: Half Cauchy prior (half-t w/ df nu = 1) ---------------------
  # Prep
    # DATA
      data("schools")
      J <- 3
      y <- schools$estimate
      sigma.y <- schools$sd
    # Prior hyperparameter for half-t distribution
      prior.scale <- 25 # scale parameter for the half-Cauchy distirbution (see text for decision)
    # Parameters definition
      parameters <- c ("theta", "mu.theta", "sigma.theta")
    data <- list ("J", "y", "sigma.y", "prior.scale")
    # Model file
    model_Gelman2006_HCP <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/BugsModels/Gelman2006-halfCauchyPrior.txt"
      # why relative path doesn't work?
  # Initial values (function)
    inits <- function (){
      list (eta      = rnorm(J), 
            mu.theta = rnorm(1), 
            xi       = rnorm(1), 
            tau.eta  = runif(1))}
  # Simulation
    schools.sim <- bugs (data, inits, parameters,
                         model.file = model_Gelman2006_HCP,
                         n.chains = 3, n.iter = 1000,
                         OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T) # very imp on mac
  # Results
  print(schools.sim)
  plot(schools.sim)

  # Access posterior simulations in R
  attributes(schools.sim$sims.list)
  hist(schools.sim$sims.list$sigma.theta,
       yaxt = "n", ylab = ".",
       xlim = c(0, 200),
       main = "Posterior for btw-group sigma (3 schools)")
  
  
  
  