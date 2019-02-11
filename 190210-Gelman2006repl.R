# Replicate findings in Gelman 2006
# Trying to replicate what was done in Gelman2006

# Sources: Gelman2006
# Set up
  library(R2OpenBUGS)
  # Set the Wine working directory and the directory to OpenBUGS
  WINE = "/usr/local/bin/wine"
  WINEPATH = "/usr/local/bin/winepath"
  OpenBUGS.pgm = "/Applications/OpenBUGS323/OpenBUGS.exe"

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
  