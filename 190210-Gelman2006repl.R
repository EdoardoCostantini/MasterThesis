# Replicate findings in Gelman 2006
# Trying to replicate what was done in Gelman2006

# Sources: Gelman2006
# Set up
  library(R2OpenBUGS)
  # Set working dir to project location (containing BugsModel folder)
  wd <- getwd()
  # Set the Wine working directory and the directory to OpenBUGS
  WINE         = "/usr/local/bin/wine"
  WINEPATH     = "/usr/local/bin/winepath"
  OpenBUGS.pgm = "/Applications/OpenBUGS323/OpenBUGS.exe"

# 5.1 Noninformative prior distribution 8-schools problem ---------------------

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
  model_Gelman2006_UNIF <- paste0(wd, "/BugsModels/Gelman2006-5_1-UNIF.txt")
  # model files should be in a BugsModels directory in your wd
  # Initial values (function)
  inits <- function(){
    list (theta       = rnorm(J, 0, 100), 
          mu.theta    = rnorm(1, 0, 100),
          sigma.theta = runif(1, 0, 100))
  }
  # Simulation
  schools.sim.UNI <- bugs (data, inits, parameters,
                           model.file = model_Gelman2006_UNIF,
                           n.chains = 3, n.iter = 6000, debug = FALSE,
                           OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T) # very imp on mac
    
#> Model2: IG(1,1) Prior ####
  model_Gelman2006_INVG11 <- paste0(wd, "/BugsModels/Gelman2006-5_1-invG11.txt")
  inits <- function (){
    list (theta     = rnorm(J, 0, 100), 
          mu.theta  = rnorm(1, 0, 100),
          tau.theta = runif(1, 0, 100)) # tau is the precision (which is 1/sigma)
                                        # this is different beacuse of the different prior
                                        # (and its parametrization)
  }
  schools.sim.GP <- bugs (data, inits, parameters,
                          model.file = model_Gelman2006_INVG11,
                          n.chains = 3, n.iter = 6000,
                          OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T)
  print(schools.sim.GP)
  #plot(schools.sim.GP)
  
#> Model3: IG(.001, .001) Prior ####
  model_Gelman2006_INVG00 <- paste0(wd, "/BugsModels/Gelman2006-5_1-invG001001.txt")
  inits <- function (){
    list (theta     = rnorm(J, 0, 100), 
          mu.theta  = rnorm(1, 0, 100),
          tau.theta = runif(1, 0, 100)) # tau is the precision (which is 1/sigma)
                                        # this is different beacuse of the different prior
                                        # (and its parametrization)
  }
  schools.sim.GP2 <- bugs (data, inits, parameters,
                          model.file = model_Gelman2006_INVG00,
                          n.chains = 3, n.iter = 6000,
                          OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T)

#> Results and plot ####
  # Results
    print(schools.sim.UNI)
    #plot(schools.sim.UNI)
    print(schools.sim.GP)
    #print(schools.sim.GP)
    print(schools.sim.GP2)
    #plot(schools.sim.GP2)
  # Mean, median and bias
    # Mean
    post_mean <- apply(cbind(schools.sim.UNI$sims.list$sigma.theta, 
                             schools.sim.GP$sims.list$sigma.theta,
                             schools.sim.GP2$sims.list$sigma.theta),
                       2, mean)
    post_median <- apply(cbind(schools.sim.UNI$sims.list$sigma.theta, 
                               schools.sim.GP$sims.list$sigma.theta,
                               schools.sim.GP2$sims.list$sigma.theta),
                         2, median)
    # To me it seems that the median is actually worst with the inv-gamma than with the 
    # mean as point estimate. In Browne and Draper it is found that the relative bias
    # is higehr with the mean than with the median as point estimate, and that it is
    # higher when using the uniform than when using the inv-gamma (see table 4, p16).
    # - Same model is used in both papers except sigma.y considered as known in Gelman
    #   (in Browne same prior is used for sigma.y and sigma.u, sigma.y is the error variance)
    # - Same priors in both paper (e = .001 in both at least for the inv-gamma)
    #   uniform prior in Gelman is U(0, 100) (I think), instead of U(0, 1000) (as in Browne).
    
  # Posterior Plots
    par(mfcol=c(3,1))
  # Uniform
  hist(schools.sim.UNI$sims.list$sigma.theta,
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "sigma.theta",
       main = "8 Schools: posterior on sigma_a given uniform prior on sigma_a",
       breaks = 40)
  # IG(1, 1)
  hist(schools.sim.GP$sims.list$sigma.theta,
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "sigma.theta",
       main = "8 Schools: posterior on sigma_a given IG(1,1) prior on sigma_a**2",
       breaks = 50)
  # IG(.001,.001)
  hist(schools.sim.GP2$sims.list$sigma.theta,
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "sigma.theta",
       main = "8 Schools: posterior on sigma_a given IG(.001,.001) prior on sigma_a**2",
       breaks = 50)
  # Conclusions:
  # 1) Uniform Prior density supports for range of values for sigma.theta 
  #    (between group variance) below 20
  # 2) Compared to the gamma priors specified, it seems that the uniform
  #    distirubtion is closer to "noninformative" as it appears to constrain 
  #    the posterior inference less.
  # 3) The IG (e,e) is not at all "noninformative" for this problem since the 
  #    resulting posterior is higly sensity to the choice of e

# 5.2 Weakly informative prior distribution 3-schools -------------------------
#> Prep ####
  # DATA
    data("schools")
    J <- 3                # 3 schools!
    y <- schools$estimate
    sigma.y <- schools$sd
  # Prior hyperparameter for half-t distribution
    prior.scale <- 25 # scale parameter for the half-Cauchy distirbution (see text for decision)
  # Parameters definition
    parameters <- c ("theta", "mu.theta", "sigma.theta")
  data <- list ("J", "y", "sigma.y", "prior.scale")
  
#> Model 1: Uniform prior ####
  model_Gelman2006_UNIF <- paste0(wd, "/BugsModels/Gelman2006-5_1-UNIF.txt")
  inits_UNIF <- function(){
    list (theta       = rnorm(J, 0, 100), 
          mu.theta    = rnorm(1, 0, 100),
          sigma.theta = runif(1,0,100))
  }
  schools.sim.UNI3 <- bugs(data, inits_UNIF, parameters,
                           model.file = model_Gelman2006_UNIF,
                           n.chains = 3, n.iter = 6000,
                           OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T) # very imp on mac
  print(schools.sim.UNI3)
  
#> Model 2: Half Cauchy prior (half-t w/ df nu = 1) ####
  model_Gelman2006_HC <- paste0(wd, "/BugsModels/Gelman2006-5_2-halfCauchyPrior.txt")
  inits_HC <- function (){
    list (eta      = rnorm(J), 
          mu.theta = rnorm(1), 
          xi       = rnorm(1), 
          tau.eta  = runif(1))}
  schools.sim.HC <- bugs (data, inits_HC, parameters,
                       model.file = model_Gelman2006_HC,
                       n.chains = 3, n.iter = 6000,
                       OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T) # very imp on mac
  print(schools.sim.HC)
  #plot(schools.sim.HC)
  attributes(schools.sim.HC$sims.list)
  
#> Results and plot ####
  # Plots
  par(mfcol=c(3,1))
  # Uniform Prior on 8 schools problem
  hist(schools.sim.UNI$sims.list$sigma.theta, 
       yaxt = "n", ylab = ".", xlim = c(0, 200), xlab = "sigma.theta",
       main = "8 Schools: posterior on sigma_a given uniform prior on sigma_a",
       breaks = 20)
  # Uniform Prior on 3 schools problem (BAD)
  hist(schools.sim.UNI3$sims.list$sigma.theta,
       yaxt = "n", ylab = ".", xlim = c(0, 200), xlab = "sigma.theta",
       main = "3 Schools: posterior on sigma_a given uniform prior on sigma_a",
       breaks = 100)
  # Half Cauchy prior (25)
  hist(schools.sim.HC$sims.list$sigma.theta,
       yaxt = "n", ylab = ".", xlim = c(0, 200), xlab = "sigma.theta",
       main = "3 Schools: posterior on sigma_a given Half_Cauchy (25) prior on sigma_a",
       breaks = 40)
  
  # Conclusions:
  # 1. The uniform prior disitbrution works well only when there are more than 3
  #    groups.
  # 2. The half cauchy prior works well even when the number of groups is small
  
  