### Project:     Master Thesis
### Object:      TEST of MCMC functions
### Description: Contains the tests for all the posterior functions you defined. Different aspects of the functions
###              (e.g. dimensionality of output) can and should be tested here
### Date:        2019-03-18

context('MCMC functions')

#setwd("/Users/Edoardo/DriveUni/MasterThesis/BayesianThesis/test")

# Packages needed
#library(testthat)
#library(lme4)
#library(ddpcr) # for quiet function

# Source Functions to be tested
  source("../R/190313-normalmod-functions.R")
  source("../R/190317-MCMC-functions.R")

# Source Objects needed for test

  # Data
  RiesbyDat <- read.table("../data/RiesbyDat.txt")
    # Make it digestable for the functions
    yvec     <- RiesbyDat$depr
    Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
    J        <- length(unique(RiesbyDat$week))
    n        <- nrow(RiesbyDat)/J
    Zi       <- cbind(rep(1, length = J), unique(Xmat[, 2]))
    samsize  <- 50
    burnin   <- 1/10

# Define tests of interest
# Test: draw_theta ####
      
# Test Functions
  quiet( out1 <- MCMC_HWprior(yvec = RiesbyDat$depr,
                       Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
                       J    = length(unique(RiesbyDat$week)), 
                       n    = nrow(RiesbyDat)/J,
                       Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
                       iniv = 0,
                       samsize = samsize,
                       burnin = burnin) )
  tocheck <- as.numeric(summary(out1)[, 1])
  benchmark <- c(samsize*ncol(Xmat), # PD_theta  same dim as defined in the MCMC function
                 samsize*n*2,        # PD_bMat   same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # PD_Psi    same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # PD_Psi_sd same dim as defined in the MCMC function
                 samsize*2)          # PD_avec   same dim as defined in the MCMC function

  # test
  test_that('TEST: MCMC_HWprior - dimensionality', {
    expect_equal(tocheck, benchmark)
  })
  
  # Test Functions
  quiet( out1 <- MCMC_matF(yvec = RiesbyDat$depr,
                           Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
                           J    = length(unique(RiesbyDat$week)),
                           n    = nrow(RiesbyDat)/J,
                           Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
                           iniv = 1,
                           B0   = list(nu=2,d=2,e=0,S0=1e3*diag(2)),
                           samsize = samsize,
                           burnin = burnin) )
  tocheck <- as.numeric(summary(out1)[, 1])
  benchmark <- c(samsize*ncol(Xmat), # PD_theta  same dim as defined in the MCMC function
                 samsize*n*2,        # PD_bMat   same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # PD_Psi    same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # PD_Psi_sd same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # Omega
                 4) # PD_Omega   same dim as defined in the MCMC function

  # test
  test_that('TEST: MCMC_matF - dimensionality', {
    expect_equal(tocheck, benchmark)
  })
  
  # Test Functions
  quiet( out1 <- MCMC_invWish(yvec = RiesbyDat$depr,
                       Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
                       J    = length(unique(RiesbyDat$week)),
                       n    = nrow(RiesbyDat)/J,
                       Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
                       iniv = 0,
                       B0   = list(nu=2,e=1,S0=diag(2)),
                       samsize = samsize,
                       burnin = burnin) )
  tocheck <- as.numeric(summary(out1)[, 1])
  benchmark <- c(samsize*ncol(Xmat), # PD_theta  same dim as defined in the MCMC function
                 samsize*n*2,        # PD_bMat   same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # PD_Psi    same dim as defined in the MCMC function
                 samsize*ncol(Zi)*2, # PD_Psi_sd same dim as defined in the MCMC function
                 3) # number of objects in a list prior (nu (useless), e, S0)

  # test
  test_that('TEST: MCMC_invWish - dimensionality', {
    expect_equal(tocheck, benchmark)
  })
  
  