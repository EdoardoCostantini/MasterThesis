### Project:     Master Thesis
### Object:      TEST of posterior draws functions
### Description: Contains the tests for all the posterior functions you defined. Different aspects of the functions
###              (e.g. dimensionality of output) can and should be tested here
### Date:        2019-03-17

context('POSTERIOR DRAWS functions')

#setwd("/Users/Edoardo/DriveUni/MasterThesis/BayesianThesis/test")

# Packages needed
#library(testthat)
#library(lme4)

# Source Functions to be tested
  source("../R/190313-normalmod-functions.R")

# Source Objects needed for test

  # Data

    RiesbyDat <- read.table("../data/RiesbyDat.txt")
      # Make it digestable for the functions
      yvec     <- RiesbyDat$depr
      Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
        xvec    <- Xmat[, 2]
        covar   <- Xmat[, 3]
        inter   <- Xmat[, 4]
      J        <- length(unique(RiesbyDat$week))
      n        <- nrow(RiesbyDat)/J
        cluster <- rep(seq(1, n), each = J)
      Zi       <- cbind(rep(1, length = J), unique(Xmat[, 2]))
      
  # Define priors
      
    bi0   <- rep(0, ncol(Zi))
    avec   <- rep(100,2)
    B0 <- 1e3*diag(2)
    B0Inv <- solve(B0)
    
  # Initial values (a good approximation of the output of an s repetition in MCMC)
    
    fit   <- lmer(yvec ~ xvec + covar + inter + (1 + xvec | cluster), REML = FALSE)
      theta  <- fixef(fit)
      Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
      PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
      sigma2 <- attr(VarCorr(fit), "sc")
      bMat   <- as.matrix(ranef(fit)$cluster)
      Omega = B0

# Define tests of interest
# Test: draw_theta ####

  theta     <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
  tocheck   <- dim(theta)
  benchmark <- c(4, 1) # it should be a vector of length 4

  # test
  test_that('TEST: theta - dimensionality', {
    expect_equal(tocheck, benchmark)
  })

# Test: draw_bMat ####

  bMat_test <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
  tocheck   <- dim(bMat_test)
  benchmark <- c(46, 2)

  # test
  test_that('TEST: bMat - dimensionality', {
    expect_equal(tocheck, benchmark)
  })

# Test: draw_sigam2_IMPprior ####

  sigma2_test    <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
  tocheck   <- length(sigma2_test)
  benchmark <- 1

  # test
  test_that('TEST: sigam2 IMPprior - dimensionality', {
    expect_equal(tocheck, benchmark)
  })

# Test: draw_PsiInv_HW ####

  outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n,nu=2,eta=1/2,Ak=10**3)
    PsiInv_test <- outDummy[[1]]
    avec_test   <- outDummy[[2]]

  tocheck   <- dim(PsiInv_test)
  benchmark <- c(2, 2)

  # test
  test_that('TEST: PsiInv HW - dimensionality', {
    expect_equal(tocheck, benchmark)
  })
  
# Test: draw_PsiInv_InvWish ####

  PsiInv_test <- draw_PsiInv_InvWish(n,bMat,S0=diag(2),e=1)

  tocheck   <- dim(PsiInv_test)
  benchmark <- c(2, 2)

  # test
  test_that('TEST: PsiInv InvWish - dimensionality', {
    expect_equal(tocheck, benchmark)
  })

# Test: draw_PsiInv_matF ####

  outDummy <- draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv,n,d=1,nu=2)
    PsiInv_test <- outDummy[[1]]
    omega_test  <- outDummy[[2]]

  tocheck   <- dim(PsiInv_test)
  benchmark <- c(2, 2)

  # test
  test_that('TEST: PsiInv matF - dimensionality', {
    expect_equal(tocheck, benchmark)
  })