### Project:     Master Thesis
### Object:      MCMC simulation (fitting the model)
### Description: Contains the functions and specification to fit the bayesian model as defined in
###              model.pdf file in text/.
### Date:        2019-03-15

# Load Sampling functions
  #setwd("/Users/Edoardo/DriveUni/MasterThesis/BayesianThesis")
  source("./R/190313-normalmod-functions.R")
  source("./R/190317-MCMC-functions.R")

# Load some data
  RiesbyDat <- read.table("./data/RiesbyDat.txt")
  
  dat_yvec     <- RiesbyDat$depr
  dat_Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
  dat_J        <- length(unique(RiesbyDat$week))
  dat_n        <- nrow(RiesbyDat)/dat_J
  dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
  dat_subjects <- RiesbyDat$id

# Define psi priors
  B0_pn  <- 1e3*diag(2) # mat-f proper neighbour guess
  B0_iW  <- diag(2) # guess of inverse Wishart prior that shou,d ressamble the IG(e,e) from Gelman
  Rstart <- solve(t(dat_Zi)%*%dat_Zi)
  
  # Educated Guess
    intercept <- rep(NA, dat_n)
    slope     <- rep(NA, dat_n)
    i <- 1
    for (id in 1:dat_n) {
      #print(paste("Before:", i))
      fit <- lm(dat_yvec[i:(i+5)] ~ scale(dat_Xmat[i:(i+5), 2]))
      intercept[[id]] <- coef(fit)[1]
      slope[[id]]     <- coef(fit)[2]
      
      i <- i+6
      #print(paste("After:", i))
    }
    vi <- var(intercept)
    vs <- var(slope)
  B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
    
# Run MCMC
  MCMC_reps   <- 2000
  MCMC_burnin <- 1/10

  output <- vector("list", 5)
  output[[1]] <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                       iniv = 0
                       )
  output[[2]] <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                       iniv = 0,
                       B0   = B0_iW
                       )
  output[[3]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                    iniv = 1,
                    B0   = B0_pn
                    )
  output[[4]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                    iniv = 1,
                    B0   = B0_ed # B_edug is the informed guess
                    )
  output[[5]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                    iniv = 1,
                    B0   = Rstar # Rstar is the emprical guess
                    )
# Results Exploration
  # Output selection
  which_r <- 1 # which result (which prior)
  which_p <- 3 # which object (which parameters type)
  which_c <- 1 # which column in object (which parameter)
  # Traceplot
    plot(1:MCMC_reps, output[[which_r]][[which_p]][,which_c],"l")
  # Estimates
    mean(output[[which_r]][[which_p]][,which_c])
    median(output[[which_r]][[which_p]][,which_c])
    # Compare with lme results
    fit
  # Posteriors
    hist(output[[which_r]][[which_p]][,which_c], breaks = 100,
         xlab = "Intercept Variance")
    abline(v = mean(output[[which_r]][[which_p]][,which_c]), col = "blue", lwd = 2)
    abline(v = median(output[[which_r]][[which_p]][,which_c]), col = "red", lwd = 2)
    plot(df(vi, 2, 1))
    varseq <- seq(0, 40, length=1e4)
    plot(varseq,
         df(varseq, 2, 1),
         col = 2)
    lambdaseq <- seq(0,10,length=1e4)
    plot(lambdaseq,
         dgamma(lambdaseq, 2, 1),
         col = 2)
