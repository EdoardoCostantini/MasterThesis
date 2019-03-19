### Project:     Master Thesis
### Object:      MCMC simulation (fitting the model)
### Description: Contains the functions and specification to fit the bayesian model as defined in
###              model.pdf file in text/.
### Date:        2019-03-15

# Packages
  library(ddpcr)

# Load Sampling functions
  #setwd("/Users/Edoardo/DriveUni/MasterThesis/BayesianThesis")
  source("./R/190313-normalmod-functions.R")
  source("./R/190317-MCMC-functions.R")
  source("./R/190318-priorplot-functions.R")

# Load some data
  RiesbyDat <- read.table("./data/RiesbyDat.txt")
  dat_ALL <- data.frame(cluster = RiesbyDat$id,    #constant dataset w/ naming that fits the loop
                        yvec    = RiesbyDat$depr,  #you only need to change the variables included here
                        xvec    = RiesbyDat$week,  #and everything will be adapted in the loop
                        cvec    = RiesbyDat$endog,
                        inter   = RiesbyDat$inter)
  
# Define psi priors
  B0_pn  <- 1e3*diag(2) # mat-f proper neighbour guess
  B0_iW  <- 1e3*diag(2) # guess of inverse Wishart prior that should ressamble the IG(e,e) from Gelman
  # Rstart <- solve(t(dat_Zi)%*%dat_Zi) this will be defined in the loop because it depedns on Zi
  
  # Educated Guess
    #   intercept <- rep(NA, dat_n)
    #   slope     <- rep(NA, dat_n)
    #   i <- 1
    #   for (id in 1:dat_n) {
    #     #print(paste("Before:", i))
    #     fit <- lm(dat_yvec[i:(i+5)] ~ scale(dat_Xmat[i:(i+5), 2]))
    #     intercept[[id]] <- coef(fit)[1]
    #     slope[[id]]     <- coef(fit)[2]
    #     
    #     i <- i+6
    #     #print(paste("After:", i))
    #   }
    #   vi <- var(intercept) # 20 # these guesses are based on the entire sample
    #   vs <- var(slope)     # 9
    vi <- 20
    vs <- 9
  B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
  
  conds <- c("ALL", "30", "20", "8", "6") # How many clusters to consider
  out1 <- out2 <- out3 <- out4 <- out5 <- vector("list", length = length(conds))
  names(out1) <- names(out2) <- names(out3) <- names(out4) <- names(out5) <- conds
  # Outputs details:
    # Store order: HW, IW, mat-F w/ prior neighbor, mat-f w/ educated guess, mat-f w/ R*
    # each out* will have as many sub lists are there are conditions (length(conds))
  
  set.seed(20190308)
  
# MCMC specs
  MCMC_reps   <- 2000
  MCMC_burnin <- 1/10
  
go <- Sys.time()

for (outps in 1:length(conds)) {
  # Reduce number of clusters
  if(conds[outps] != "ALL"){
    n_goal        <- as.numeric(conds[outps])                       # how many clusters do you want to work with?
    clusters_goal <- sample(unique(dat_ALL$cluster), n_goal)        # sample that many from full dataset.
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal, ]  
  } else {
    dat           <- dat_ALL                                        # First condtion uses all
  }
  
  # Define Data for function
  dat_yvec     <- dat$yvec
  dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
  dat_J        <- length(unique(dat$xvec))
  dat_n        <- nrow(dat)/dat_J
  dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
  dat_subjects <- dat$cluster
  
  # Define prior (only prior element that needs dat_* elements, so in loop)
  Rstar <- solve(t(dat_Zi)%*%dat_Zi)
  
  # Perform MCMC
  print(paste0("n = ", conds[outps],"; prior 1: Huang and Wand Prior") )
  quiet(
    out1[[outps]] <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                                  iniv = 1)
  )
  Sys.time()
  print(paste0("n = ", conds[outps],"; prior 2: inv-Wishart"))
  quiet(
    out2[[outps]] <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                                  iniv = 1,     # the starting values option (1 = lme4 estimates)
                                  B0   = B0_iW) # the prior guess you are using
  )
  Sys.time()
  print(paste0("n = ", conds[outps],"; prior 3: mat-F w/ prior neighbor"))
  quiet(
    out3[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                               iniv = 1,
                               B0   = B0_pn) # prior neighbour
  )
  Sys.time()
  print(paste0("n = ", conds[outps],"; prior 4: mat-f w/ educated guess"))
  quiet(
    out4[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                               iniv = 1,
                               B0   = B0_ed) # B_edug is the informed guess
  )
  Sys.time()
  print(paste0("n = ", conds[outps],"; prior 5: mat-f w/ R*"))
  quiet(
    out5[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                               iniv = 1,
                               B0   = Rstar) # Rstar is the emprical guess
  )
  Sys.time()
}
  
stop <- Sys.time()
timetaken <- stop - go

# Save Results of analysis
  output <- list(out1, out2, out3, out4, out5)
  names(output) <- c("prior_HW", "prior_IW", "prior_MFPN", "prior_MFEG", "prior_MFR*")
  outname <- paste0("./output/", "output", "_rep", MCMC_reps, ".rds")
  
  comment(output) <- paste("time:", round(as.numeric(timetaken), 2), "mins")
  str(output)
  
  saveRDS(output, paste0("./output/", "riesbydata", "_rep", MCMC_reps, ".rds"))
  # loadedresults <- readRDS("./output/riesbydata_rep10000.rds")
  
# Results Exploration
  # Output selection
  which_pri <- 4 # which prior
  which_con <- 1 # which condition
  which_par <- 4 # which object (which parameters type)
  which_col <- 1 # which column in object (which parameter)
  
  output[[which_pri]][[which_con]][[which_par]][,which_col] #out$prior_HW$ALL[[which_par]][,which_col]
  
  # Traceplot
    plot(1:MCMC_reps, output[[which_pri]][[which_con]][[which_par]][,which_col],"l")
  # Estimates
    mean(output[[which_pri]][[which_con]][[which_par]][,which_col])
    median(output[[which_pri]][[which_con]][[which_par]][,which_col])
    # Compare with lme results
    fit
  # Posteriors
    # Histograms
    hist(output[[which_pri]][[which_con]][[which_par]][,which_col], breaks = 100,
         xlab = "Intercept Variance")
    abline(v = mean(output[[which_pri]][[which_con]][[which_par]][,which_col]), col = "blue", lwd = 2)
    abline(v = median(output[[which_pri]][[which_con]][[which_par]][,which_col]), col = "red", lwd = 2)
    
    # Density
    varseq <- seq(0, 100, length = 1000)
    vi_priorG <- c(NA, B0_iW[1, 1], B0_pn[1, 1], B0_ed[1, 1], Rstar[1, 1])
    sdseq <- sqrt(varseq)
    
    plot(density(output[[which_pri]][[which_con]][[which_par]][,which_col]), xlim = c(0, 10), ylim=c(0, 1))
    # Plot prior for variance
      lines(varseq, df_freeb(varseq, nu = 2, d = 1, b = vi_priorG[which_pri]), type = "l", lty = 2)
    # Plot prior for sd
      lines(sdseq, df_freeb_SD(sdseq, nu = 2, d = 1, b = sqrt(vi_priorG[which_pri])), type = "l", lty = 2)
    
    plot(density(output[[which_pri]][[which_con]][[which_par]][,which_col]), xlim = c(0, 10))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# SCRAPS ------------------------------------------------------------------
  
  output <- vector("list", 5)
  output[[1]] <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                       iniv = 1
                       )
  output[[2]] <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                       iniv = 1,
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
  which_r <- 4 # which result (which prior)
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
    # Histograms
    hist(output[[which_r]][[which_p]][,which_c], breaks = 100,
         xlab = "Intercept Variance")
    abline(v = mean(output[[which_r]][[which_p]][,which_c]), col = "blue", lwd = 2)
    abline(v = median(output[[which_r]][[which_p]][,which_c]), col = "red", lwd = 2)
    
    # Density
    plot(density(output[[which_r]][[which_p]][,which_c]), xlim = c(0, max(output[[which_r]][[which_p]][,which_c])+1))
    varseq <- seq(0, 100, length = 1000)
    sdseq <- sqrt(varseq)
    # Plot prior for variance
      lines(varseq, df_freeb(varseq, nu = 2, d = 1, b = vi), type = "l")
    # Plot prior for sd
      lines(sdseq, df_freeb_SD(sdseq, nu = 2, d = 1, b = sqrt(vs)), type = "l", lty = 2)
    
    plot(density(output[[which_r]][[which_p]][,which_c]), xlim = c(0, 10))
    