### Project:     Master Thesis
### Object:      MCMC simulation (fitting the model)
### Description: Contains the functions and specification to fit the bayesian model as defined in
###              model.pdf file in text/.
### Date:        2019-03-15

# Packages
  library(ddpcr)
  library(parallel)

# Load Data ---------------------------------------------------------------

  # Riesbt dat
  RiesbyDat <- read.table("./data/RiesbyDat.txt")
  dat_ALL <- data.frame(cluster = RiesbyDat$id,    #constant dataset w/ naming that fits the loop
                        yvec    = RiesbyDat$depr,  #you only need to change the variables included here
                        xvec    = RiesbyDat$week,  #and everything will be adapted in the loop
                        cvec    = RiesbyDat$endog,
                        inter   = RiesbyDat$inter)
  dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
  dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition

# Load Sampling functions
  #setwd("/Users/Edoardo/DriveUni/MasterThesis/BayesianThesis")
  source("./R/190313-normalmod-functions.R")
  source("./R/190317-MCMC-functions.R")
  which.model <- "normal_rep"; source("./R/190330-hyperparameters.R")
    # priors are defined in this R script for ease
    # by sourcing this code, you get a list for each type of prior
    # (ie inverse-Wishart, matrix-F, HW does not need it)
    # where the objects contained in the list are the actual priors
    # usually to be included in the "B0 = " argument of your sampling functions

  
# Set up  ---------------------------------------------------------------------

  allCONDs <- list(
    n_cond = list("46" = unique(dat_ALL$cluster),
                  "8" = c(101, 117, 505, 302, 335, 338, 319, 353), # 319 was 350
                  "4" = c(101, 117, 505, 302)), #504, 319, 328, 353))
    J_cond = list("6" = unique(dat_ALL$xvec),
                  "4" = c(1,2,3,4),
                  "3" = c(1,2,3))
  )
  
  conds.index <- matrix(c(46, 6,
                          8, 6,
                          8, 4,
                          8, 3,  # not 2 because w/ lme random effects are notidentifiable
                          4, 6,
                          4, 4,
                          4, 3), # not smaller than 4 because otherwise:
                        ncol = 2, byrow = T)
  
  nconds <- nrow(conds.index)
  
#####################################################
  
  # LMEr estimates ####
  lme4_loop_COND <- vector("list", length = nconds)
  names(lme4_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec     <- dat$yvec
    dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
    dat_J        <- length(unique(dat$xvec))
    dat_n        <- nrow(dat)/dat_J
    dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
    dat_subjects <- dat$cluster
    
    # # Standard estimation
    lmefit <- lmer(yvec ~ xvec + cvec + inter + (1 + xvec | cluster), data = dat, REML = FALSE)
    Psi    <- matrix(VarCorr(lmefit)[[1]][1:4], ncol = 2) # as in MulderPericchi2018 code
    Psi_sd <- matrix(c(attributes(VarCorr(lmefit)[[1]])$stddev[1],         # save sd of intercepts
                       rep(attributes(VarCorr(lmefit)[[1]])$correlation[1,2], 2),  # save correlation
                       attributes(VarCorr(lmefit)[[1]])$stddev[2]),         # save sd of slopes
                     ncol = 2)
    sigma2 <- attr(VarCorr(lmefit), "sc")
    bMat   <- as.matrix(ranef(lmefit)$cluster)
    theta  <- fixef(lmefit)
    lme4_loop_COND[[outps]] <- list(lme_fixed  = theta,
                                    lme_bMat   = bMat,
                                    lme_Psi    = Psi,
                                    lme_Psi_sd = Psi_sd,
                                    lme_sigma2 = sigma2)
  }
  str(lme4_loop_COND)
  
  # Sampling Repetitions
  MCMC_reps   <- 1e4
  MCMC_burnin <- 1/10
  
  # IW SAMPLING ####
  set.seed(19044)
  IW_loop_COND <- vector("list", length = nconds)
    names(IW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  for (outps in 1:nconds) {
    output_loop_prior <- vector("list", length = length(IW_PR))
      names(output_loop_prior) <- names(IW_PR)
    #outps <- 2
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec     <- dat$yvec
    dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
    dat_J        <- length(unique(dat$xvec))
    dat_n        <- nrow(dat)/dat_J
    dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
    dat_subjects <- dat$cluster
    
    for (PV in 1:length(IW_PR)) {
      quiet(
        try(
          output_loop_prior[[PV]] <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1,
                                                  B0 = IW_PR[[PV]]), silent = TRUE
        )
      )
    }
    IW_loop_COND[[outps]] <- output_loop_prior
  }
    str(IW_loop_COND)
  
  # MF SAMPLING ####
  set.seed(19044)
  MF_loop_COND <- vector("list", length = nconds)
    names(MF_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  for (outps in 1:nconds) {
    #outps <- 2
    output_loop_prior <- vector("list", length = length(MF_PR))
      names(output_loop_prior) <- names(MF_PR)
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec     <- dat$yvec
    dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
    dat_J        <- length(unique(dat$xvec))
    dat_n        <- nrow(dat)/dat_J
    dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
    dat_subjects <- dat$cluster
    for (PV in 1:length(MF_PR)) {
        quiet(
      try(
          output_loop_prior[[PV]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1, 
                                               B0 = MF_PR[[PV]] ), silent = TRUE
        )
          )
    }
    MF_loop_COND[[outps]] <- output_loop_prior
  }
    str(MF_loop_COND)
    # 
    # par <- 1
    # x <- MF_loop_COND$`n_4:J_3`[[4]]$PD_Psi_sd[,c(1,4)]
    #       # Posteriro Plot
    #       plot(density(x[,par]),
    #            xlim = c(0, 20), ylim = c(0, .5))
    #       
    # plot(1:5000, x[,par],"l",
    #    xlim = c(0, 5000), ylim = c(0, 50))
    
    
  # HW SAMPLING ####
  set.seed(19044)
  HW_loop_COND <- vector("list", length = nconds)
    names(HW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  for (outps in 1:nconds) {
    #outps <- 2
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec     <- dat$yvec
    dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
    dat_J        <- length(unique(dat$xvec))
    dat_n        <- nrow(dat)/dat_J
    dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
    dat_subjects <- dat$cluster
    
    quiet(
      try(
        HW_loop_COND[[outps]] <- MCMC_HWprior(yvec = dat_yvec, 
                                              Xmat = dat_Xmat, 
                                              Zi = dat_Zi, 
                                              J = dat_J, 
                                              n = dat_n, 
                                              samsize = MCMC_reps, 
                                              burnin = MCMC_burnin, 
                                              iniv = 1), 
        silent = TRUE
      )
    )
     
  }
    str(HW_loop_COND)
    
    
#####################################################

  # Save output
  output <- list(out_lme4 = lme4_loop_COND,
                 out_IW  = IW_loop_COND,
                 out_MF  = MF_loop_COND,
                 out_HW  = HW_loop_COND)
  saveRDS(output, paste0("./output/", "riesbydata",Sys.Date(), "-ncond_", nconds, "-rep_", MCMC_reps, ".rds"))
  
  # Checks and stitching
  plot(density(MF_loop_COND$`n_4:J_3`$MF_e.1$PD_Psi_sd[,1]))
  
  output_final <- readRDS("./output/riesbydata2019-04-09-ncond_7-rep_10000.rds")
  str(output_final$out_lme4)
  str(lme4_loop_COND)
  str(output_final$out_IW)
  str(IW_loop_COND)
  str(output_final$out_HW)
  str(HW_loop_COND)
  str(output_final$out_MF)
  str(MF_loop_COND)
  
  # Join Results
  output_final$out_lme4 <- lme4_loop_COND
  output_final$out_IW   <- IW_loop_COND
  output_final$out_MF   <- MF_loop_COND
  output_final$out_HW   <- HW_loop_COND
  
  output <- output_final
  saveRDS(output, paste0("./output/", "riesbydata",Sys.Date(), "-ncond_", nconds, "-rep_", MCMC_reps,"DEFINITIVE", ".rds"))
  
  # SAVE
  
  str(output$out_IW)

# mcapply version ---------------------------------------------------------

  nconds <- nrow(conds.index)
  output_temp <- vector("list", nconds)
  names(output_temp) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  out_lme4 <- out_IW <- out_HW <- out_MF <- output_temp
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(190401)
  
begtime <- Sys.time()
# LMEr estimates
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec     <- dat$yvec
    dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
    dat_J        <- length(unique(dat$xvec))
    dat_n        <- nrow(dat)/dat_J
    dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
    dat_subjects <- dat$cluster

    # Sampling Repetitions
    MCMC_reps   <- 5e3
    MCMC_burnin <- 1/10

    # # Standard estimation
    lmefit <- lmer(yvec ~ xvec + cvec + inter + (1 + xvec | cluster), data = dat, REML = TRUE) # with small samples regular ML produces biased estiamtes, therefore advisable to use REML (see Goldstein andn Browene Draper summary of Goldstein)
    Psi    <- matrix(VarCorr(lmefit)[[1]][1:4], ncol = 2) # as in MulderPericchi2018 code
    Psi_sd <- matrix(c(attributes(VarCorr(lmefit)[[1]])$stddev[1],         # save sd of intercepts
                       rep(attributes(VarCorr(lmefit)[[1]])$correlation[1,2], 2),  # save correlation
                       attributes(VarCorr(lmefit)[[1]])$stddev[2]),         # save sd of slopes
                     ncol = 2)
    sigma2 <- attr(VarCorr(lmefit), "sc")
    bMat   <- as.matrix(ranef(lmefit)$cluster)
    theta  <- fixef(lmefit)
    out_lme4[[outps]] <- list(lme_fixed  = theta,
                              lme_bMat   = bMat,
                              lme_Psi    = Psi,
                              lme_Psi_sd = Psi_sd,
                              lme_sigma2 = sigma2)
    # #IW - fitting the model using the IW prior (will be done for all the prior specifications in the IW_PR list)
    # out_IW[[outps]] <- mcmapply(MCMC_invWish, IW_PR,
    #                    MoreArgs = list(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1),
    #                    SIMPLIFY = FALSE)
    # #HW
    # quiet(
    #   out_HW[[outps]] <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
    #                          iniv = 1)
    # )
    #Mat-f - fitting the model using the maf-f prior (will be done for all the prior specifications in the MF_PR list)
    # out_MF[[outps]] <- mcmapply(MCMC_matF, MF_PR,
    #                    MoreArgs = list(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1),
    #                    SIMPLIFY = FALSE)
  }

out <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1, B0 = MF_PR[[3]] )

colMeans(out_MF$`n_4:J_6`[[1]]$PD_Psi_sd)
endtime <- Sys.time()
begtime-endtime

  str(out_lme4)
    out_lme4$n8_J4$lme_Psi_sd
  str(out_IW)
    out_IW$n46_J6$IW_e1$PD_Psi_sd[, 1]
    out_IW$n4_J4$IW_e1$PD_Psi
  str(out_HW)
    out_HW$n46_J6$PD_Psi_sd[, 1]
    out_HW$n20_J3$PD_Psi_sd
  str(out_MF)
    out_MF$n46_J6$MF_e1$PD_Psi_sd[, 1]
    out_MF$n8_J6$MF_pn$PD_Psi_sd
    out_MF$n8_J3$MF_e0.1$hyperpar
# - for IW figure out why the prior you thought was uninformative, is actually more informative than the other one
# - e1 should be ok, why not? See again how this prior should be the regular nu = 2, d = 1 and S0 = B0eg
# Maybe same thing goes wrong w/ last prior
# Run again after you figure these details out
  plot(1:MCMC_reps, out_MF$n46[[1]]$PD_Psi_sd[,1],"l")

# Save Results of analysis ------------------------------------------------

  output_final <- list(out_lme4=out_lme4,out_IW=out_IW,out_HW=out_HW,out_MF=out_MF)
  str(output_final)
  saveRDS(output_final, paste0("./output/", "riesbydata",Sys.Date(), "-ncond_", nconds, "-rep_", MCMC_reps, ".rds"))
  output_final <- readRDS("./output/riesbydata2019-04-04-ncond_7-rep_5000.rds")
  
  saveRDS(output_final, paste0("./output/", "TEMP190402", ".rds"))
  OLDoutput <- readRDS(paste0("./output/", "TEMP190402", ".rds"))
  
  str(OLDoutput$out_IW)
    
  # Check consistency with resulst of priors with old fitting method Previous results of prior
  # Fit with invWishart Prior
  # CHECKed: to obtain same prior as in old fit use the following prior list(e=1,S0=diag(2))
  out_IW <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1,
               B0 = list(e=1,S0=diag(2)))
  par(mfcol = c(2, 5))
  plot(density(out_IW$`n_46:J_6`$IW_e.01$PD_Psi_sd[,1]),xlim = c(0, 20), ylim = c(0,.6))
  plot(density(out_IW$PD_Psi_sd[,4]),xlim = c(0, 20), ylim = c(0,.6))
  plot(density(out_IW$PD_Psi_sd[,2]),xlim = c(-1, 1), ylim = c(0,2))
  
  # Fit with HW prior
  # CHECKed: same results as with previous set up
        quiet(
  out_HW <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                         iniv = 1)
        )
  str(out_HW)
  
  par(mfcol = c(2, 5))
  plot(density(out_HW$PD_Psi_sd[,1]),xlim = c(0, 20), ylim = c(0,2))
  plot(density(out_HW$PD_Psi_sd[,4]),xlim = c(0, 20), ylim = c(0,2))
  plot(density(out_HW$PD_Psi_sd[,2]),xlim = c(-1, 1), ylim = c(0,2))
  
  # Fit with mat-f prior
  # CHECKed: to obtain same priors as in old fit use list(nu=2,d=1,e=0,S0=1e3*diag(2)) (or other prior guess for S0)
  outPN <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 0,
            B0 = list(nu=2,d=1,e=0,S0=1e3*diag(2)))
  outPN <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 0,
            B0 = list(nu=2,d=1,e=0,S0=matrix(c(20, 0, 0, 9), ncol = 2)))
  par(mfcol = c(2, 5))
  plot(density(out_MF$n46[[1]]$PD_Psi_sd[,1]),xlim = c(0, 20), ylim = c(0,1))
  plot(density(outPN$PD_Psi_sd[,4]),xlim = c(0, 20), ylim = c(0,.3))
  plot(density(outPN$PD_Psi_sd[,2]),xlim = c(-1, 1), ylim = c(0,2))
  
  out_MF <- mcmapply(MCMC_matF, MF_PR,
                     MoreArgs = list(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 0),
                     SIMPLIFY = FALSE)
  str(out_MF)
  
  
  
# OLD ---------------------------------------------------------------------


  
# Define psi priors
  B0_pn  <- 1e3*diag(2) # mat-f proper neighbour guess
  B0_iW  <- diag(2)     # guess of inverse Wishart prior that should ressamble the IG(e,e) from Gelman
                        # for now follow Gelman 2014 (p.73) indication of non-informative choice
  Rstar <- solve(t(dat_Zi)%*%dat_Zi) #this will be defined in the loop because it depedns on Zi
  
  # Educated Guess
      # intercept <- rep(NA, dat_ALLn)
      # slope     <- rep(NA, dat_ALLn)
      # i <- 1
      # for (id in 1:dat_ALLn) {
      #   #print(paste("Before:", i))
      #   fit <- lm(dat_ALL$yvec[i:(i+5)] ~ scale(dat_ALL$xvec[i:(i+5)]))
      #   intercept[[id]] <- coef(fit)[1]
      #   slope[[id]]     <- coef(fit)[2]
      # 
      #   i <- i+6
      #   #print(paste("After:", i))
      # }
      # vi <- var(intercept) # 20 # these guesses are based on the entire sample
      # vs <- var(slope)     # 9
    vi <- 20
    vs <- 9
  B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
  
  set.seed(20190308)
  
  conds <- c(dat_ALLn, 30, 20, 8, 4) # How many clusters to consider
  clusters4cond  <- list(unique(dat_ALL$cluster),
                         sample(unique(dat_ALL$cluster), conds[2]),
                         sample(unique(dat_ALL$cluster), conds[3]),
                         c(101, 117, 505, 302, 335, 338, 350, 353),
                         c(504, 319, 328, 337))
  
  out1 <- out2 <- out3 <- out4 <- out5 <- vector("list", length = length(conds))
  names(out1) <- names(out2) <- names(out3) <- names(out4) <- names(out5) <- conds
  # Outputs details:
    # Store order: HW, IW, mat-F w/ prior neighbor, mat-f w/ educated guess, mat-f w/ R*
    # each out* will have as many sub lists are there are conditions (length(conds))
  
# MCMC specs
  MCMC_reps   <- 2000
  MCMC_burnin <- 1/10
  
# Fit model ---------------------------------------------------------------
  
go <- Sys.time()

for (outps in 2:length(conds)) {
  
  clusters_goal <- clusters4cond[[outps]]
  dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal, ]  
  
  # Define Data for function
  dat_yvec     <- dat$yvec
  dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
  dat_J        <- length(unique(dat$xvec))
  dat_n        <- nrow(dat)/dat_J
  dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
  dat_subjects <- dat$cluster
  
  print(paste0("n = ", conds[outps],"; prior 1: mat-F w/ prior neighbor"))
  quiet(
    out1[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1,
                               B0   = MF_PR$MF_e1)
  )
  Sys.time()
  print(paste0("n = ", conds[outps],"; prior 2: mat-f w/ educated guess"))
  quiet(
    out2[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                               iniv = 1,
                               B0   = MF_PR$MF_e0.5) # B_edug is the informed guess
  )
  Sys.time()
  print(paste0("n = ", conds[outps],"; prior 3: mat-f w/ R*"))
  quiet(
    out3[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                               iniv = 1,
                               B0   = MF_PR$MF_e0.1) # Rstar is the emprical guess
  )
  print(paste0("n = ", conds[outps],"; prior 4: mat-f w/ R*"))
  quiet(
    out4[[outps]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                               iniv = 1,
                               B0   = MF_PR$MF_e0.01) # Rstar is the emprical guess
  )
  Sys.time()
}
  
stop <- Sys.time()
timetaken <- stop - go

str(out2)
str(out_MF)

# Save Results of analysis ------------------------------------------------

  pg <- data.frame(vi_pg = c(10**5, B0_iW[1, 1], B0_pn[1, 1], B0_ed[1, 1], Rstar[1, 1]), # list of prior gruesses used
                   vs_pg = c(10**5, B0_iW[2, 2], B0_pn[2, 2], B0_ed[2, 2], Rstar[2, 2]))
  row.names(pg) <- c("prior_HW", "prior_IW", "prior_MFPN", "prior_MFEG", "prior_MFR*")
  output <- list(out1, out2, out3, out4, out5, pg)
  names(output) <- c("prior_HW", "prior_IW", "prior_MFPN", "prior_MFEG", "prior_MFR*", "prior_guesses")

  comment(output) <- paste("time:", round(as.numeric(timetaken), 2), "mins")
  ls(output)
  
  saveRDS(output, paste0("./output/", "riesbydata","-cond_", paste(conds,collapse="_"), "-rep_", MCMC_reps,"-HWnu0_1", ".rds"))
  #output <- readRDS("./output/riesbydata_rep10000.rds")
  
# Results Exploration -----------------------------------------------------
  fit   <- lmer(yvec ~ xvec + cvec + inter + (1 + xvec | cluster), REML = FALSE, data = dat_ALL)
  # Output selection
  which_pri <- 1 # which prior
  which_con <- 1 # which condition
  which_par <- 4 # which object (which parameters type)
  which_col <- 1 # which column in object (which parameter)
  
  #output[[which_pri]][[which_con]][[which_par]][,which_col] #out$prior_HW$ALL[[which_par]][,which_col]
  
  # Traceplot
    plot(1:MCMC_reps, output[[which_pri]][[which_con]][[which_par]][,which_col],"l")
  # Estimates
    mean(output[[which_pri]][[which_con]][[which_par]][,which_col])
    median(output[[which_pri]][[which_con]][[which_par]][,which_col])
  #   # Compare with lme results
  #   fit
  # # Posteriors
  #   # Histograms
  #   hist(output[[which_pri]][[which_con]][[which_par]][,which_col], breaks = 100,
  #        xlab = "Intercept Variance")
  #   abline(v = mean(output[[which_pri]][[which_con]][[which_par]][,which_col]), col = "blue", lwd = 2)
  #   abline(v = median(output[[which_pri]][[which_con]][[which_par]][,which_col]), col = "red", lwd = 2)
  #   
    # Density
    varseq <- seq(0, 10000, length = 1000)
    vi_priorG <- c(NA, B0_iW[1, 1], B0_pn[1, 1], B0_ed[1, 1], Rstar[1, 1])
    sdseq <- sqrt(varseq)

    plot(density(output[[which_pri]][[which_con]][[which_par]][,which_col]), xlim = c(0, 100), ylim=c(0, 1))
    # Plot prior for variance
      #lines(varseq, df_freeb(varseq, nu = 2, d = 1, b = vi_priorG[which_pri]), type = "l", lty = 2)
    # Plot prior for sd
      lines(sdseq, df_freeb_SD(sdseq, nu = 2, d = 1, b = sqrt(vi_priorG[which_pri])), type = "l", lty = 2)
    
    plot(density(output[[which_pri]][[which_con]][[which_par]][,which_col]), xlim = c(0, 10))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# SCRAPS ------------------------------------------------------------------
  # # Schiz dat
  # schizdata <- read.table("./data/schizdata_noNA.txt")
  #   head(schizdata)
  # dat_ALL <- data.frame(cluster = schizdata$id,     #constant dataset w/ naming that fits the loop
  #                       yvec    = schizdata$SevIll, #you only need to change the variables included here
  #                       xvec    = schizdata$week,   #and everything will be adapted in the loop
  #                       cvec    = schizdata$drug,
  #                       inter   = schizdata$inter)
  # dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
  # dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition
  
  # system.time(
  #   outMcap <- mcmapply(MCMC_invWish, IW_PRIORS,
  #                       MoreArgs = list(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1),
  #                       SIMPLIFY = FALSE)
  # )
  # str(outMcap)
  # 
  # # With a loop
  # system.time(
  #   for (i in 1:2) {
  #     outLoop[[i]] <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
  #                                  iniv = 1,     # the starting values option (1 = lme4 estimates)
  #                                  B0   = B0_iW[[i]]) # the prior guess you are using
  #   }
  # )
  # str(outLoop)
  
  # Define Priors for Psi
  # B0_pn  <- 1e3*diag(2) # mat-f proper neighbour guess
  # B0_iW  <- diag(2)     # guess of inverse Wishart prior that should ressamble the IG(e,e) from Gelman
  #                       # for now follow Gelman 2014 (p.73) indication of non-informative choice
  # Rstar <- solve(t(dat_Zi)%*%dat_Zi) #this will be defined in the loop because it depedns on Zi
  # B0_ed  <- matrix(c(20, 0, 0, 9), ncol = 2) # guess based on data exploration and knowledge


    
    
    
#   output <- vector("list", 5)
#   output[[1]] <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
#                        iniv = 1
#                        )
#   output[[2]] <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
#                        iniv = 1,
#                        B0   = B0_iW
#                        )
#   output[[3]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
#                     iniv = 1,
#                     B0   = B0_pn
#                     )
#   output[[4]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
#                     iniv = 1,
#                     B0   = B0_ed # B_edug is the informed guess
#                     )
#   output[[5]] <- MCMC_matF(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
#                     iniv = 1,
#                     B0   = Rstar # Rstar is the emprical guess
#                     )
# # Results Exploration
#   # Output selection
#   which_r <- 4 # which result (which prior)
#   which_p <- 3 # which object (which parameters type)
#   which_c <- 1 # which column in object (which parameter)
#   # Traceplot
#     plot(1:MCMC_reps, output[[which_r]][[which_p]][,which_c],"l")
#   # Estimates
#     mean(output[[which_r]][[which_p]][,which_c])
#     median(output[[which_r]][[which_p]][,which_c])
#     # Compare with lme results
#     fit
#   # Posteriors
#     # Histograms
#     hist(output[[which_r]][[which_p]][,which_c], breaks = 100,
#          xlab = "Intercept Variance")
#     abline(v = mean(output[[which_r]][[which_p]][,which_c]), col = "blue", lwd = 2)
#     abline(v = median(output[[which_r]][[which_p]][,which_c]), col = "red", lwd = 2)
#     
#     # Density
#     plot(density(output[[which_r]][[which_p]][,which_c]), xlim = c(0, max(output[[which_r]][[which_p]][,which_c])+1))
#     varseq <- seq(0, 100, length = 1000)
#     sdseq <- sqrt(varseq)
#     # Plot prior for variance
#       lines(varseq, df_freeb(varseq, nu = 2, d = 1, b = vi), type = "l")
#     # Plot prior for sd
#       lines(sdseq, df_freeb_SD(sdseq, nu = 2, d = 1, b = sqrt(vs)), type = "l", lty = 2)
#     
#     plot(density(output[[which_r]][[which_p]][,which_c]), xlim = c(0, 10))
#     