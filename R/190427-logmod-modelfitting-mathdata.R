### Project:     Master Thesis
### Object:      Log Model fitting
### Description: Contains the functions and specification to fit the bayesian logistic model
### Date:        2019-04-27

# Packages and functions
  library(ddpcr)
  library(lme4)
  
# Load Data and Conditions ------------------------------------------------
  
  {# Schiz data
  schizdata.noNA <- read.table("./data/schizdata_noNA_dic.txt")
  # Source: Hedeker Gibbons 2006, Longitudinal Data Analysis (p.183)
  #         [original: National Institute of Mental Health Schizophrenia Collaborative Study]
  # Y: Severity of Illness (Schizophrenia),
  #    ordinal (0-7, not_ill-extremely_ill)
  #    dichotomized (< 4 -> 0 ; >= 4 -> 1) (book used only 3 and 4?)
  # X: 0 placebo, 1 drug
  # t: 0, 1, 3, 6 weeks (final)
  # J = 4 repeated measures (weeks: 0, 1, 3, 6)
  # N: 312 individuals (final)
  dat_ALL <- data.frame(cluster = schizdata.noNA$id,
                        yvec    = schizdata.noNA$SevIll,
                        xvec    = schizdata.noNA$week, 
                        cvec    = schizdata.noNA$drug,
                        inter   = schizdata.noNA$inter)
  dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
  dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition
  # Conditions
  allCONDs <- list(
    n_cond = list("312" = unique(dat_ALL$cluster),
                  "8" = c(1129, 1306, 2302, 3106, 1103, 1113, 1109, 2316),
                  "4" = c(1129, 1306, 2316, 1103)),
    J_cond = list("4" = unique(dat_ALL$xvec),
                  "3" = c(1,3,4))
  )
  conds.index <- matrix(c(312, 4,
                          8, 4,
                          8, 3,
                          4, 4,  # not 2 because w/ lme random effects are notidentifiable
                          4, 3), # not smaller than 4 because otherwise:
                        ncol = 2, byrow = T)
  nconds <- nrow(conds.index)}
  
  {# GPA data
  gpadata <- read.table("./data/gpa2_dic.txt"); head(gpadata)
  dat_ALL <- data.frame(cluster = gpadata$cluster,
                        yvec    = gpadata$yvec,
                        xvec    = gpadata$xvec, 
                        cvec    = gpadata$cvec,
                        inter   = gpadata$inter)
  dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
  dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition
  # Conditions  
  allCONDs <- list(
    n_cond = list("200" = unique(dat_ALL$cluster),
                  "8" = c(1, 2, 3, 4, 5, 6, 7, 8),
                  "4" = c(1, 2, 3, 4)),
    J_cond = list("6" = unique(dat_ALL$xvec),
                  "3" = unique(dat_ALL$xvec)[c(1, 3, 4)])
  )
  conds.index <- matrix(c(200, 6,
                          8, 6,
                          8, 3,
                          4, 6,  
                          4, 3), 
                        ncol = 2, byrow = T)
  nconds <- nrow(conds.index)
  }
  
  {# Simulated data
  simdata <- readRDS("./Data/simdata.rds"); head(simdata)
  
  dat_ALL <- data.frame(cluster = simdata$cluster,
                        yvec    = simdata$yvec,
                        xvec    = simdata$xvec, 
                        cvec    = simdata$cvec,
                        inter   = simdata$inter)
  dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
  dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition
  # Conditions  
  allCONDs <- list(
    n_cond = list("400" = unique(dat_ALL$cluster),
                  "8" = c(1, 2, 3, 4, 5, 6, 7, 8),
                  "4" = c(1, 2, 3, 4)),
    J_cond = list("7" = unique(dat_ALL$xvec),
                  "3" = unique(dat_ALL$xvec)[c(1, 3, 4)])
  )
  conds.index <- matrix(c(400, 7,
                          8, 7,
                          8, 3,
                          4, 7,  
                          4, 3), 
                        ncol = 2, byrow = T)
  nconds <- nrow(conds.index)}
  
  # Sampling Functions and specs
  source("./R/190427-logmod-MCMC-functions.R")
  which.model <- "logistic"; source("./R/190330-hyperparameters.R") # which model defines the way the educated guess is ocmputed in the hyperpar script
  
  MCMC_reps   <- 1e4
  MCMC_burnin <- 1/10
  
# LMEr estimates ####
  glme4_loop_COND <- vector("list", length = nconds)
  names(glme4_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
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
    
    # Standard estimation
    glmefit <- glmer(yvec ~ xvec + cvec + inter + (1 + xvec|cluster), family=binomial, data = dat)
    Psi    <- matrix(VarCorr(glmefit)[[1]][1:4], ncol = 2) # as in MulderPericchi2018 code
    Psi_sd <- matrix(c(attributes(VarCorr(glmefit)[[1]])$stddev[1],         # save sd of intercepts
                       rep(attributes(VarCorr(glmefit)[[1]])$correlation[1,2], 2),  # save correlation
                       attributes(VarCorr(glmefit)[[1]])$stddev[2]),         # save sd of slopes
                     ncol = 2)
    sigma2 <- attr(VarCorr(glmefit), "sc")
    bMat   <- as.matrix(ranef(glmefit)$cluster)
    theta  <- fixef(glmefit)
    glme4_loop_COND[[outps]] <- list(glme_fixed  = theta,
                                     glme_bMat   = bMat,
                                     glme_Psi    = Psi,
                                     glme_Psi_sd = Psi_sd,
                                     glme_sigma2 = sigma2)
  }
  str(glme4_loop_COND)
  glme4_loop_COND$`n_200:J_6`$glme_Psi_sd
  
  
# IW SAMPLING ####
set.seed(19044)
IW_loop_COND <- vector("list", length = nconds)
  names(IW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
for (outps in 1:nconds) {
  output_loop_prior <- vector("list", length = length(IW_PR))
    names(output_loop_prior) <- names(IW_PR)
  print(paste("Condition", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
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
        output_loop_prior[[PV]] <- MCMC_mixlogreg_InvWish(yvec    = dat_yvec, 
                                                          Xmat    = dat_Xmat, 
                                                          Zi      = dat_Zi, 
                                                          J       = dat_J, 
                                                          n       = dat_n, 
                                                          samsize = MCMC_reps, 
                                                          burnin  = MCMC_burnin, 
                                                          iniv    = 1,
                                                          B0      = IW_PR[[PV]]), 
        silent = TRUE
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
          output_loop_prior[[PV]] <- MCMC_mixlogreg_matF(yvec = dat_yvec, 
                                                         Xmat = dat_Xmat, 
                                                         Zi = dat_Zi, 
                                                         J = dat_J, 
                                                         n = dat_n, 
                                                         samsize = MCMC_reps, 
                                                         burnin = MCMC_burnin, 
                                                         iniv = 1, 
                                                         B0 = MF_PR[[PV]] ), 
          silent = TRUE
        )
      )
    }
    MF_loop_COND[[outps]] <- output_loop_prior
  }
    str(MF_loop_COND)
  
  # HW SAMPLING ####
  set.seed(19044)
  HW_loop_COND <- vector("list", length = nconds)
    names(HW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
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
    
    quiet(
      try(
        HW_loop_COND[[outps]] <- MCMC_mixlogreg_HW(yvec = dat_yvec, 
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

# Save results
  output <- list(out_lme4 = glme4_loop_COND,
                 out_IW   = IW_loop_COND,
                 out_MF   = MF_loop_COND,
                 out_HW   = HW_loop_COND)
  str(output)
  saveRDS(output, paste0("./output/", "gpadata",Sys.Date(), "-ncond_", nconds, "-rep_", MCMC_reps, ".rds"))
    
# Check results (posterior and traceplots)
  str(output_loop_prior$IW_e1)
  x <- output_loop_prior$IW_e1[[1]]$PD_Psi_sd
  par <- 2
  sum(x[,par] < 0)
  # Posterior
  plot(density(x[,par]))
  median(x[,par])
  hist(x[,par],
       xlim = c(0, 20), breaks = 20)
  hist(x[,par],
       xlim = c(-1, 1), breaks = 20)
  # Traceplot
  x <- out$PD_Psi_sd[,par]
  plot(1:samsize, x,"l",
       xlim = c(0, samsize))
  # Point estimates
  glme4_loop_COND$`n_200:J_6`$glme_Psi_sd
  colMeans(out$PD_Psi_sd)
    