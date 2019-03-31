### Project:     Master Thesis
### Object:      MCMC simulation (fitting the model)
### Description: Contains the functions and specification to fit the bayesian model as defined in
###              model.pdf file in text/.
### Date:        2019-03-15

# Packages
  library(ddpcr)
  library(parallel)

# Load data
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
  source("./R/190330-hyperparameters.R")
    # priors are defined in this R script for ease
    # by sourcing this code, you get a list for each type of prior
    # (ie inverse-Wishart, matrix-F, HW does not need it)
    # where the objects contained in the list are the actual priors
    # usually to be included in the "B0 = " argument of your sampling functions
  source("./R/190318-priorplot-functions.R")

  
# NEW ---------------------------------------------------------------------
  
# Conditions
  n_cond <- c(46, 20, 8, 4)
  J_cond <- c(6, 4, 3)
  conds <- expand.grid(n_cond, J_cond)[order(expand.grid(n_cond, J_cond)[,1], decreasing = T),][-c(2,3), ]
  #conds <- expand.grid(n_cond, J_cond)[order(expand.grid(n_cond, J_cond)[,1], decreasing = T),][c(1,4), ]
  
# Storing Objects
  output <- vector("list", nrow(conds))
  names(output) <- c(paste0("n", conds[,1],"_", "J", conds[,2]))
  out_IW <- out_HW <- out_MF <- output
  # you have to fix the sampling functions to accomodate for different J sies!
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(190331) 
  
  for (outps in 1:nrow(conds)) {
    # Selection of clusters and observations
    n_goal <- conds[outps, 1]
    J_goal <- conds[outps, 2]
    if(n_goal == 46){
      dat <- dat_ALL
    } else{
      clusters_goal <- sample(unique(dat_ALL$cluster), n_goal)
      obs_goal <- unique(dat_ALL$xvec)[1:J_goal]
      dat <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    }
    
    # Define Data for function
    dat_yvec     <- dat$yvec
    dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
    dat_J        <- length(unique(dat$xvec))
    dat_n        <- nrow(dat)/dat_J
    dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
    dat_subjects <- dat$cluster
    
    # Sampling Repetitions
    MCMC_reps   <- 50
    MCMC_burnin <- 1/10
    # IW
    out_IW[[outps]] <- mcmapply(MCMC_invWish, IW_PR,
                       MoreArgs = list(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1),
                       SIMPLIFY = FALSE)
    # HW
    quiet(
      out_HW[[outps]] <- MCMC_HWprior(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin,
                             iniv = 1)
    )
    out_MF[[outps]] <- mcmapply(MCMC_matF, MF_PR,
                       MoreArgs = list(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 0),
                       SIMPLIFY = FALSE)
  }
  
  str(out_IW)
    out_IW$n46_J6$IW_e2$PD_Psi_sd[, 1]
    
  str(out_HW)
    out_HW$n46_J6$PD_Psi_sd[, 1]
    
  str(out_MF)
    out_MF$n46_J6$MF_e1$PD_Psi_sd[, 1]
  
    
  # Check consistency with resulst of priors with old fitting method Previous results of prior
  # Fit with invWishart Prior
  # CHECKed: to obtain same prior as in old fit use the following prior list(e=1,S0=diag(2))
  out_IW <- MCMC_invWish(yvec = dat_yvec, Xmat = dat_Xmat, Zi = dat_Zi, J = dat_J, n = dat_n, samsize = MCMC_reps, burnin = MCMC_burnin, iniv = 1,
               B0 = list(e=1,S0=diag(2)))
  par(mfcol = c(2, 5))
  plot(density(out_IW$PD_Psi_sd[,1]),xlim = c(0, 20), ylim = c(0,.6))
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
  plot(density(outPN$PD_Psi_sd[,1]),xlim = c(0, 20), ylim = c(0,.3))
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
  MCMC_reps   <- 1e4
  MCMC_burnin <- 1/10
  
# Fit model ---------------------------------------------------------------
  
go <- Sys.time()

for (outps in 1:length(conds)) {
  # Reduce number of clusters
  #outps <- 5
  # if(conds[outps] != "ALL"){
  #   if(conds[outps] == "8"){
  #     dat           <- dat_ALL[dat_ALL$cluster %in% c(101, 117, 505, 302, 335, 338, 350, 353), ]
  #   }
  #   if(conds[outps] == "4"){
  #     dat           <- dat_ALL[dat_ALL$cluster %in% c(504, 319, 328, 337), ]
  #   } else {
  #     n_goal        <- as.numeric(conds[outps])                       # how many clusters do you want to work with?
  #     clusters_goal <- sample(unique(dat_ALL$cluster), n_goal)        # sample that many from full dataset.
  #     dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal, ]
  #   }
  # } else {
  #   dat           <- dat_ALL                                        # First condtion uses all
  # }
  clusters_goal <- clusters4cond[[outps]]
  dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal, ]  
  
  # Define Data for function
  dat_yvec     <- dat$yvec
  dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
  dat_J        <- length(unique(dat$xvec))
  dat_n        <- nrow(dat)/dat_J
  dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
  dat_subjects <- dat$cluster
  
  # Define prior (only prior element that needs dat_* elements, so in loop)
  #Rstar <- solve(t(dat_Zi)%*%dat_Zi)
  
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