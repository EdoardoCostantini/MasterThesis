### Project:     Master Thesis
### Object:      Posterior and prior plots for the components of the random effects vairance covairance matrix
### Description: Contains the functions and calls for the plots. You can sue this syntax to plot the results of 
###              you saved from the modelfitting file
### Date:        2019-03-20

# Set up ####
library(MCMCpack)
library(logspline)
library(VarianceGamma)
library(ddpcr)
source("./R/190318-priorplot-functions.R")

# Load Results
  output <- readRDS("./output/riesbydata2019-04-09-ncond_7-rep_10000.rds")
  # NEW! This is with more repetitions, and correct R*
  output <- readRDS("./output/riesbydata2019-04-15-ncond_7-rep_10000.rds")
  output <- readRDS("./output/riesbydata2019-04-27-ncond_7-rep_10000.rds")
  output <- readRDS("./output/riesbydata2019-05-03-ncond_7-rep_10000.rds") # with corss sectional adapted function
  # "best" point estiamtes
  best_est <- output$out_lme4$`n_46:J_6`$lme_Psi_sd
  sd_intercept  <- best_est[1,1]
  sd_slope      <- best_est[2,2]
  cor_int_slope <- best_est[1,2]
    
# SD posteriors and traceplot ---------------------------------------------
  ncond     <- length(output[[1]]) # how many n of clusters conditions are there
  nprior    <- length(output$out_IW$`n_46:J_6`) + length(output$out_MF$`n_46:J_6`) + 1 # how many priors will be considered
  sdseq  <- sqrt(seq(0, 3000, length = 1e5)) # sequence of plausible values of sd to plot HW and mat-F priors
  corseq  <- seq(-1, 1, length = 1e2) # sequence of plausible values of cor to plot HW prior
  
  # Output list levels
  PT  <- 2 # range: 3
  CON <- 2 # range: 7
  PV  <- 4 # range: 2, 1, and 5
  pat <- 4 # range: 1 (just sd is of interest)
  par <- 1 # range: 3 (sd intercept, corr, corr, sd slope)
  
# Inv-Wish ----------------------------------------------------------------
  
  output$out_IW$`n_46:J_6`$IW_eg$hypepar
  output$out_IW$`n_46:J_6`$IW_e1$hypepar
  output$out_IW$`n_46:J_6`$IW_e01$hypepar
  output$out_IW$`n_46:J_6`$IW_e0001$hypepar
  
  #> SDs ####
  PT      <- which(names(output) == "out_IW")
  npar    <- ncol(output$out_IW$`n_46:J_6`$IW_eg$PD_Psi_sd) - 2
  npriors <- length(output$out_IW$`n_46:J_6`)
  
  pdf(paste0("./output/graphs/riesbydata","_SDs_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
      height = 10, width = 17.5)
  
  par(mfcol = c(npar, npriors))
  for (CON in 1:ncond) {
    for (PV in 1:npriors) {
      for (par in 1:2) {
        x <- output[[PT]][[CON]][[PV]][[4]][, c(1,4)]
        plot(density(x[,par]),
             xlim = c(0, 15), ylim = c(0,.5),
             main = paste("w/", names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
             xlab = colnames(x)[par])
        
        # Prior Plot
        prior <- output[[PT]][[CON]][[PV]]$hypepar
        if(PV == 1){
          IW_prior_draws <- sqrt(rinvgamma(1e7,
                                           shape = (2-1+prior$e -1)/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                           scale = prior$S0[par, par]/2))
          h <- hist(IW_prior_draws[order(IW_prior_draws)][1:8500000], # range(IW_prior_draws[order(IW_prior_draws)][1:3600000])
                                                                      # selects an interval that goes from a min close to 0 to a max close to 25
                    breaks = 100, plot = FALSE)
          lines(h$mids, h$density, type = "l", col = "Gray")
        }
        if(PV == 2){
          IW_prior_draws <- sqrt(rinvgamma(1e7,
                                           shape = (2-1+prior$e-1)/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                           scale = prior$S0[par, par]/2))
          h <- hist(IW_prior_draws[order(IW_prior_draws)][1:9700000], #[1:9700000] for 1
                    breaks = 1000, plot = FALSE)
          lines(h$mids, h$density, type = "l", col = "Gray")
        }
        if(PV == 3){
          IW_prior_draws <- sqrt(rinvgamma(1e7,
                                           shape = (2-1+prior$e-1)/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                           scale = prior$S0[par, par]/2))
          h <- hist(IW_prior_draws[order(IW_prior_draws)][1:3600000], #[1:3600000] for 1
                    breaks = 100, plot = FALSE)
          lines(h$mids, h$density, type = "l", col = "Gray")
        }
        if(PV == 4){
          IW_prior_draws <- sqrt(rinvgamma(1e7,
                                           shape = (2-2+prior$e)/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                           scale = prior$S0[par, par]/2))
          h <- hist(IW_prior_draws[order(IW_prior_draws)][1:67500], #[1:8500000] for 1
                    breaks = 100, plot = FALSE)
          lines(h$mids, h$density, type = "l", col = "Gray")
        }
        if(PV == 5){
          IW_prior_draws <- sqrt(rinvgamma(1e7,
                                           shape = prior$nu/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                           scale = prior$S0[par, par]/2))
          h <- hist(IW_prior_draws[order(IW_prior_draws)][1:9950000], #[1:8500000] for 1
                    breaks = 100, plot = FALSE)
          lines(h$mids, h$density, type = "l", col = "Gray")
        }
        # Point Estimates
        MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(1,4)]
        post_mean <- colMeans(x)
        post_median <- apply(x, 2, median)
        points(c(MLE[par], post_mean[par], post_median[par], best_est[par, par]), 
               c(0, 0, 0, 0), 
               pch = c(1, 16, 15, 4), cex = 2)
      }
    }
  }
  
  dev.off()
  
  #> Correlations ####
  PT      <- which(names(output) == "out_IW")
  npar    <- ncol(output$out_IW$`n_46:J_6`$IW_eg$PD_Psi_sd) - 3
  npriors <- length(output$out_IW$`n_46:J_6`)
  
  pdf(paste0("./output/graphs/riesbydata","_COR_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
      height = 5, width = 20)
  
  par(mfcol = c(npar, npriors))
  for (CON in 1:ncond) {
    for (PV in 1:npriors) {
        x <- output[[PT]][[CON]][[PV]][[4]][, c(2)]
        # Posteriro Plot
        # Hist version
        hist(x, # what is the difference when using density instead of logspline. Plots are similar but logspline has soo many more spikes
             freq = F,
             xlim = c(-1, 1), ylim = c(0,2),
             main = paste(names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
             xlab = paste0("RI-", "RS ", colnames(output[[PT]][[CON]][[PV]][[4]])[2]))
        # Prior Plot
        prior <- output[[PT]][[CON]][[PV]]$hypepar
        nu = 2 - 1 + prior$e; S0 = prior$S0 # the parameters for the Inv-Wishart, will use to sample the matrix itself, then select the corr
        reps <- 1e4; prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
        if(PV != 4){
          for (i in 1:reps) {
            Psi <- matrixsampling::rinvwishart(n = 1, nu = nu, Omega = S0)[,,1]
            prior_draws[i, ] <- c(Psi)
          }
          
          prior_draws_corr <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
          h <- hist(prior_draws_corr, breaks = 15, plot = FALSE)
          lines(h$mids, h$density, type = "l", xlim = c(-1, 1), col = "Black", lty = 1)
        }
          
        # Point Estimates
        MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(2)]
        post_mean <- mean(x)
        post_median <- median(x)
        points(c(MLE, post_mean, post_median, best_est[1, 2]), 
               c(-.05, -.05, -.05, -.05), 
               pch = c(1, 16, 15, 4), cex = 2)
    }
  }
  
  dev.off()
    
# Mat-f -------------------------------------------------------------------

  # Priors used
  output$out_MF$`n_46:J_6`$MF_pn$hyperpar
  output$out_MF$`n_46:J_6`$`MF_e=1`$hyperpar
  output$out_MF$`n_46:J_6`$`MF_e=0.5`$hyperpar
  output$out_MF$`n_46:J_6`$`MF_e=0.1`$hyperpar
  output$out_MF$`n_46:J_6`$`MF_e=0.01`$hyperpar
  output$out_MF$`n_46:J_6`$`MF_R*`$hyperpar
  
  str(output$out_MF)
  output$out_MF$`n_8:J_4`$MF_e.5$PD_Psi
  
  PT <- 3 # selects the matf prior list
  npriors <- length(output$out_MF$`n_46:J_6`)
  
  #> SDs ####
  npar <- 2 # intercepts and slopes
  pdf(paste0("./output/graphs/riesbydata","_SDs_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
      height = 10, width = 20)
  for (CON in 1:ncond) {
    par(mfcol = c(npar, npriors))
    for (PV in 1:npriors) {
      for (par in 1:2) {
        if(is.null(output[[PT]][[CON]][[PV]]) == TRUE){
          plot(seq(0,15, length.out = 10), seq(0,.5, length.out = 10),
               main = "ERROR")
        } else {
          x <- output[[PT]][[CON]][[PV]][[4]][, c(1,4)]
          # Posteriro Plot
          plot(density(x[,par]),
               xlim = c(0, 15), ylim = c(0, .5),
               main = paste(names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
               xlab = colnames(x)[par])
          # Prior Plot
          prior <- output[[PT]][[CON]][[PV]]$hyperpar
          lines(sdseq,
                df_freeb_SD(sdseq, nu = prior$nu, d = prior$d, 
                            b = prior$S0[par, par]), # the prior guess goes here! and the guess is not squared!!
                                                     # You verified with plotting the prior by sampling.
                type = "l", col = "Gray")
          
          # Point Estimates
          MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(1,4)]
          post_mean <- colMeans(x)
          post_median <- apply(x, 2, median)
          points(c(MLE[par], post_mean[par], post_median[par], best_est[par, par]), 
               c(0, 0, 0, 0), 
               pch = c(1, 16, 15, 4), cex = 2)
        }
      }
    }
  }
  dev.off()
  
  #> Correlations ####
  PT <- 3
  npriors <- length(output$out_MF$`n_46:J_6`)
  k <- 2; npar <- 1
  
  # Samples for priors
    set.seed(20190426)
    prior_draws_corr <- vector("list", npriors)
    for (PV in 1:npriors) {
      prior <- output[[PT]][[CON]][[PV]]$hyperpar
      B0 <- prior$S0; d <- prior$d; nu0 <- prior$nu; k <- ncol(prior$S0); reps <- 1e4
      prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
      for (ss in 1:reps) {
        Omega <- matrixsampling::rwishart(1, nu = nu0, Sigma = B0)[,,1] + 1e-6*diag(2)
        prior_draws[ss, ] <- matrixsampling::rinvwishart(1, nu = d + k - 1, Omega = Omega)[,,1]
      }
      prior_draws_corr[[PV]] <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
    }
  str(prior_draws_corr)

  pdf(paste0("./output/graphs/riesbydata","_COR_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
        height = 5, width = 20)
  par(mfcol = c(npar, npriors))
  for (CON in 1:ncond) {
    for (PV in 1:npriors) {
      x <- output[[PT]][[CON]][[PV]][[4]][, c(2)]
      # Posteriro Plot
      hist(x, # what is the difference when using density instead of logspline. Plots are similar but logspline has soo many more spikes
           freq = F,
           xlim = c(-1, 1), ylim = c(0,3),
           main = paste(names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
           xlab = paste0("RI-", "RS ", colnames(output[[PT]][[CON]][[PV]][[4]])[2]),
           breaks = 20)
      # Prior Plot
      h <- hist(prior_draws_corr[[PV]], breaks = 20, plot = F)
      lines(h$mids, h$density, type = "l", col = "Black", lty = 1)
      # Point Estimates
      MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(2)]
      post_mean <- mean(x)
      post_median <- median(x)
      points(c(MLE, post_mean, post_median, best_est[1, 2]), 
             c(-.05, -.05, -.05, -.05), 
             pch = c(1, 16, 15, 4), cex = 2)
    }
  }
  dev.off()
  
  output$out_MF$`n_46:J_6`$`MF_R*`$hyperpar$S0
  output$out_MF$`n_8:J_6`$`MF_R*`$hyperpar$S0
  output$out_MF$`n_8:J_3`$`MF_R*`$hyperpar$S0
  output$out_MF$`n_4:J_6`$`MF_R*`$hyperpar$S0
  output$out_MF$`n_4:J_3`$`MF_R*`$hyperpar$S0
  
  # HW ####
  #> SDs ####
  PT <- 4
  npar <- 2
  CON <- 1
  pdf(paste0("./output/graphs/riesbydata","_SDs_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
        height = 7.5, width = 21.5)
  par(mfcol = c(npar, ncond))
  for (CON in 1:ncond) {
    for (par in 1:2) {
      x <- output[[PT]][[CON]][[4]][, c(1,4)]
      # Posteriro Plot
      plot(density(x[,par]),
           xlim = c(0, 15), ylim = c(0,.5),
           main = paste(names(output[[PT]])[[CON]], names(output)[PT]),
           xlab = colnames(x)[par])
      # Prior Plot
      lines(sdseq, dt_folded_sigma(sdseq, nu = 2, mu = 0, sigma = 10**5),
            type = "l", col = "Gray")
      # Point Estimates
      MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(1,4)]
      post_mean <- colMeans(x)
      post_median <- apply(x, 2, median)
      points(c(MLE[par], post_mean[par], post_median[par], best_est[par, par]), 
               c(0, 0, 0, 0), 
               pch = c(1, 16, 15, 4), cex = 2)
    }
  }
  dev.off()
  
  #> Correlations ####
  PT <- 4; npar <- 1
  pdf(paste0("./output/graphs/riesbydata","_COR_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
        height = 6, width = 40)
  par(mfcol = c(npar, ncond))
  for (CON in 1:ncond) {
      x <- output[[PT]][[CON]][[4]][, c(2)]
      # Posterior Plot
      # Hist version
      hist(x, # what is the difference when using density instead of logspline. Plots are similar but logspline has soo many more spikes
           freq = F,
           xlim = c(-1, 1), ylim = c(0,2),
           main = paste(names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
           xlab = paste0("RI-", "RS ", colnames(output[[PT]][[CON]][[PV]][[4]])[2]))
      # Prior Plot
      lines(corseq, dunif(corseq, -1, 1),
            type = "l", col = "Gray")
      # Point Estimates
      MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(2)]
      post_mean <- mean(x)
      post_median <- median(x)
      points(c(MLE, post_mean, post_median, best_est[1, 2]), 
             c(-.05, -.05, -.05, -.05), 
             pch = c(1, 16, 15, 4), cex = 2)
  }
  dev.off()
  
# Traceplot (Flexible) ####
  str(output$out_MF$`n_46:J_6`$MF_pn$PD_Psi_sd[, 4])
  output$out_MF$`n_46:J_6`$MF_pn$hyperpar
  par(mfcol = c(2, 2))
  PT  <- 3 # Prior type (2 = IW, 3 = MF, 4 = HW)
  CON <- 1 # Condition (1 to 7)
  PV  <- 1 # Prior version (depends on prior type)
  par <- 4 # Parameter of interest (intercept sd, corr, corr, slope sd)
  ylim_ob <- matrix(c(0, 25, -1, 1, -1, 1, 0, 3), ncol = 2, byrow = TRUE)
  
  x <- output[[PT]][[CON]][[PV]][[4]][, par]
  #x <- output_loop_prior$MF_e.1$PD_Psi_sd[, par]
  plot(1:1e4, x,"l",
       main = paste(names(output)[PT], "w/", names(output[[PT]])[[CON]], "clusters"),
       xlab = colnames(x),
       xlim = c(0, 1e4), ylim = ylim_ob[par, ])

  
  output <- readRDS("./output/gpadata2019-04-30-ncond_5-rep_10000.rds")
  x <- output$out_MF$`n_200:J_6`$MF_pn$PD_Psi_sd[, 4]
  plot(1:1e4, x,"l",
       xlab = colnames(x),
       xlim = c(0, 1e4), ylim = c(0,5))