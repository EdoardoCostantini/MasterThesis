### Project:     Master Thesis
### Object:      Posterior and prior plots for the components of the random effects vairance covairance matrix
### Description: Contains the functions and calls for the plots. You can sue this syntax to plot the results of 
###              you saved from the modelfitting file
### Date:        2019-03-20

# Set up
library(MCMCpack)
library(logspline)
library(VarianceGamma)
library(ddpcr)
source("./R/190318-priorplot-functions.R")

# Load Results
  output <- readRDS("./output/riesbydata2019-04-03-ncond_7-rep_5000.rds")
  output <- readRDS("./output/riesbydata2019-04-04-ncond_7-rep_5000_Frank2.rds")
  output <- readRDS("./output/riesbydata2019-04-04-ncond_7-rep_10000_Frank2.rds")
  output <- readRDS("./output/riesbydata2019-04-05-ncond_7-rep_5000.rds")
  output <- readRDS("./output/riesbydata2019-04-05-ncond_7-rep_5000correctIW.rds")
  
  output <- readRDS("./output/riesbydata2019-04-07-ncond_7-rep_5000_DEFINITIVE.rds")
  

    output_OLD <- readRDS(output_filename)[-6] #exclude last object because it contains the priors (load separately)
      str(output)
    pg <- readRDS(output_filename)[[6]]    #last object contains the priors

    # COMPARE
    output_OLD$prior_HW$`46`$PD_Psi_sd
    output$out_HW$`n_46:J_6`$PD_Psi_sd
    output_OLD$prior_IW$`46`$PD_Psi_sd
    output$out_IW$`n_46:J_6`$IW_eg
    output[[2]][[1]][[2]][[4]]
    output[[PT]][[CON]][[PV]][[PAR]]
    # PRIOR TYPE # CONDITION # PRIOR VERSION # PARAMETER CHOICE
    
# SD posteriors and traceplot ---------------------------------------------
  ncond     <- length(output[[1]]) # how many n of clusters conditions are there
  nprior    <- length(output$out_IW$`n_46:J_6`) + length(output$out_MF$`n_46:J_6`) + 1 # how many priors will be considered
  sdseq  <- sqrt(seq(0, 3000, length = 1e5)) # sequence of plausible values of sd to plot HW and mat-F priors
  corseq  <- seq(-1, 1, length = 1e2) # sequence of plausible values of cor to plot HW prior
  
  # Output list levels
  PT  <- 2 # range: 3
  CON <- 6 # range: 7
  PV  <- 4 # range: 2, 1, and 5
  pat <- 4 # range: 1 (just sd is of interest)
  par <- 1 # range: 3 (sd intercept, corr, corr, sd slope)
  
  #Inv-Wish ####
  PT <- 2
  output$out_IW$`n_46:J_6`$IW_eg$hypepar
  output$out_IW$`n_46:J_6`$IW_e1$hypepar
  output$out_IW$`n_46:J_6`$IW_e01$hypepar
  output$out_IW$`n_46:J_6`$IW_e0001$hypepar
  
  #> SDs ####
  npar <- 2
  npriors <- 4
  
  pdf(paste0("./output/graphs/riesbydata","_SDs_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
      height = 10, width = 17.5)
  
  par(mfcol = c(npar, npriors))
  for (CON in 1:ncond) {
    for (PV in 1:npriors) {
      for (par in 1:2) {
        
        x <- output[[PT]][[CON]][[PV]][[4]][, c(1,4)]
        # Posteriro Plot
        plot(density(x[,par]),
             xlim = c(0, 20), ylim = c(0,.5),
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
        # Point Estimates
        MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(1,4)]
        post_mean <- colMeans(x)
        post_median <- apply(x, 2, median)
        points(c(MLE[par], post_mean[par], post_median[par]), c(0, 0, 0), pch = c(1, 16, 15), cex = 2)
      }
    }
  }
  
  dev.off()
  
  #> Correlations ####
  PT   <- 2
  npar <- 1
  nprior <- 4
  
  CON <- 1
  PV <- 3
  
  pdf(paste0("./output/graphs/riesbydata","_COR_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
      height = 10, width = 20)
  par(mfcol = c(npar, npriors))
  
  for (CON in 1:ncond) {
    for (PV in 1:npriors) {
        x <- output[[PT]][[CON]][[PV]][[4]][, c(2)]
        # Posteriro Plot
        plot(density(x),
             xlim = c(-1, 1), ylim = c(0,2),
             main = paste(names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
             xlab = paste0("RI-", "RS ", colnames(output[[PT]][[CON]][[PV]][[4]])[2]))
        # Prior Plot
        prior <- output[[PT]][[CON]][[PV]]$hypepar
        if((2 - 1 + prior$e) >= 2){
          nu = 2 - 1 + prior$e; S0 = prior$S0 # the parameters for the Inv-Wishart, will use to sample the matrix itself, then select the corr
          reps <- 1e4; prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
          for (i in 1:reps) {
            Psi   <- riwish(nu, S0)
            prior_draws[i, ] <- c(Psi)
          }
          prior_draws_corr <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
          h <- hist(prior_draws_corr, breaks = 10, plot = FALSE)
          lines(h$mids, h$density, type = "l", xlim = c(-1, 1), col = "Gray", lty = 1)
        }
        # Point Estimates
        MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(2)]
        post_mean <- mean(x)
        post_median <- median(x)
        points(c(MLE, post_mean, post_median), c(0, 0, 0), pch = c(1, 16, 15), cex = 2)
    }
  }
  dev.off()

  # Mat-f ####
  PT <- 4 # selects the matf prior list
  # Priors used
  output$out_MF$`n_46:J_6`$MF_pn$hyperpar
  output$out_MF$`n_46:J_6`$`MF_eg(e1)`$hyperpar
  output$out_MF$`n_46:J_6`$MF_e.5$hyperpar
  output$out_MF$`n_46:J_6`$MF_e.1$hyperpar
  output$out_MF$`n_46:J_6`$`MF_R*`$hyperpar
  str(output$out_MF)
  #> SDs ####
  npar <- 2 # intercepts and slopes
  pdf(paste0("./output/graphs/riesbydata","_SDs_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
      height = 10, width = 20)
  for (CON in 1:ncond) {
    par(mfcol = c(npar, 5))
    for (PV in 1:5) {
      for (par in 1:2) {
        if(is.null(output[[PT]][[CON]][[PV]]) == TRUE){
          plot(seq(0,20, length.out = 10), seq(0,.5, length.out = 10),
               main = "ERROR")
        } else {
          x <- output[[PT]][[CON]][[PV]][[4]][, c(1,4)]
          # Posteriro Plot
          plot(density(x[,par]),
               xlim = c(0, 20), ylim = c(0, .5),
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
          points(c(MLE[par], post_mean[par], post_median[par]), c(0, 0, 0), pch = c(1, 16, 15), cex = 2)
        }
      }
    }
  }
  dev.off()
  
  #> Correlations ####
  CON <- 1; PV <- 5;
  k <- 2
  npar <- 1
  pdf(paste0("./output/graphs/riesbydata","_COR_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
        height = 5, width = 20)
  par(mfcol = c(npar, 5))
  for (CON in 1:ncond) {
    for (PV in 1:5) {
      if(is.null(output[[PT]][[CON]][[PV]]) == TRUE){
        plot(seq(-1,1, length.out = 10), seq(0,2, length.out = 10),
             main = "ERROR")
      } else {
        x <- output[[PT]][[CON]][[PV]][[4]][, c(2)]
        # Posteriro Plot
        plot(density(x),
             xlim = c(-1, 1), ylim = c(0,5),
             main = paste(names(output[[PT]])[[CON]], names(output[[PT]][[CON]])[[PV]]),
             xlab = paste0("RI-", "RS ", colnames(output[[PT]][[CON]][[PV]][[4]])[2]))
        # Prior Plot
        prior <- output[[PT]][[CON]][[PV]]$hyperpar
        reps <- 1e4; Omega <- B <- prior$S0 # matrix(c(21, 0, 0, 9), ncol = 2) # guess based on data exploration and knowledge
        draws <- matrix(rep(NA, reps*4), ncol = 4)
        if(((prior$d + k -1) >= 2) == TRUE){
          for (i in 1:reps) {
            # Psi|Omega,.
            ScaleMatrix = Omega
            PsiInvDraw = rwish(v = prior$d + k - 1,
                               S = solve(ScaleMatrix))
            # Omega|Psi,.
            ScaleOmega = (PsiInvDraw + solve(B))
            Omega = rwish(v = prior$nu + prior$d + k - 1,
                          S = solve(ScaleOmega))
            draws[i, ] <- c(solve(PsiInvDraw))
          }
          prior_draws_corr <- draws[, 2] / (sqrt(draws[, 1]) * sqrt(draws[, 4]))
          h <- hist(prior_draws_corr, breaks = 30, plot = T)
          lines(h$mids, h$density, type = "l", col = "Gray", lty = 1)
        }
        # Point Estimates
        MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(2)]
        post_mean <- mean(x)
        post_median <- median(x)
        points(c(MLE, post_mean, post_median), c(0, 0, 0), pch = c(1, 16, 15), cex = 2)
      }
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
  PT <- 3
  npar <- 2
  pdf(paste0("./output/graphs/riesbydata","_SDs_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
        height = 5, width = 20)
  par(mfcol = c(npar, ncond))
  for (CON in 1:ncond) {
    for (par in 1:2) {
      x <- output[[PT]][[CON]][[4]][, c(1,4)]
      # Posteriro Plot
      plot(density(x[,par]),
           xlim = c(0, 20), ylim = c(0,.5),
           main = paste(names(output[[PT]])[[CON]], names(output)[3]),
           xlab = colnames(x)[par])
      # Prior Plot
      lines(sdseq, dt_folded_sigma(sdseq, nu = 2, mu = 0, sigma = 10**5),
            type = "l", col = "Gray")
      # Point Estimates
      MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(1,4)]
      post_mean <- colMeans(x)
      post_median <- apply(x, 2, median)
      points(c(MLE[par], post_mean[par], post_median[par]), c(0, 0, 0), pch = c(1, 16, 15), cex = 1.5)
    }
  }
  dev.off()
  
  #> Correlations ####
  npar <- 1
  pdf(paste0("./output/graphs/riesbydata","_COR_",names(output)[PT], ".pdf"), # COMMENT OUT THIS PART IF YOU ONLY WANT TO SEE THE PLOTS
        height = 5, width = 20)
  par(mfcol = c(npar, ncond))
  for (CON in 1:ncond) {
      x <- output[[PT]][[CON]][[4]][, c(2)]
      # Posteriro Plot
      plot(density(x),
           xlim = c(-1, 1), ylim = c(0,2),
           main = paste(names(output[[PT]])[[CON]], names(output)[3]),
           xlab = colnames(output[[PT]][[CON]][[4]][, ])[2])
      # Prior Plot
      lines(corseq, dunif(corseq, -1, 1),
            type = "l", col = "Gray")
      # Point Estimates
      MLE <- c(output$out_lme4[[CON]]$lme_Psi_sd)[c(2)]
      post_mean <- mean(x)
      post_median <- median(x)
      points(c(MLE, post_mean, post_median), c(0, 0, 0), pch = c(1, 16, 15), cex = 2)
  }
  dev.off()
  
# Traceplot (Flexible) ####
  par(mfcol = c(2, 2))
  PT  <- 4 # Prior type (2 = IW, 3 = HW, 4, MF)
  CON <- 1 # Condition (1 to 7)
  PV  <- 2 # Prior version (depends on prior type)
  par <- 1 # Parameter of interest (intercept sd, corr, corr, slope sd)
  ylim_ob <- matrix(c(0, 25, -1, 1, -1, 1, 0, 50), ncol = 2, byrow = TRUE)
  
  x <- output[[PT]][[CON]][[PV]][[4]][, par]
  #x <- output_loop_prior$MF_e.1$PD_Psi_sd[, par]
  plot(1:5000, x,"l",
       main = paste(names(output)[PT], "w/", names(output[[PT]])[[CON]], "clusters"),
       xlab = colnames(x),
       xlim = c(0, 5000), ylim = ylim_ob[par, ])