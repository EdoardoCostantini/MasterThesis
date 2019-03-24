### Project:     Master Thesis
### Object:      Posterior and prior plots for the components of the random effects vairance covairance matrix
### Description: Contains the functions and calls for the plots. You can sue this syntax to plot the results of 
###              you saved from the modelfitting file
### Date:        2019-03-20

# Set up
library(MCMCpack)
source("./R/190318-priorplot-functions.R")

# Load Results
  
  # Riesbydata (depresion)
    #output <- readRDS("./output/riesbydata_rep10000.rds")
  MCMC_reps <- 10000
  IW_nu0    <- 1
  output_filename <- paste0("./output/riesbydata-cond_46_30_20_8_4-rep_",
                            MCMC_reps,
                            "-IWnu0_",
                            IW_nu0,
                            ".rds")
  #output_filename <- "./output/riesbydata-cond_46_30_20_8_4-rep_10000-HWnu0_1.rds"
    # this for the HW prior with nu0 = 1
# Load output
    output <- readRDS(output_filename)[-6] #exclude last object because it contains the priors
      str(output)
    pg <- readRDS(output_filename)[[6]] #last object contains the priors

# SD posteriors and traceplot ---------------------------------------------

# parameter of interest selection (eg. do you want to see the variance or the sd? correlation?)
  #which_pri <- 5 # which prior
  #which_con <- 5 # which condition
  which_par <- 4  # which object (which parameters type)

  ncond  <- length(output[[1]]) # how many n of clusters conditions are there
  npar   <- 2 # how many parameters are of interest (random slope and random intercept)
  nprior <- length(output) # how many priors will be considered
  
  sdseq  <- sqrt(seq(0, 3000, length = 100000)) # important to plot the piror
  IW_prior_draws <- rinvgamma(100000,
                              shape = IW_nu0/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                              scale = pg[2, 1]/2)
  plot(density(IW_prior_draws),
       xlim = c(0, 100),#, ylim = c(0,2),
       type = "l", col = "Gray")
  
  par(mfcol = c(npar, nprior))
  
# > Posterior Distirbutions ####  
  for (cond in 1:ncond) {
    #cond <- 1
    for (prior in 1:nprior) {
      #prior <- 1      
      for (param in 1:npar) {
        #param <- 1        
        x <- output[[prior]][[cond]][[which_par]][,-c(2,3)]
    
        plot(density(x[,param]),
             xlim = c(0, 10), ylim = c(0,2),
             main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
             xlab = colnames(x)[param])
        # abline(v = mean(x[,param]), col = "blue", lwd = 1)
        # abline(v = median(x[,param]), col = "red", lwd = 1)
        if(prior == 1){
          lines(sdseq,
                dt_folded(sdseq,
                          nu = 2,
                          mu = 0,
                          sigma2 = pg[prior, param]), # the prior guess goes here! Arbitrarly large values induce arbitrarly weak priors
                                                      # see Huang and Wand 2013 p.3 and property 2 p.4
                type = "l", col = "Gray")
        }
        if(prior == 2){
          # IW prior (on sd)
          lines(density(IW_prior_draws),
                type = "l", col = "Gray")
        }
        if(prior >= 3){
          lines(sdseq,
                df_freeb_SD(sdseq, 
                            nu = 2, 
                            d = 1, 
                            b = sqrt(pg[prior, param])), # the prior guess goes here!
                type = "l", col = "Gray")
        }
      }
    }
  }
# > Traceplots ####
  for (cond in 1:ncond) {
    #cond <- 1
    for (prior in 1:nprior) {
      #prior <- 1      
      for (param in 1:npar) {
        #param <- 1        
        x <- output[[prior]][[cond]][[which_par]][,-c(2,3)]
        
        plot(1:MCMC_reps, x[, param],"l",
             main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
             xlab = colnames(x)[param],
             xlim = c(0, MCMC_reps), ylim = c(0,100))
      }
    }
  }

# > Estimates ####
  estimates_sd <- vector("list", ncond)
  names(estimates_sd) <- c("cond46","cond30","cond20","cond8","cond4")
  for (cond in 1:ncond) {
    #cond <- 1
    store_estiamtes <- matrix(rep(NA, 20), ncol = 5, nrow = 4)
    for (prior in 1:nprior) {
      #prior <- 3
      rownames(store_estiamtes) <- c("mean_RI", "mean_RS", "median_RI", "median_RS")
      colnames(store_estiamtes) <- c("HW", "IW", "MFPN", "MFEG", "MFR*")
      x <- output[[prior]][[cond]][[which_par]][,-c(2,3)]
      store_estiamtes[1:2, prior] <- colMeans(x)  
      store_estiamtes[3:4, prior] <- apply(x, 2, median)
    }
    estimates_sd[[cond]] <- round(store_estiamtes, 3)
  }
  
# Corr posteriors and traceplot -------------------------------------------

# parameter of interest selection (eg. do you want to see the variance or the sd? correlation?)
  #which_pri <- 5 # which prior
  #which_con <- 5 # which condition
  which_par <- 4  # which object (which parameters type)
  ncond  <- length(output[[1]]) # how many n of clusters conditions are there
  npar   <- 1 # how many parameters are of interest (random slope and random intercept)
  nprior <- 5#length(output) # how many priors will be considered
  par(mfcol = c(npar, nprior))
  
  corseq  <- seq(-1, 1, length = 100) # important to plot the piror
  
# > Posterior Distirbutions ####  
  for (cond in 1:ncond) {
#cond <- 1
    for (prior in 1:nprior) {
#prior <- 1      
        x <- output[[prior]][[cond]][[which_par]][,2]
        
        plot(density(x),
             xlim = c(-1, 1), ylim = c(0,2),
             main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
             xlab = colnames(x))
        # abline(v = mean(x), col = "blue", lwd = 1)
        # abline(v = median(x), col = "red", lwd = 1)
        if(prior == 1){
          lines(corseq,
               dunif(corseq, -1, 1),
               type = "l", col = "Gray")
        }
        # if(prior >= 3){
        #   lines(sdseq,
        #         df_freeb_SD(sdseq, 
        #                     nu = 2, d = 1, 
        #                     b = sqrt(pg[prior, param])), # the prior guess goes here!
        #         type = "l", col = "Gray")
        # }
    }
  }
  
# > Traceplots ####
  for (cond in 1:ncond) {
    #cond <- 1
    for (prior in 1:nprior) {
      
      #prior <- 1      
      x <- output[[prior]][[cond]][[which_par]][,2]
      
      plot(1:MCMC_reps, x,"l",
           main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
           xlab = colnames(x),
           xlim = c(0, MCMC_reps), ylim = c(-1,1))
    }
  }
  
# > Estimates ####
  estimates_corr <- vector("list", ncond)
  names(estimates_corr) <- c("cond46","cond30","cond20","cond8","cond4")
  for (cond in 1:ncond) {
    #cond <- 1
    store_estiamtes <- matrix(rep(NA, 10), ncol = 5, nrow = 2)
    rownames(store_estiamtes) <- c("mean_CORR", "median_CORR")
    colnames(store_estiamtes) <- c("HW", "IW", "MFPN", "MFEG", "MFR*")
    for (prior in 1:nprior) {
      #prior <- 1      
      x <- output[[prior]][[cond]][[which_par]][,2]
      store_estiamtes[1, prior] <- mean(x)
      store_estiamtes[2, prior] <- median(x)
    }
    estimates_corr[[cond]] <- round(store_estiamtes, 3)
  }
  
# Experiments
  prior <- 2 # which prior
  cond  <- 1 # which condition
  which_par <- 4 # which object (which parameters type)
  x <- output[[prior]][[cond]][[which_par]][,-c(2,3)]
        param <- 1
        plot(density(x[,param]),
             xlim = c(0, 10), ylim = c(0,2),
             main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
             xlab = colnames(x)[param])
        # HW prior (on sd)
        lines(sdseq,
              dt_folded(sdseq,
                        nu = 1,
                        mu = 0,
                        sigma2 = pg[which_pri, param]), # the prior guess goes here! Arbitrarly large values induce arbitrarly weak priors
              # see Huang and Wand 2013 p.3 and property 2 p.4
              type = "l", col = "Gray")
        # IW prior (on sd)
        IW_prior_draws <- rinvgamma(5000,
                                    shape = 2,
                                    scale = 1)
        lines(density(IW_prior_draws),
              type = "l", col = "Gray")
        
        # Univaraite F
        lines(sdseq,
                df_freeb_SD(sdseq, 
                            nu = 2, d = 1, 
                            b = sqrt(pg[which_pri, param])), # the prior guess goes here!
                type = "l", col = "Gray")
