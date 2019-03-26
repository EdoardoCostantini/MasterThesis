### Project:     Master Thesis
### Object:      Posterior and prior plots for the components of the random effects vairance covairance matrix
### Description: Contains the functions and calls for the plots. You can sue this syntax to plot the results of 
###              you saved from the modelfitting file
### Date:        2019-03-20

# Set up
library(MCMCpack)
library(ddpcr)
source("./R/190318-priorplot-functions.R")

# Load Results
  
  # Riesbydata (depresion) output
  IW_nu0 <- 2; MCMC_reps <- 10000; 
  output_filename <- paste0("./output/riesbydata-cond_46_30_20_8_4-rep_",MCMC_reps,
                            "-IWnu0_",IW_nu0,
                            ".rds")
  #output_filename <- "./output/riesbydata-cond_46_30_20_8_4-rep_10000-HWnu0_1.rds"
    # this for the HW prior with nu0 = 1
# Load output
    output <- readRDS(output_filename)[-6] #exclude last object because it contains the priors (load separately)
      str(output)
    pg <- readRDS(output_filename)[[6]]    #last object contains the priors

# SD posteriors and traceplot ---------------------------------------------
  npar      <- 2 # how many parameters are of interest (random slope and random intercept)
  ncond     <- length(output[[1]]) # how many n of clusters conditions are there
  nprior    <- length(output) # how many priors will be considered
  obj_param <- 4  # which parameters "type"
                  # usually interested in the standardized Psi (sd and correlations, instead of var and covar)
  
  sdseq  <- sqrt(seq(0, 3000, length = 100000)) # sequence of plausible values of sd to plot HW and mat-F priors

  IW_prior_draws <- sqrt(rinvgamma(1e7,
                                   shape = IW_nu0/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                   scale = pg[2, 1]/2))
  h <- hist(IW_prior_draws, col = "Gray",
          prob = TRUE,
          xlim = c(0, 20), breaks = 1e4)
        
  par(mfcol = c(npar, nprior))
  
# > Posterior Distirbutions ####  
  for (cond in 1:ncond) {
    for (prior in 1:nprior) {
      for (param in 1:npar) {
        x <- output[[prior]][[cond]][[obj_param]][,-c(2,3)]
        # Posteriro Plot
          plot(density(x[,param]),
               xlim = c(0, 20), ylim = c(0,1),
               main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
               xlab = colnames(x)[param])
        # Point Estimates
          # abline(v = mean(x[,param]), col = "blue", lwd = 1)
          # abline(v = median(x[,param]), col = "red", lwd = 1)
        # Prior Plots
          # HW prior
          if(prior == 1){
            lines(sdseq, dt_folded_sigma(sdseq, nu = 2, mu = 0, sigma = 100),
                  type = "l", col = "Gray")
            # This prior does look like a t disitrbution when zooming in
            # Since it's very vague (arbitrarly high A_k = 100), we have 
            # very little density even at the peak
          }
          # IW prior (on sd)
            if(prior == 2){
              lines(h$mids, h$density, type = "l", col = "Gray", xlim = c(0, 20))
            }
          # Mat-F
          if(prior >= 3){
            lines(sdseq,
                  df_freeb_SD(sdseq, nu = 2, d = 1, 
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

  npar      <- 1                   # number of parameters of interest: 1 (RI-RS correlation)
  ncond     <- length(output[[1]]) # how many conditions are there?
  nprior    <- length(output)      # how many priors will be considered
  obj_param <- 4                   # which object (which parameters type)
  
  par(mfcol = c(npar, nprior))     # graph setting
  
  corseq <- seq(-1,1,length=100)   # sequence to plot HW prior
  
# > Posterior Distirbutions ####  
  for (cond in 1:ncond) { # cond = number of clusters and repeated measures w/in
    for (prior in 1:nprior) {
        x <- output[[prior]][[cond]][[obj_param]][,2] # selects: prior > condition > standardized version > column2
                                                      # e.g.   : HW    > all clust > sd and correlations  > correlations
        # Posterior plots
          plot(density(x),
               xlim = c(-1, 1), ylim = c(0,2),
               main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
               xlab = colnames(x))
        # Point Estimates
          # abline(v = mean(x), col = "blue", lwd = 1)
          # abline(v = median(x), col = "red", lwd = 1)
        # Prior Plots
          if(prior == 1){
            lines(corseq, dunif(corseq, -1, 1), type = "l", col = "Gray")
          }
          if(prior == 2){
            nu = 2; S0 = matrix(c(1, 0, 0, 1), ncol = 2) # guess based on data exploration and knowledge
            reps <- 10000; prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
            for (i in 1:reps) {
              Psi   <- riwish(nu, S0)
              prior_draws[i, ] <- c(Psi)
            }
            prior_draws_corr <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
            h <- hist(prior_draws_corr, breaks = 10, plot = FALSE)
            lines(h$mids, h$density, type = "l", xlim = c(-1, 1), col = "Gray", lty = 1)
          }
          if(prior >= 3){
            nu = 2; delta = 1
            B  = matrix(c(pg[prior, 1], 0, 0, pg[prior, 2]), ncol = 2) # guess based on data exploration and knowledge
            reps <- 10000
            prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
            for (i in 1:reps) {
              Omega <- rwish(nu, B)
              Psi   <- riwish(delta+2-1, Omega)
              prior_draws[i, ] <- c(Psi)
            }
            prior_draws_corr <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
            h <- hist(prior_draws_corr, breaks = 10, plot = FALSE)
            lines(h$mids, h$density, type = "l", xlim = c(-1, 1), col = "Gray", lty = 1)
          }
    }
  }
  
# > Traceplots ####
  for (cond in 1:ncond) {
    for (prior in 1:nprior) {
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
  
# Experiments ####
  prior <- 4 # which prior
  cond  <- 1 # which condition
  which_par <- 4 # which object (which parameters type)
  x <- output[[prior]][[cond]][[which_par]][,-c(2,3)]
  param <- 1
  plot(density(x[,param]),
       xlim = c(0, 20), ylim = c(0,2),
       main = paste(names(output)[prior], "w/", names(output[[prior]])[[cond]], "clusters"),
       xlab = colnames(x)[param])
# SD of random effects
  # HW prior (on sd)
    plot(sdseq,
         dt_folded(sdseq,
                   nu = 1,
                   mu = 0,
                   sigma2 = 100), # the prior guess goes here! Arbitrarly large values induce arbitrarly weak priors
         # see Huang and Wand 2013 p.3 and property 2 p.4
         type = "l", col = "Gray")
  # IW prior (on sd)
    nu0 <- 1.001
    S0 <- 1.001
    IW_prior_draws <- sqrt(rinvgamma(1e6,
                                     shape = nu0/2, # you chose nu0 = 2 in the IW, IG on sd is nu0/2
                                     scale = S0/2))
    plot(density(IW_prior_draws),
         type = "l", col = "Gray",
         xlim = c(0, 20))
    h <- hist(IW_prior_draws, col = "Gray",
              prob = TRUE,
              xlim = c(0, 20), breaks = 1e6)
    plot(h$mids, h$density, type = "l", col = "Gray", xlim = c(0, 20)) # <- THE METHOD!
    
    # When in doubt on the efficacy of the method you have chosen 
    # to plot the prior in the IW case, check that it does what you
    # want with the uniform case
    rnorm_draws <- rnorm(1e6, mean = 0, sd = 1)
    plot(density(rnorm_draws),
         type = "l", col = "black", xlim = c(-5, 5), ylim = c(0, 0.5))
    h <- hist(rnorm_draws, col = "Gray",
              prob = TRUE,
              xlim = c(-5, 5), ylim = c(0, 0.5),
              breaks = 100)
    lines(h$mids, h$density, type = "l", xlim = c(-5, 5), ylim = c(0, 0.5),
          col = "Gray", lty = 2)
    
  # Univaraite F
    plot(sdseq,
          df_freeb_SD(sdseq, 
                      nu = 2, d = 1, 
                      b = sqrt(pg[prior, param])), # the prior guess goes here!
          type = "l", col = "Gray")
# Correlation
  # matrix F
    nu    = 2
    delta = 1
    B     = matrix(c(20, 0, 0, 9), ncol = 2) # guess based on data exploration and knowledge
    reps <- 10000
    prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
    for (i in 1:reps) {
      Omega <- rwish(nu, B)
      Psi   <- riwish(delta+2-1, Omega)
      prior_draws[i, ] <- c(Psi)
    }
    # sd
    prior_draws_sqrt <- sqrt(prior_draws[, c(1,4)])
    viVar <- prior_draws_sqrt[prior_draws_sqrt[,1] < 1000, 1]
    h <- hist(viVar, breaks = 500, xlim = c(0, 50))
    plot(h$mids, h$density, type = "l", xlim = c(0, 50), col = "Gray", lty = 1)
    # correlation
    prior_draws_corr <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
    h <- hist(prior_draws_corr, breaks = 100, xlim = c(-1, 1))
    lines(h$mids, h$density, type = "l", xlim = c(-1, 1), col = "Gray", lty = 1)
    
  # inverse Wishart
    nu    = 2
    S0     = matrix(c(1, 0, 0, 1), ncol = 2) # guess based on data exploration and knowledge
    reps <- 10000
    prior_draws <- matrix(rep(NA, reps*4), ncol = 4)
    for (i in 1:reps) {
      Psi   <- riwish(nu, S0)
      prior_draws[i, ] <- c(Psi)
    }
    # sd
    prior_draws_sqrt <- sqrt(prior_draws[, c(1,4)])
    viVar <- prior_draws_sqrt[prior_draws_sqrt[,1] < 1000, 1]
    h <- hist(viVar, breaks = 500, xlim = c(0, 50))
    plot(h$mids, h$density, type = "l", xlim = c(0, 50), col = "Gray", lty = 1)
    # correlation
    prior_draws_corr <- prior_draws[, 2] / (sqrt(prior_draws[, 1]) * sqrt(prior_draws[, 4]))
    h <- hist(prior_draws_corr, breaks = 10, xlim = c(-1, 1))
    plot(h$mids, h$density, type = "l", xlim = c(-1, 1), col = "Gray", lty = 1)