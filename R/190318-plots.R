### Project:     Master Thesis
### Object:      Posterior and prior plots for the components of the random effects vairance covairance matrix
### Description: Contains the functions and calls for the plots
### Date:        2019-03-20

# Set up
#library(MCMCpack)
source("./R/190318-priorplot-functions.R")

# Load Results
  
  # Riesbydata (depresion)
    #output <- readRDS("./output/riesbydata_rep10000.rds")
    output <- readRDS("./output/riesbydata-cond_46_30_20_8_4-rep_10000v2.rds")[-6]
      str(output)
    pg <- readRDS("./output/riesbydata-cond_46_30_20_8_4-rep_10000v2.rds")[[6]] #last object contains the priors
    MCMC_reps <- 10000 # info in the file name

# SD posteriors and traceplot ---------------------------------------------

# parameter of interest selection (eg. do you want to see the variance or the sd? correlation?)
  #which_pri <- 5 # which prior
  #which_con <- 5 # which condition
  which_par <- 4  # which object (which parameters type)

  ncond  <- length(output[[1]]) # how many n of clusters conditions are there
  npar   <- 2 # how many parameters are of interest (random slope and random intercept)
  nprior <- length(output) # how many priors will be considered
  par(mfcol = c(npar, nprior))
  
  sdseq  <- sqrt(seq(0, 3000, length = 100000)) # important to plot the piror
  
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
        if(prior == 1){
          lines(sdseq,
                dt_folded(sdseq,
                          nu = 2,
                          mu = 0,
                          sigma2 = pg[prior, param]), # the prior guess goes here! Arbitrarly large values induce arbitrarly weak priors
                                                      # see Huang and Wand 2013 p.3 and property 2 p.4
                type = "l", col = "Gray")
        }
        if(prior >= 3){
          lines(sdseq,
                df_freeb_SD(sdseq, 
                            nu = 2, d = 1, 
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

# Corr posteriors and traceplot -------------------------------------------

# parameter of interest selection (eg. do you want to see the variance or the sd? correlation?)
  #which_pri <- 5 # which prior
  #which_con <- 5 # which condition
  which_par <- 4  # which object (which parameters type)
  length(output[[1]])
  ncond  <- length(output[[1]]) # how many n of clusters conditions are there
  npar   <- 1 # how many parameters are of interest (random slope and random intercept)
  nprior <- 5#length(output) # how many priors will be considered
  par(mfcol = c(npar, nprior))
  
  sdseq  <- sqrt(seq(0, 3000, length = 100000)) # important to plot the piror
  
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
        # if(prior == 1){
        #   lines(sdseq,
        #         dt_folded(sdseq,
        #                   nu = 2,
        #                   mu = 0,
        #                   sigma2 = pg[prior, param]), # the prior guess goes here! Arbitrarly large values induce arbitrarly weak priors
        #         # see Huang and Wand 2013 p.3 and property 2 p.4
        #         type = "l", col = "Gray")
        # }
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
  
  
# Experiments
  which_pri <- 1 # which prior
  which_con <- 4 # which condition
  which_par <- 4 # which object (which parameters type)
  x <- output[[which_pri]][[which_con]][[which_par]][,-c(2,3)]
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
        lines(sdseq,
              dinvgamma(sdseq,
                        shape = 2,
                        scale = 1), # the prior guess goes here! Arbitrarly large values induce arbitrarly weak priors
                                         # see Huang and Wand 2013 p.3 and property 2 p.4
              type = "l", col = "Gray")
