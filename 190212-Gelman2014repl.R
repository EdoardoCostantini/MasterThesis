# Title:       Replicate findings in Gelman et al 2014 (adn Gelman 2006, )appendix C)
# Description: Trying to replicate to use Stan for modelling

library("rstan")

# C.2 Fitting a hierarchical model in Stan --------------------------------
set.seed(0)
# Data
  schools_dat <- list(J     = 8,
                      y     = c(28, 8, -3, 7, -1, 1, 18, 12),
                      sigma = c(15, 10, 16, 11, 9, 11, 10, 18)) 
  # schools_dat <- schools_dat3 <- list(J     = 3,
  #                                     y     = c(28, 8, -3),
  #                                     sigma = c(15, 10, 16)) 
# Models
  eightschools.unif <- "
  data {
    int<lower=0> J;          // number of schools
    real y[J];               // vector of estimated treatment effects
    real<lower=0> sigma[J];  // s.e.’s of effect estimates
  }
  parameters {
    real mu;                 // population mean
    real<lower=0> tau;       // population sd
    vector[J] eta;           // school-level errors
  }
  transformed parameters {
    vector[J] theta;         // school effects
    theta = mu + tau*eta;
  }
  model {
    eta ~ normal(0, 1);
    y ~ normal(theta, sigma);
  }"
  
  eightschools.gamm <- "
  data {
    int<lower=0> J;          // number of schools
    real y[J];               // vector of estimated treatment effects
    real<lower=0> sigma[J];  // s.e.’s of effect estimates
  }
  parameters {
    real mu;                 // population mean
    real<lower=0> tau;       // population sd
    vector[J] eta;           // school-level errors
  }
  transformed parameters {
    vector[J] theta;         // school effects
    theta = mu + tau*eta;
  }
  model {
    eta ~ normal(0, 1);
    y ~ normal(theta, sigma);
  }"
  
  eightschools.hC <- "
    data {
      int<lower=0> J;          // number of schools
      real y[J];               // vector of estimated treatment effects
      real<lower=0> sigma[J];  // s.e.’s of effect estimates
    }
    parameters {
      real mu;                 // population mean
      real<lower=0> tau;       // population sd
      vector[J] eta;           // school-level errors
      real<lower=1> nu;
    }
    transformed parameters {
      vector[J] theta;         // school effects
      theta = mu + tau*eta;
    }
    model {
      eta ~ student_t(nu,0,1);
      y ~ normal(theta, sigma);
    }"
  
# Fit Model
  schools_fit <- stan(model_code = eightschools.hC,
                      data = schools_dat, iter = 2000, chains = 4)
  print(schools_fit)
  plot(schools_fit)
# Accessing the posterior simulations in R
  schools_sim <- extract(schools_fit, permuted=TRUE)
  attributes(schools_sim)
  # Plot
  hist(schools_sim$tau)
  # Posterior probability that the effect is larger in school A than in school C
  mean(schools_sim$theta[,1] > schools_sim$theta[,3])
  
# Posterior predictive simulations and graphs in R
  # Replicated data in the existing schools
  n_sims <- length(schools_sim$lp__)
  y_rep <- array (NA, c(n_sims, schools_dat$J))
  for (s in 1:n_sims){
    y_rep[s,] <- rnorm (schools_dat$J, schools_sim$theta[s,], schools_dat$sigma)
  }
  # Plot posterior predictive check
  par (mfrow=c(5,4), mar=c(4,4,2,2))
  hist (schools_dat$y, xlab = "", main = "y")
  for (s in 1:19){
    hist (y_rep[s,], xlab="", main=paste("y_rep",s))
  }
  #upper-left histogram displays the observed data, and the other 19 histograms
  #are posterior predictive replications
  
  
  