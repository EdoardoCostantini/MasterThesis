# Title:       Replicate findings in Gelman et al 2014 appendix C (and Gelman 2006)
# Description: Trying to replicate to use Stan for modelling

library("rstan")

# Models ------------------------------------------------------------------

  model.unif <- "
  data {
    int<lower=0> J;          // number of schools (restricted to be positve)
    real y[J];               // vector of estimated treatment effects
    real<lower=0> sigma[J];  // s.e.’s of effect estimates (restricted to be positve)
  }
  parameters {
    real mu;                 // population mean     (mu)
    real<lower=0> tau;       // population sd       (tau)
    vector[J] eta;           // school-level errors (eta)
  }
  transformed parameters {
    vector[J] theta;         // school effects      (theta.j)
    theta = mu + tau*eta;
  }
  model {
    eta ~ normal(0, 1);      // Stan parametrize N(mean, sd)
    y ~ normal(theta, sigma);
  }"
  
  model.halfC <- "
  data {
    int<lower=0> J;          // number of schools (restricted to be positve)
    real y[J];               // vector of estimated treatment effects
    real<lower=0> sigma[J];  // s.e.’s of effect estimates (restricted to be positve)
  }
  parameters {
    real mu;                 // population mean     (mu)
    real<lower=0> tau;       // population sd       (tau)
    vector[J] eta;           // school-level errors (eta)
  }
  transformed parameters {
    vector[J] theta;         // school effects      (theta.j)
    theta = mu + tau*eta;
  }
  model {
    eta ~ normal(0, 1);      // Stan parametrize N(mean, sd)
    y ~ normal(theta, sigma);
    tau ~ cauchy(0,25);
  }"
  
  model.gamm11 <- "
  data {
    int<lower=0> J;          // number of schools
    real y[J];               // vector of estimated treatment effects
    real<lower=0> sigma[J];  // s.e.’s of effect estimates
  }
  parameters {
    real mu;                 // population mean
    real<lower=0> tau;       // population sd
    vector[J] theta;         // school effects
  }
  model {
    y ~ normal(theta, sigma);
    theta ~ normal(mu, tau);
    mu ~ normal(0, 1000);
    tau ~ inv_gamma(1, 1);
  }"
  
  model.gamm00 <- "
  data {
    int<lower=0> J;          // number of schools
    real y[J];               // vector of estimated treatment effects
    real<lower=0> sigma[J];  // s.e.’s of effect estimates
  }
  parameters {
    real mu;                 // population mean
    real<lower=0> tau;       // population sd
    vector[J] theta;         // school effects
  }
  model {
    y ~ normal(theta, sigma);
    theta ~ normal(mu, tau);
    mu ~ normal(0, 1000);
    tau ~ inv_gamma(.001, .001);
  }"
  

# Fitting the models ------------------------------------------------------
  
  #> 8 Shcools Example ####
  schools_dat <- list(J     = 8,
                      y     = c(28, 8, -3, 7, -1, 1, 18, 12),
                      sigma = c(15, 10, 16, 11, 9, 11, 10, 18)) 
  # UNIF (8 schools)
  schools_fit_8_UNIF <- stan(model_code = model.unif,
                             data = schools_dat, 
                             warmup = 2000, iter = 6000, chains = 1)
  schools.sim.8.UNIF <- extract(schools_fit_8_UNIF, permuted=TRUE)
  # INVG 1, 1 (8 schools)
  schools_fit_8_INV11 <- stan(model_code = model.gamm11,
                              data = schools_dat, 
                              warmup = 2000, iter = 6000, chains = 1)
  schools.sim.8.INV11 <- extract(schools_fit_8_INV11, permuted=TRUE)
  # INVG .001, .001 (8 schools)
  schools_fit_8_INV00 <- stan(model_code = model.gamm00,
                              data = schools_dat,
                              warmup = 2000, iter = 6000, chains = 1)
  schools.sim.8.INV00 <- extract(schools_fit_8_INV00, permuted=TRUE)
  
  # Point Estimates
  data.frame(mean =   colMeans(cbind(schools.sim.8.UNIF$tau,
                                     schools.sim.8.INV11$tau,
                                     schools.sim.8.INV00$tau)),
             median = apply((cbind(schools.sim.8.UNIF$tau,
                                   schools.sim.8.INV11$tau,
                                   schools.sim.8.INV00$tau)), 2, median)
  )

  # Plots
  par(mfrow=c(3,1))
  hist(schools.sim.8.UNIF$tau,
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "sigma.theta",
       main = "3 Schools: sigma_a posterior (uniform prior on sigma_a)",
       breaks = 20)
  hist(schools.sim.8.INV11$tau,
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "sigma.theta",
       main = "3 Schools: sigma_a posterior (uniform prior on sigma_a)",
       breaks = 50)
  hist(schools.sim.8.INV00$tau,
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "sigma.theta",
       main = "3 Schools: sigma_a posterior (uniform prior on sigma_a)",
       breaks = 50)
  
  #> 3 Schools Example ####
  schools_dat_3 <- list(J     = 3,
                      y     = c(28, 8, -3),
                      sigma = c(15, 10, 16)) 
  # Unif (3 schools)
  schools_fit_3_UNIF <- stan(model_code = model.unif,
                             data = schools_dat_3,
                             warmup = 1000, iter = 6000, chains = 4)
  schools.sim.3.UNI <- extract(schools_fit_3_UNIF, permuted=TRUE)
  
  # Half Cauchy
  schools_fit_3_HC <- stan(model_code = model.halfC,
                           data = schools_dat_3,
                           warmup = 1000, iter = 6000, chains = 4,
                           control = list(adapt_delta = 0.8)
                           )
  schools.sim.3.HC <- extract(schools_fit_3_HC, permuted=TRUE)
  
  # Plots
  par(mfrow=c(3,1))
  hist(schools.sim.3.UNI$tau,
       yaxt = "n", ylab = ".", xlim = c(0, 200), xlab = "sigma.theta",
       main = "3 Schools: sigma_a posterior (uniform prior on sigma_a)",
       breaks = 26)
  hist(schools.sim.3.HC$tau,
       yaxt = "n", ylab = ".", xlim = c(0, 200), xlab = "sigma.theta",
       main = "3 Schools: sigma_a posterior (half-Cauchy (0, 25) prior on sigma_a)",
       breaks = 50)
  pairs(schools_fit_3_HC)
  
# Scraps ####
  # # Fit Model
  # schools_fit_UNIF <- stan(model_code = model.unif,
  #                          data = schools_dat, iter = 2000, chains = 4)
  # print(schools_fit_UNIF)
  # plot(schools_fit_UNIF)
  # # Accessing the posterior simulations in R
  # schools.sim.UNI <- extract(schools_fit_UNIF, permuted=TRUE)
  # attributes(schools.sim.UNI)
  # # Plot
  # hist(schools.sim.UNI$tau)
  # # Posterior probability that the effect is larger in school A than in school C
  # mean(schools.sim.UNI$theta[,1] > schools.sim.UNI$theta[,3])
  # 
  # # Posterior predictive simulations and graphs in R
  # # Replicated data in the existing schools
  # n_sims <- length(schools.sim.UNI$lp__)
  # y_rep <- array (NA, c(n_sims, schools_dat$J))
  # for (s in 1:n_sims){
  #   y_rep[s,] <- rnorm (schools_dat$J, schools.sim.UNI$theta[s,], schools_dat$sigma)
  # }
  # # Plot posterior predictive check
  # par (mfrow=c(5,4), mar=c(4,4,2,2))
  # hist (schools_dat$y, xlab = "", main = "y")
  # for (s in 1:19){
  #   hist (y_rep[s,], xlab="", main=paste("y_rep",s))
  # }
  # #upper-left histogram displays the observed data, and the other 19 histograms
  # #are posterior predictive replications