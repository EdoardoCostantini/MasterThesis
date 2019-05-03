### Project:     Master Thesis
### Object:      Prior Definition
### Description: Contains the code I used to define the prior guesses. You can source this code to get the lists cointaining
###              the priors of interest in your environment. The call to the data is kept in the comments in case one wants
###              to run this code independently. However, the data will usually already be in the environment.
### Date:        2019-03-30

# Normal model:
# y_{ij} = normal(p_{ij})
# y_{ij} = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)

# # Riesbt dat
#   RiesbyDat <- read.table("./data/RiesbyDat.txt")
#   dat_ALL <- data.frame(cluster = RiesbyDat$id,    #constant dataset w/ naming that fits the loop
#                         yvec    = RiesbyDat$depr,  #you only need to change the variables included here
#                         xvec    = RiesbyDat$week,  #and everything will be adapted in the loop
#                         cvec    = RiesbyDat$endog,
#                         inter   = RiesbyDat$inter)
#   dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
#   dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition

library(lme4)

  dimens <- 2   # k
  
# Educated Guess ----------------------------------------------------------
  if(which.model == "normal_rep"){
    # Normal Model repeated measures (use this to define the educated guessmodel)
    # Define the educated guess by fitting a regression line for each group and getting
    # the variance of the parameters (intercept and slope)
    intercept <- rep(NA, dat_ALLn)
    slope <- rep(NA, dat_ALLn)
    fit <- vector("list", dat_ALLn)
    i <- 1
    for (id in 1:dat_ALLn) {
      #print(paste("Before:", i))
      fit[[id]] <- lm(dat_ALL$yvec[i:(i+5)] ~ scale(dat_ALL$xvec[i:(i+5)]))
      intercept[id] <- coef(fit[[id]])[1]
      slope[id]     <- coef(fit[[id]])[2]
      
      i <- i+6
      #print(paste("After:", i))
    }
    vi <- round(var(intercept), 0) # 21 # these guesses are based on the entire sample
    vs <- round(var(slope), 0)     # 9
    
    plot(unique(dat_ALL$xvec), predict(fit[[1]]),
         main = "Depression Scale data: Within-group regression lines (normal model)",
         xlim = c(0, 5), ylim = c(0, 40),
         ylab = "depression", xlab = "time",
         type = "l")
    for (id in 2:dat_ALLn) {
      lines(unique(dat_ALL$xvec), predict(fit[[id]]))
    }
    
    B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
  }
  
  if(which.model == "normal_cross"){
    # Normal Model repeated measures (use this to define the educated guessmodel)
    # Define the educated guess by fitting a regression line for each group and getting
    # the variance of the parameters (intercept and slope)
    intercept <- rep(NA, dat_ALLn)
    slope <- rep(NA, dat_ALLn)
    fit <- vector("list", dat_ALLn)
    for (i in 1:dat_ALLn) {
      current_dat <- hsb[hsb$id == unique(hsb$id)[i], ]
      
      fit[[i]] <- lm(mathach ~ ses, data = current_dat)
      intercept[i] <- coef(fit[[i]])[1]
      slope[i]     <- coef(fit[[i]])[2]
    }
    vi <- round(var(intercept), 0) # 21 # these guesses are based on the entire sample
    vs <- round(var(slope), 0)     # 9
    
    plot(hsb[hsb$id == unique(hsb$id)[1], "ses"], predict(fit[[1]]),
         main = "Math Achivement data: Within-group regression lines (normal model)",
         xlim = c(-3, 3), ylim = c(0, 25),
         ylab = "Math Achivement", xlab = "individual SES",
         type = "l")
    for (i in 2:dat_ALLn) {
      lines(hsb[hsb$id == unique(hsb$id)[i], "ses"], predict(fit[[i]]))
    }
    
    B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
  }
  
  if(which.model == "logistic"){
    # Logistic Model
    # Use this code to find the logistic model educated guess
    intercept <- rep(NA, dat_ALLn)
    slope <- rep(NA, dat_ALLn)
    fit <- vector("list", dat_ALLn)
    i <- 1
    for (id in 1:dat_ALLn) {
      fit[[id]] <- glm(dat_ALL$yvec[i:(i+5)] ~ scale(dat_ALL$xvec[i:(i+5)]), family = binomial)
      intercept[id] <- coef(fit[[id]])[1]
      slope[id]     <- coef(fit[[id]])[2]
      i <- i+6
    }
    vi <- round(var(intercept[abs(intercept) < 5]), 2)
    vs <- round(var(slope[abs(slope) < 20]), 2)
    # Plot all lines
    plot(unique(dat_ALL$xvec), predict(fit[[1]]),
         main = "GPA data: Within-group regression lines (logit scale)",
         xlim = c(-1.5, 1.5), ylim = c(-10, 10),
         ylab = "logit", xlab = "time",
         type = "l")
    for (id in 2:dat_ALLn) {
      lines(unique(dat_ALL$xvec), predict(fit[[id]]))
    }

    B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
  }
  
# Inverse Wishart Priors --------------------------------------------------
# IW uninformative
  IW_PR <- NULL
  IW_PR[[length(IW_PR)+1]] <- list(nu=dimens-1+1,    e=1,    S0=B0_ed) # Educated prior guess (computed above) with a given e
  IW_PR[[length(IW_PR)+1]] <- list(nu=dimens-1+1,    e=1,    S0=(dimens - 2 + 1)*diag(2))
  IW_PR[[length(IW_PR)+1]] <- list(nu=dimens-1+.1,   e=.1,   S0=(dimens - 2 + .1)*diag(2))
  IW_PR[[length(IW_PR)+1]] <- list(nu=dimens-1+.001, e=.001, S0=(dimens - 2 + .001)*diag(2))
  IW_PR[[length(IW_PR)+1]] <- list(nu=2,             e=1,    S0="R*") # R* from Kass Natarajan
  
  # Define names
  names(IW_PR) <- c("IW_eg","IW_e1", "IW_e01", "IW_e0001", "IW_R*")
  
  # Check the structure of the list
  str(IW_PR)


# Matrix F ----------------------------------------------------------------
  MF_PR <- NULL
  MF_PR[[length(MF_PR)+1]] <- list(nu=2,     d=1,  e=1,  S0=10**3*diag(2)) # mat-f proper neighbour guess
  MF_PR[[length(MF_PR)+1]] <- list(nu=2,     d=1,  e=1,  S0=B0_ed)         # Educated guess e = 1
  MF_PR[[length(MF_PR)+1]] <- list(nu=1.1,   d=.1, e=.1, S0=B0_ed)         # Educated guess e =.1
  MF_PR[[length(MF_PR)+1]] <- list(nu=2,     d=1,  e=0,  S0="R*")          # R* [actual prior is computed in the sampling function when S0 == R* is true]
  MF_PR[[length(MF_PR)+1]] <- list(nu=1.001, d=1,  e=1,  S0=10**3*diag(2)) # flat: mat-f proper neighbour to p(sigma) prop to 1
  # MF_PR[[length(MF_PR)+1]] <- list(nu=1.5,d=.5,e=.5, S0=B0_ed)           # plug in educated guess
  # MF_PR[[length(MF_PR)+1]] <- list(nu=1.01,d=.01,e=.01,S0=B0_ed)         # plug in educated guess
  
  # Define names
  names(MF_PR) <- c("MF_pn", "MF_e=1","MF_e=0.1","MF_R*","MF_Flat")
  
  #Check structure
  str(MF_PR)
  
  
  
  

# EXTRA -------------------------------------------------------------------

# # Studing Matrix F prior 
# # Here you can change the values of epsilon, nu and delta to see how they effect
# # a matrix F distribution
# df_freeb_SD <- function(x, nu, d, b){
#   ( 2*gamma( (d+nu)/2 ) / ( gamma(nu/2)*gamma(d/2)*b**(nu/2) ) ) * (x)**(nu - 1) * (1 + (x**2)/b)**(-(nu+d)/2)
# }
# k <- 2 # number of random effects
# epsilon <- 1 # small quantity
# nu <- k - 1 + epsilon
# delta <- epsilon
# 
# # Sampling from matrix-F using full conditonals
# reps <- 10000
# Omega <- B <- matrix(c(21, 0, 0, 9), ncol = 2) # guess based on data exploration and knowledge
# draws <- matrix(rep(NA, reps*4), ncol = 4)
# for (i in 1:reps) {
#   # Psi|Omega,.
#     ScaleMatrix = Omega
#   PsiInvDraw = rwish(v = delta + k - 1,
#                      S = solve(ScaleMatrix))
#   # Omega|Psi,.
#     ScaleOmega = (PsiInvDraw + solve(B))
#   Omega = rwish(v = nu + delta + k - 1,
#                 S = solve(ScaleOmega))
#   draws[i, ] <- c(solve(PsiInvDraw))
# }
# # For precise density
# x <- seq(0, 50, by = .1)
# # Sampling density
# par(mfcol=c(2,3))
# prior_draws_sd <- sqrt(draws[, 1]) #intercept sd
# # h <- hist(prior_draws_sd, breaks = 1e5, xlim = c(0,100))
# #   plot(h$mids, h$density, type = "l", xlim = c(0,50),ylim=c(0,.5), col = "Gray", lty = 1)
# #   lines(x, df_freeb_SD(x, nu = nu, d = delta, b = 21), type = "l", lty = 2)
# plot(x, df_freeb_SD(x, nu = nu, d = delta, b = 21), ylim=c(0,.5), type = "l", lty = 2)
# # prior_draws_sd <- sqrt(draws[, 4]) #slope sd
# # h <- hist(prior_draws_sd, breaks = 1e5, xlim = c(0,100))
# #   plot(h$mids, h$density, type = "l", xlim = c(0,50),ylim=c(0,.5), col = "Gray", lty = 1)
# #   lines(x, df_freeb_SD(x, nu = nu, d = delta, b = 9), type = "l", lty = 2)
# plot(x, df_freeb_SD(x, nu = nu, d = delta, b = 9), ylim=c(0,.5), type = "l", lty = 2)
#   # correlation
# prior_draws_corr <- draws[, 2] / (sqrt(draws[, 1]) * sqrt(draws[, 4]))
# h <- hist(prior_draws_corr, breaks = 10, xlim = c(-1,1))
#   plot(h$mids, h$density, type = "l", col = "Gray", lty = 1)
#   # Cannot plot prior for correaltion for small e (less than 1)
#   # as the wishart sampling function works only for v >= k

# param <- .001
# plot(sdseq, dinvgamma(sdseq, shape = param, scale = param), type = "l", xlim = c(0,20), ylim = c(0, 2))