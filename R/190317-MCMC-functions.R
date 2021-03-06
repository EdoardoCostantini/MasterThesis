### Project:     Master Thesis
### Object:      MCMC functions for model fitting
### Description: Contains the functions used to fit all the models of interest. The scripts to actually fit the model is in a different file
###              Comments in the beginning could be used to run this script independently.
### Date:        2019-03-17

# Normal model:
# y_{ij} = normal(p_{ij})
# y_{ij} = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)

# # Set up
#   post_draws_fun_filename <- "./R/190313-normalmod-functions.R"
#   source(post_draws_fun_filename)
# 
# Data
# Depression Scale (n = 46, J = 6)
# RiesbyDat <- read.table("./data/RiesbyDat.txt")
# 
# # Make it digestable for the functions
# yvec     <- RiesbyDat$depr
# Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
# J        <- length(unique(RiesbyDat$week))
# n        <- nrow(RiesbyDat)/J
# Zi       <- cbind(rep(1, length = J), unique(Xmat[, 2]))
# 
# yvec <- dat_yvec     <- dat$yvec
# Xmat <- dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
# J    <- dat_J        <- length(unique(dat$xvec))
# n    <- dat_n        <- nrow(dat)/dat_J
# Zi   <- dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))
# 
# 
# # Test Functions
#   out1 <- MCMC_invWish(yvec = RiesbyDat$depr,
#                        Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
#                        J    = length(unique(RiesbyDat$week)),
#                        n    = nrow(RiesbyDat)/J,
#                        Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
#                        iniv = 0,
#                        B0   = list(nu=2,e=1,S0=diag(2)),
#                        samsize = 2000,
#                        burnin = 1/10)
#   out2 <- MCMC_HWprior(yvec = RiesbyDat$depr,
#                        Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
#                        J    = length(unique(RiesbyDat$week)),
#                        n    = nrow(RiesbyDat)/J,
#                        Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
#                        iniv = 0,
#                        samsize = 2000,
#                        burnin = 1/10)
#   out1 <- MCMC_matF(yvec = RiesbyDat$depr,
#                     Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
#                     J    = length(unique(RiesbyDat$week)),
#                     n    = nrow(RiesbyDat)/J,
#                     Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
#                     iniv = 1,
#                     B0   = 1e3*diag(2),
#                     samsize = 2000,
#                     burnin = 1/10)
# 
# # Results Exploration
#   # Output selection
#   which_o <- 5 # which Object in list
#   which_c <- 1 # which Column in object
#   # Traceplot
#     plot(1:samsize,out1[[which_o]][,which_c],"l")
#   # Estimates
#     mean(out1[[which_o]][,which_c])
#     median(out1[[which_o]][,which_c])
#     # Compare with lme results
#     fit
#     colMeans(ranef(fit)$cluster)
#   # Posteriors
#     hist(out1[[which_o]][,which_c], breaks = 100,
#          xlab = "Intercept Variance")
#     abline(v = mean(out1[[which_o]][,which_c]), col = "blue", lwd = 2)
#     abline(v = median(out1[[which_o]][,which_c]), col = "red", lwd = 2)

# MCMC Functions ----------------------------------------------------------

# Model fitting with inv-Wish standard prior
# Logic: given the data (in a particular, see above) you can define the initial value option (1 = lme result, 0 = fixed)
# "B0 = " arguments requires a list containing the value of epsilon and the scale matrix to be used in the full conditional
# of Psi. Also can define differen samsize and burn-in period (=0 if you do not want it)

MCMC_invWish = function(yvec, Xmat, Zi, J, n, iniv = 1, B0 = list(e=1,S0=diag(2)), samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec    <- Xmat[, 2]
  covar   <- Xmat[, 3]
  inter   <- Xmat[, 4]
  cluster <- rep(seq(1, n), each = J)
  # Define priors
  bi0   <- rep(0, ncol(Zi)) # prior for bi
  e     <- B0$e # small quantity to add to prior nu
  # Define priors
  S0 <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    lm_fit    <- lm(yvec ~ xvec + covar + inter)
    weightMat <- sigma(lm_fit)^(-1) * diag(J)
    Rstar     <- solve(t(Zi)%*%weightMat%*%Zi)
    S0        <- 2*Rstar # should this be k * Rstar?
    S0Inv     <- solve(S0)
    # Derivation of the weight vector is in your notes on Natarajan and Kass 2018
  }
  # Define initial Values
  if (iniv == 1){
    fit   <- lmer(yvec ~ xvec + covar + inter + (1 + xvec | cluster), REML = FALSE)
    Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + .01*diag(2)# as in MulderPericchi2018 code
    PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
    sigma2 <- attr(VarCorr(fit), "sc")
    bMat   <- as.matrix(ranef(fit)$cluster)
  } else {
    Psi    <- matrix(c(1,0,0,1), ncol = 2)
    PsiInv <- solve(Psi)
    sigma2 <- 1
    bMat   <- matrix(rep(0, n*2), ncol = 2)
  }
  
  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
  colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr
  colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    PsiInv   <- draw_PsiInv_InvWish(n,bMat,S0 = S0, e = e)
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    PsiInv   <- draw_PsiInv_InvWish(n,bMat,S0 = S0, e = e)
    
    PD_theta[ss,] = theta
    PD_bMat[ss,]  = c(bMat)
    PD_Psi[ss,]   = c(solve(PsiInv))
    sdRI    <- sqrt(PD_Psi[ss,1])       # save the satandard deviation version
    sdRS    <- sqrt(PD_Psi[ss,4])
    corRIRS <- PD_Psi[ss,2]/(sdRI*sdRS) # save the correlation
    PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
  }
  B0$S0 <- S0
  out <- list(PD_theta = PD_theta,  #1
              PD_bMat = PD_bMat,   #2
              PD_Psi = PD_Psi,    #3
              PD_Psi_sd = PD_Psi_sd,
              hypepar = B0) #4
  return(out)
}
    
# Model fitting with inv-Wish standard prior
# Logic: given the data (in a particular, see above) you can define the initial value option (1 = lme result, 0 = fixed)
# "B0 = " arguments requires a list containing the value of nu, d, epsilon, and the scale matrix to be used in the full conditional
# of Psi. Also can define differen samsize and burn-in period (=0 if you do not want it)

MCMC_matF = function(yvec, Xmat, Zi, J, n, iniv = 1, B0 = list(nu=2,d=1,e=0,S0=1e3*diag(2)), samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec    <- Xmat[, 2]
  covar   <- Xmat[, 3]
  inter   <- Xmat[, 4]
  cluster <- rep(seq(1, n), each = J)
  # Define priors
  bi0   <- rep(0, ncol(Zi)) # prior for bi
  S0 <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    lm_fit    <- lm(yvec ~ xvec + covar + inter)
    #betaHat   <- c(coefficients(lm_fit))
    weightMat <- sigma(lm_fit)^(-1) * diag(J)
    Rstar     <- solve(t(Zi)%*%weightMat%*%Zi)
    S0        <- 2*Rstar # should this be k * Rstar?
    S0Inv     <- solve(S0)
    # Derivation of the weight vector is in your notes on Natarajan and Kass 2018
  }
  # nu, and d are provided in order by the prior list
  # the first three prior objects have nu = k - 1 + e, d = e, for e = c(1, .5, .1)
  # The function specificaion below follows this priors because when the prior object 
  # is nu = 1.5, d = e and e = .5, we have the prior nu = k -1 + e and delta = e for e = .5
  
  # Define initial Values
  if (iniv == 1){
    fit   <- lmer(yvec ~ xvec + covar + inter + (1 + xvec | cluster), REML = FALSE)
    #theta  <- fixef(fit)
    Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)# as in MulderPericchi2018 code
    PsiInv <- solve(Psi + 1e-6*diag(2))        
    sigma2 <- attr(VarCorr(fit), "sc")
    bMat   <- as.matrix(ranef(fit)$cluster)
  } else {
    #theta  <- t(rep(0, ncol(Xmat)))
    Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
    PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
    sigma2 <- 1 #attr(VarCorr(fit), "sc")
    bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
  }
  Omega = S0
  
  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Omg    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
  colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr
  colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv = S0Inv,n,d=B0$d,nu=B0$nu)
    PsiInv <- outDummy[[1]]
    Omega  <- outDummy[[2]]
    
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv = S0Inv,n,d=B0$d,nu=B0$nu)
    PsiInv <- outDummy[[1]]
    Omega  <- outDummy[[2]]
    
    PD_theta[ss,] = theta
    PD_bMat[ss,]  = c(bMat)
    PD_Psi[ss,]   = c(solve(PsiInv))
    PD_Omg[ss,]   = Omega
    sdRI    <- sqrt(PD_Psi[ss,1])       # save the satandard deviation version
    sdRS    <- sqrt(PD_Psi[ss,4])
    corRIRS <- PD_Psi[ss,2]/(sdRI*sdRS) # save the correlation
    PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
  }
  B0$S0 <- S0
  out <- list(PD_theta=PD_theta,
              PD_bMat=PD_bMat,
              PD_Psi=PD_Psi,
              PD_Psi_sd=PD_Psi_sd,
              PD_Omg=PD_Omg,
              hyperpar=B0)
  return(out)
}
# 
# plot(1:2000, out$PD_Psi_sd[, 1],"l",
#          xlim = c(0, 2000))
# max(out$PD_Psi[, 3])
# plot(density(out$PD_Psi_sd[, 1]))
    
# Model fitting with HW prior
# Logic: given the data (in a particular, see above) you can define the initial value option (1 = lme result, 0 = fixed)
# Compared to the previous functions, there is no need to specify the prior values. They are defined in the functions as
# I'm working with only one version of this prior

MCMC_HWprior = function(yvec, Xmat, Zi, J, n, iniv = 1, samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec    <- Xmat[, 2]
  covar   <- Xmat[, 3]
  inter   <- Xmat[, 4]
  cluster <- rep(seq(1, n), each = J)
  # Define priors
  bi0   <- rep(0, ncol(Zi))
  # Define initial Values
  if (iniv == 1){
    fit   <- lmer(yvec ~ xvec + covar + inter + (1 + xvec | cluster), REML = FALSE)
    #theta  <- fixef(fit)
    Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)# as in MulderPericchi2018 code
    PsiInv <- solve(Psi + 1e-6*diag(2))  
    sigma2 <- attr(VarCorr(fit), "sc")
    bMat   <- as.matrix(ranef(fit)$cluster)
  } else {
    Psi    <- matrix(c(1,0,0,1), ncol = 2)
    PsiInv <- solve(Psi)                   
    sigma2 <- 1
    bMat   <- matrix(rep(0, n*2), ncol = 2)
  }
  avec   <- rep(100,2)
  
  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
  colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr
  colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  PD_avec   = matrix(0, nrow = samsize, ncol = 2)
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
      PsiInv <- outDummy[[1]]
      avec   <- outDummy[[2]]
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
      PsiInv <- outDummy[[1]]
      avec   <- outDummy[[2]]
    
    PD_theta[ss,] = theta
    PD_bMat[ss,]  = c(bMat)
    PD_avec[ss,]  = avec
    PD_Psi[ss,]   = c(solve(PsiInv))
    sdRI    <- sqrt(PD_Psi[ss,1])       # save the satandard deviation version
    sdRS    <- sqrt(PD_Psi[ss,4])
    corRIRS <- PD_Psi[ss,2]/(sdRI*sdRS) # save the correlation
    PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
  }
  out <- list(PD_theta,  #1
              PD_bMat,   #2
              PD_Psi,    #3
              PD_Psi_sd, #4
              PD_avec)
  names(out) <- c("PD_theta", "PD_bMat", "PD_Psi", "PD_Psi_sd", "PD_avec")
  return(out)
}
    