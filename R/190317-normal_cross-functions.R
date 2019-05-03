### Project:     Master Thesis
### Object:      MCMC functions for model fitting w/ crossectional HSB data
### Description: Contains all functions and model fitting for the normal 
###              model to corss sectional data
### Date:        2019-04-29

# Normal model:
# y_{ij} = normal(p_{ij})
# y_{ij} = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)

library(lme4)
library(mvtnorm)
library(MCMCpack)
library(ddpcr)
  
# Sampling Functions ------------------------------------------------------

draw_theta = function(yvec,Xmat,Zi,bMat,sigma2,n){
  
  bizij <- NULL
  for (i in 1:n) {
    bizij <- c(bizij, Zi[[i]] %*% bMat[i,])
  }
  ytilde    <- yvec - bizij
  # thetaCovm <- solve(solve(Sigma0) + t(Xmat)%*%Xmat/sigma2)
  # thetaMean <- thetaCovm %*% ( solve(Sigma0)%*%theta0 + t(Xmat)%*%ytilde/sigma2)
  thetaCovm <- solve(t(Xmat)%*%Xmat/sigma2)
  thetaMean <- thetaCovm %*% (t(Xmat)%*%ytilde/sigma2)
  thetaDraw <- rmvnorm(1, thetaMean, thetaCovm)
  return(t(thetaDraw))
  
}
  
draw_bMat = function(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n){
  
  ytilde = yvec - c(Xmat%*%theta)
  for(i in 1:n){
    ytilde_i <- ytilde[cluster == unique(cluster)[i]]
    bi_Covm  <- solve(PsiInv + t(Zi[[i]])%*%Zi[[i]]/sigma2)
    bi_mean  <- bi_Covm %*% (PsiInv %*% bi0 + t(Zi[[i]])%*%ytilde_i/sigma2)
    bMat[i,] <- c(rmvnorm(1, bi_mean, bi_Covm))
  }
  return(bMat)
  
}
  
draw_sigam2_IMPprior = function(yvec,Xmat,Zi,theta,bMat,n,J){
  # sigma^(-2) prior
  
  bizij <- NULL
  for (i in 1:n) {
    bizij <- c(bizij, Zi[[i]] %*% bMat[i,])
  }
  
  SSR <- t(yvec - Xmat %*% theta - bizij) %*% (yvec - Xmat %*% theta - bizij)
  sigma2Draw <- rinvgamma(1, (length(yvec))/2, (SSR)/2)
  return(sigma2Draw)
  
}

draw_sigam2_IGprior = function(yvec,Xmat,Zi,theta,bMat,n,J,nu0=.001,sigma20=.001){
  # IG(nu0/2, nu0*S20/2) prior
  
  bizij <- NULL
  for (i in 1:n) {
    bizij <- c(bizij, Zi[[i]] %*% bMat[i,])
  }
  
  SSR <- t(yvec - Xmat %*% theta - bizij) %*% (yvec - Xmat %*% theta - bizij)
  sigma2Draw <- rinvgamma(1, (nu0 + length(yvec))/2, (nu0*sigma20 + SSR)/2)
  return(sigma2Draw)
  
}

draw_PsiInv_InvWish = function(n,bMat,S0=diag(2),e=1){
  # Prior definition
    k <- ncol(bMat)
    nu0 <- k - 1 + e
  # Posterior samples
    ScaleMatrix = t(bMat)%*%bMat + S0
    PsiInvDraw = rwish(v = n + nu0,             # nu0 = 2
                                                # usually nu0 = 2, v = n + 2; needs to be nu = k - 1 + e
                       S = solve(ScaleMatrix))  # distirbution is just Wishart, not inverse! Therefore, the guess priovaded is invterted!
  return(PsiInvDraw)
}

draw_PsiInv_matF = function(bMat,PsiInv,Omega,B0Inv,n,d=1,nu=2){
    k  <- ncol(bMat)
  # Psi|Omega,.s
      ScaleMatrix = t(bMat)%*%bMat + Omega + 1e-6*diag(2) # current data estimation of random effect covariance matrix + Omega = Previuos draw (or initial guess)
    PsiInvDraw = rwish(v = (n) + (d + k - 1),             # (n) + (d + k - 1) = (posterior) + (prior contribution) [original: n + 2]
                       S = solve(ScaleMatrix)) + 1e-6*diag(2)
    # It's an inverse draw because I would like to get a draw from IW(v, S),
    # but I get a draw from W(v, solve(S)), which becomes a draw from the IW(v, S), when I take its inverse
    # IW(v, S) <=> solve( draw from W(v, solve(S)) )
  # Omega|Psi,.
      ScaleOmega = solve(PsiInvDraw + B0Inv)
    Omega = rwish(v = nu + d + k - 1,           # nu + d + k - 1 = 4 when? (all prior: no information directly available on Omega) [original: 4]
                  S = ScaleOmega)
  return(list(PsiInvDraw,Omega))
}

draw_PsiInv_HW = function(PsiInv,avec,bMat,n,nu=2,eta=1/2,Ak=10**5){
  # Prior Parameters
    # nu  <- 2; eta <- 1/2; Ak  <- 10**5
    k   <- ncol(bMat)
  # Posterior samples
    ScaleMatrix = t(bMat)%*%bMat +               # sum of ui %*% ui'
                  2*nu*diag(1/avec)              
    PsiInvDraw = rwish(v = n + nu + k - 1,       # e.g.: nu = 2, k = 2, n = 30 (number of clusters), 
                       S = solve(ScaleMatrix))
    
    scale_avec = 2*diag(PsiInvDraw) + 1/Ak**2    # a_k scale parameter in Haung Wand 2013 (section 4 full conditionals)
    avec = rinvgamma(2,                          # k = 2: 1 random intercept, 1 random slope
                     shape = eta*(nu + k),       # (nu + k) /2, where k = 2, nu = 2
                     scale = scale_avec)
  return(list(PsiInvDraw,avec))
}

# MCMC Functions ----------------------------------------------------------
# > Crossectional data - Raudenbush Byrk model ############################

MCMC_cross_invWish = function(yvec, Xmat, Zi, cluster, J, n, B0 = list(e=1,S0=diag(2)), samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec <- Xmat[, 2]
  lvl2_cov1 <- Xmat[, 3]
  lvl2_cov2 <- Xmat[, 4]
  int_xvec_lvl2_cov1 <- Xmat[, 5]
  int_xvec_lvl2_cov2 <- Xmat[, 6]
  
  # Define priors
  bi0   <- rep(0, ncol(Zi[[1]])) # prior for bi
  e     <- B0$e # small quantity to add to prior nu
  S0 <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    Rstar <- vector("list", n)
    lm_fit    <- lm(yvec ~ xvec + lvl2_cov1 + lvl2_cov2 + int_xvec_lvl2_cov1 + int_xvec_lvl2_cov2)
    for (i in 1:n) {
      weightMat <- sigma(lm_fit)^(-1) * diag(J[i])
      Rstar[[i]]     <- solve(t(Zi[[i]])%*%weightMat%*%Zi[[i]])
    }
    S0 <- 2*Reduce("+", Rstar)/n
    S0Inv     <- solve(S0)
  }
  
  # Define initial Values
  fit    <- lmer(yvec ~ xvec*lvl2_cov1 + xvec*lvl2_cov2 + (1 + xvec | cluster), REML = T)
  Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + 1e-6*diag(2)# as in MulderPericchi2018 code
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$cluster)

  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS var and cov
    colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS sd and corr
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    #sigma2   <- draw_sigam2_IGprior(yvec,Xmat,Zi,theta,bMat,n,J) # 
    PsiInv   <- draw_PsiInv_InvWish(n,bMat,S0 = S0, e = e)
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
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
  out <- list(PD_theta  = PD_theta,  #1
              PD_bMat   = PD_bMat,   #2
              PD_Psi    = PD_Psi,    #3
              PD_Psi_sd = PD_Psi_sd,
              hypepar   = B0) #4
  return(out)
}

MCMC_cross_matF = function(yvec, Xmat, Zi, J, n, cluster, B0 = list(nu=2,d=1,e=0,S0=1e3*diag(2)), samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec <- Xmat[, 2]
  lvl2_cov1 <- Xmat[, 3]
  lvl2_cov2 <- Xmat[, 4]
  int_xvec_lvl2_cov1 <- Xmat[, 5]
  int_xvec_lvl2_cov2 <- Xmat[, 6]
  
  # Define priors
  bi0   <- rep(0, ncol(Zi[[1]])) # prior for bi
  e     <- B0$e # small quantity to add to prior nu
  S0 <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    Rstar <- vector("list", n)
    lm_fit    <- lm(yvec ~ xvec + lvl2_cov1 + lvl2_cov2 + int_xvec_lvl2_cov1 + int_xvec_lvl2_cov2)
    for (i in 1:n) {
      weightMat <- sigma(lm_fit)^(-1) * diag(J[i])
      Rstar[[i]]     <- solve(t(Zi[[i]])%*%weightMat%*%Zi[[i]])
    }
    S0 <- 2*Reduce("+", Rstar)/n
    S0Inv     <- solve(S0)
  }
  Omega = S0
  
  # Define initial Values
  fit    <- lmer(yvec ~ xvec*lvl2_cov1 + xvec*lvl2_cov2 + (1 + xvec | cluster), REML = T)
  Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + 1e-6*diag(2)# as in MulderPericchi2018 code
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$cluster)
  
  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Omg    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS var and cov
    colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS sd and corr
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
    
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_matF(bMat,PsiInv,Omega,B0Inv = S0Inv,n,d=B0$d,nu=B0$nu)
      PsiInv <- outDummy[[1]]
      Omega  <- outDummy[[2]]
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_matF(bMat,PsiInv,Omega,B0Inv = S0Inv,n,d=B0$d,nu=B0$nu)
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
  out <- list(PD_theta = PD_theta,
              PD_bMat  = PD_bMat,
              PD_Psi   = PD_Psi,
              PD_Psi_sd= PD_Psi_sd,
              PD_Omg   = PD_Omg,
              hyperpar = B0)
  return(out)
}
    
MCMC_cross_HWprior = function(yvec, Xmat, Zi, J, n, cluster, samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec <- Xmat[, 2]
  lvl2_cov1 <- Xmat[, 3]
  lvl2_cov2 <- Xmat[, 4]
  int_xvec_lvl2_cov1 <- Xmat[, 5]
  int_xvec_lvl2_cov2 <- Xmat[, 6]
  
  # Define priors
  bi0   <- rep(0, ncol(Zi[[1]]))
  avec   <- rep(100,2)
  
  # Define initial Values
  fit    <- lmer(yvec ~ xvec*lvl2_cov1 + xvec*lvl2_cov2 + (1 + xvec | cluster), REML = T)
  Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + 1e-6*diag(2)# as in MulderPericchi2018 code
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$cluster)

  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS var and cov
    colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS sd and corr
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  PD_avec   = matrix(0, nrow = samsize, ncol = 2)
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
      PsiInv <- outDummy[[1]]
      avec   <- outDummy[[2]]
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
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

# > Repeated obs data - 1 iv-lvl1, 1 cov-lvl2, 1 inter #####################

MCMC_longi_invWish = function(yvec, Xmat, Zi, cluster, J, n, B0 = list(e=1,S0=diag(2)), samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec <- Xmat[, 2]
  lvl2_cov1 <- Xmat[, 3]
  int_xvec_lvl2_cov1 <- Xmat[, 4]
  
  # Define priors
  bi0   <- rep(0, ncol(Zi[[1]])) # prior for bi
  e     <- B0$e # small quantity to add to prior nu
  S0 <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    Rstar <- vector("list", n)
    lm_fit    <- lm(yvec ~ xvec*lvl2_cov1)
    for (i in 1:n) {
      weightMat <- sigma(lm_fit)^(-1) * diag(J[i])
      Rstar[[i]]     <- solve(t(Zi[[i]])%*%weightMat%*%Zi[[i]])
    }
    S0 <- 2*Reduce("+", Rstar)/n
    S0Inv     <- solve(S0)
  }
  
  # Define initial Values
  fit    <- lmer(yvec ~ xvec*lvl2_cov1 + (1 + xvec | cluster), REML = T)
  Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + 1e-6*diag(2)# as in MulderPericchi2018 code
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$cluster)

  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS var and cov
    colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS sd and corr
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n=n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n=n,J=J)
    #sigma2   <- draw_sigam2_IGprior(yvec,Xmat,Zi,theta,bMat,n,J) # 
    PsiInv   <- draw_PsiInv_InvWish(n,bMat,S0 = S0, e = e)
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n=n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n=n,J=J)
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
  out <- list(PD_theta  = PD_theta,  #1
              PD_bMat   = PD_bMat,   #2
              PD_Psi    = PD_Psi,    #3
              PD_Psi_sd = PD_Psi_sd,
              hypepar   = B0) #4
  return(out)
}

MCMC_longi_matF = function(yvec, Xmat, Zi, J, n, cluster, B0 = list(nu=2,d=1,e=0,S0=1e3*diag(2)), samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec <- Xmat[, 2]
  lvl2_cov1 <- Xmat[, 3]
  int_xvec_lvl2_cov1 <- Xmat[, 4]
  
  # Define priors
  bi0   <- rep(0, ncol(Zi[[1]])) # prior for bi
  e     <- B0$e # small quantity to add to prior nu
  S0 <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    Rstar <- vector("list", n)
    lm_fit    <- lm(yvec ~ xvec*lvl2_cov1)
    for (i in 1:n) {
      weightMat <- sigma(lm_fit)^(-1) * diag(J[i])
      Rstar[[i]]     <- solve(t(Zi[[i]])%*%weightMat%*%Zi[[i]])
    }
    S0 <- 2*Reduce("+", Rstar)/n
    S0Inv     <- solve(S0)
  }
  Omega = S0
  
  # Define initial Values
  fit    <- lmer(yvec ~ xvec*lvl2_cov1 + (1 + xvec | cluster), REML = T)
  Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + 1e-6*diag(2)# as in MulderPericchi2018 code
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$cluster)
  
  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Omg    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS var and cov
    colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS sd and corr
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
    
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_matF(bMat,PsiInv,Omega,B0Inv = S0Inv,n,d=B0$d,nu=B0$nu)
      PsiInv <- outDummy[[1]]
      Omega  <- outDummy[[2]]
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_matF(bMat,PsiInv,Omega,B0Inv = S0Inv,n,d=B0$d,nu=B0$nu)
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
  out <- list(PD_theta = PD_theta,
              PD_bMat  = PD_bMat,
              PD_Psi   = PD_Psi,
              PD_Psi_sd= PD_Psi_sd,
              PD_Omg   = PD_Omg,
              hyperpar = B0)
  return(out)
}

MCMC_longi_HWprior = function(yvec, Xmat, Zi, J, n, cluster, samsize = 2000, burnin = 1/10){
  # Data Prep
  xvec <- Xmat[, 2]
  lvl2_cov1 <- Xmat[, 3]
  int_xvec_lvl2_cov1 <- Xmat[, 4]
  
  # Define priors
  bi0   <- rep(0, ncol(Zi[[1]]))
  avec   <- rep(100,2)
  
  # Define initial Values
  fit    <- lmer(yvec ~ xvec*lvl2_cov1 + (1 + xvec | cluster), REML = T)
  Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2) + 1e-6*diag(2)# as in MulderPericchi2018 code
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$cluster)

  # Run MCMC
  PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
  PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
  PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS var and cov
    colnames(PD_Psi) <- c("RI_var", "cov", "cov", "RS_var")
  PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi[[1]])*2) # store RI, RS sd and corr
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  PD_avec   = matrix(0, nrow = samsize, ncol = 2)
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
    sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
    outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
      PsiInv <- outDummy[[1]]
      avec   <- outDummy[[2]]
  }
  #sampling
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2,n=n)
    bMat     <- draw_bMat(yvec,Xmat,Zi,theta,cluster,bi0,bMat,sigma2,PsiInv,n)
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

# Traceplots
  plot(1:1e4, out$PD_Psi_sd[, 4],"l",
           xlim = c(0, 1e4), ylim = c(0,3))
  max(out$PD_Psi[, 3])
  plot(density(out$PD_Psi_sd[, 1]))

# Check results (posterior and traceplots)
  str(out)
  x <- out[[3]]
  matrix(apply(x, 2, mean), ncol = 2) # very similar to what reported in Raudenbush, then prior, repetitions and summary play a role
  x <- out[[3]]
  par <- 4
  sum(x[,par] < 0)
  # Posterior
  plot(density(x[,par]))
  median(x[,par])
  hist(x[,par],
       xlim = c(0, 20), breaks = 20)
  hist(x[,par],
       xlim = c(-1, 1), breaks = 20)
  # Traceplot
  trax <- out[[3]][,par]
  plot(1:samsize, trax, "l",
       xlim = c(0, samsize), ylim = c(0, 5))
  # Point estimates
  glme4_loop_COND$`n_200:J_6`$glme_Psi_sd
  colMeans(out$PD_Psi_sd)
  
  # Fixed efffects
  x <- out$PD_theta; head(x)
  apply(x, 2, median)
  par <- 2
  # Posterior
  plot(density(x[,par]))
  median(x[,par])
  hist(x[,par],
       xlim = c(0, 20), breaks = 20)

# Model Fitting -----------------------------------------------------------
  
  which.model <- "normal_rep"; source("./R/190330-hyperparameters.R")
  # Sampling Repetitions
  MCMC_reps   <- 1e4
  MCMC_burnin <- 1/10
  
{
  
# > Longitudinal data #######################################################
  # Load Repeated measurements Data
  RiesbyDat <- read.table("./data/RiesbyDat.txt")
  head(RiesbyDat)
  yvec    = RiesbyDat$depr
  xvec    = RiesbyDat$week
  cvec    = RiesbyDat$endog
  inter   = RiesbyDat$inter
  
  dat_ALL <- data.frame(cluster = RiesbyDat$id,    #constant dataset w/ naming that fits the loop
                        yvec    = RiesbyDat$depr,  #you only need to change the variables included here
                        xvec    = RiesbyDat$week,  #and everything will be adapted in the loop
                        cvec    = RiesbyDat$endog,
                        inter   = RiesbyDat$inter)
  dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
  dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition
  # Define conditions
  allCONDs <- list(
    n_cond = list("46" = unique(dat_ALL$cluster),
                  "8" = c(101, 117, 505, 302, 335, 338, 319, 353), # 319 was 350
                  "4" = c(101, 117, 505, 302)), #504, 319, 328, 353))
    J_cond = list("6" = unique(dat_ALL$xvec),
                  "4" = c(1,2,3,4),
                  "3" = c(1,2,3))
  )
  
  conds.index <- matrix(c(46, 6,
                          8, 6,
                          8, 4,
                          8, 3,  # not 2 because w/ lme random effects are notidentifiable
                          4, 6,
                          4, 4,
                          4, 3), # not smaller than 4 because otherwise:
                        ncol = 2, byrow = T)
  
  nconds <- nrow(conds.index)
  
  # LMEr estimates ####
  lme4_loop_COND <- vector("list", length = nconds)
  names(lme4_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # # Standard estimation
    lmefit <- lmer(yvec ~ xvec + cvec + inter + (1 + xvec | cluster), data = dat, REML = FALSE)
    Psi    <- matrix(VarCorr(lmefit)[[1]][1:4], ncol = 2) # as in MulderPericchi2018 code
    Psi_sd <- matrix(c(attributes(VarCorr(lmefit)[[1]])$stddev[1],         # save sd of intercepts
                       rep(attributes(VarCorr(lmefit)[[1]])$correlation[1,2], 2),  # save correlation
                       attributes(VarCorr(lmefit)[[1]])$stddev[2]),         # save sd of slopes
                     ncol = 2)
    sigma2 <- attr(VarCorr(lmefit), "sc")
    bMat   <- as.matrix(ranef(lmefit)$cluster)
    theta  <- fixef(lmefit)
    lme4_loop_COND[[outps]] <- list(lme_fixed  = theta,
                                    lme_bMat   = bMat,
                                    lme_Psi    = Psi,
                                    lme_Psi_sd = Psi_sd,
                                    lme_sigma2 = sigma2)
  }
  str(lme4_loop_COND)
  
  # IW SAMPLING ####
  set.seed(19044)
  # Create objects for storing results of all conditions
  IW_loop_COND <- vector("list", length = nconds)
    names(IW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
    
  # Fit the model looping over all conditions
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Create object for storing resutls of condition loop
    output_loop_prior <- vector("list", length = length(IW_PR))
      names(output_loop_prior) <- names(IW_PR)
      
    # Select data for this condition
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec <- dat$yvec
    dat_Xmat <- matrix(c(rep(1, length(dat$yvec)),
                         dat$xvec,
                         dat$cvec,
                         dat$inter),
                       ncol = 4)
    dat_n  <- length(unique(dat$cluster)) # <- n
    dat_Zi <- vector("list", dat_n)       # <- Zi
    dat_J  <- rep(NA, dat_n)              # <- J
    cluster <- unique(dat$cluster)
    for (i in 1:dat_n) {
      dat_Zi[[i]] <- matrix(c(rep(1, sum(dat$cluster == cluster[i])),
                              dat[dat$cluster == cluster[i], c("xvec")]),
                            ncol = 2)
      dat_J[i] <- sum(dat$cluster == cluster[i])
    }
    dat_subjects <- dat$cluster
    
    # Run MCMC for condition "outps" for all priors specified in the IW_PR
    for (PV in 1:length(IW_PR)) {
      quiet(
        try(
          output_loop_prior[[PV]] <- MCMC_longi_invWish(yvec    = dat_yvec, 
                                                        Xmat    = dat_Xmat, 
                                                        Zi      = dat_Zi,
                                                        cluster = dat_subjects,
                                                        J       = dat_J, 
                                                        n       = dat_n, 
                                                        samsize = MCMC_reps, 
                                                        burnin  = MCMC_burnin, 
                                                        B0 = IW_PR[[PV]]), 
          silent = TRUE
        )
      )
    }
    IW_loop_COND[[outps]] <- output_loop_prior
  }
    str(IW_loop_COND)
  
  # MF SAMPLING ####
  set.seed(19044)
  # Create objects for storing results of all conditions
  MF_loop_COND <- vector("list", length = nconds)
    names(MF_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  # Fit the model looping over all conditions
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Create object for storing resutls obtained with different priors in condition "outps"
    output_loop_prior <- vector("list", length = length(MF_PR))
      names(output_loop_prior) <- names(MF_PR)
      
    # Select data for this condition
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec <- dat$yvec
    dat_Xmat <- matrix(c(rep(1, length(dat$yvec)),
                         dat$xvec,
                         dat$cvec,
                         dat$inter),
                       ncol = 4)
    dat_n  <- length(unique(dat$cluster)) # <- n
    dat_Zi <- vector("list", dat_n)       # <- Zi
    dat_J  <- rep(NA, dat_n)              # <- J
    cluster <- unique(dat$cluster)
    for (i in 1:dat_n) {
      dat_Zi[[i]] <- matrix(c(rep(1, sum(dat$cluster == cluster[i])),
                              dat[dat$cluster == cluster[i], c("xvec")]),
                            ncol = 2)
      dat_J[i] <- sum(dat$cluster == cluster[i])
    }
    dat_subjects <- dat$cluster
    
    # Run MCMC for condition "outps" for all priors specified in the MF_PR
    for (PV in 1:length(MF_PR)) {
      quiet(
        try(
          output_loop_prior[[PV]] <- MCMC_longi_matF(yvec    = dat_yvec, 
                                                     Xmat    = dat_Xmat, 
                                                     Zi      = dat_Zi, 
                                                     J       = dat_J, 
                                                     n       = dat_n,
                                                     cluster = dat_subjects,
                                                     samsize = MCMC_reps, 
                                                     burnin  = MCMC_burnin, 
                                                     B0      = MF_PR[[PV]]), 
          silent = TRUE
        )
      )
    }
    MF_loop_COND[[outps]] <- output_loop_prior
  }
    str(MF_loop_COND)
    
  # HW SAMPLING ####
  set.seed(19044)
  # Create objects for storing results of all conditions
  HW_loop_COND <- vector("list", length = nconds)
    names(HW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  # Fit the model looping over all conditions
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Select data for this condition
    clusters_goal <- allCONDs$n_cond[[which(names(allCONDs$n_cond) == conds.index[outps, 1])]]
    obs_goal      <- allCONDs$J_cond[[which(names(allCONDs$J_cond) == conds.index[outps, 2])]]
    dat           <- dat_ALL[dat_ALL$cluster %in% clusters_goal & dat_ALL$xvec %in% obs_goal, ]
    
    # Define Data for function
    dat_yvec <- dat$yvec
    dat_Xmat <- matrix(c(rep(1, length(dat$yvec)),
                         dat$xvec,
                         dat$cvec,
                         dat$inter),
                       ncol = 4)
    dat_n  <- length(unique(dat$cluster)) # <- n
    dat_Zi <- vector("list", dat_n)       # <- Zi
    dat_J  <- rep(NA, dat_n)              # <- J
    cluster <- unique(dat$cluster)
    for (i in 1:dat_n) {
      dat_Zi[[i]] <- matrix(c(rep(1, sum(dat$cluster == cluster[i])),
                              dat[dat$cluster == cluster[i], c("xvec")]),
                            ncol = 2)
      dat_J[i] <- sum(dat$cluster == cluster[i])
    }
    dat_subjects <- dat$cluster
    
    # Run MCMC for condition "outps"
    quiet(
      try(
        HW_loop_COND[[outps]] <- MCMC_longi_HWprior(yvec    = dat_yvec, 
                                                    Xmat    = dat_Xmat, 
                                                    Zi      = dat_Zi, 
                                                    J       = dat_J, 
                                                    n       = dat_n, 
                                                    cluster = dat_subjects,
                                                    samsize = MCMC_reps, 
                                                    burnin  = MCMC_burnin), 
        silent = TRUE
      )
    )
     
  }
    str(HW_loop_COND)

  # Save output
  output <- list(out_lme4 = lme4_loop_COND,
                 out_IW   = IW_loop_COND,
                 out_MF   = MF_loop_COND,
                 out_HW   = HW_loop_COND)
  saveRDS(output, paste0("./output/", "riesbydata", Sys.Date(), "-ncond_", nconds, "-rep_", MCMC_reps, ".rds"))
}  
  
# > Cross-sectional data ####################################################
  
# Data (High School and Beyond data)
  hsb <- as.data.frame(haven::read_dta(file = "./Data/hsb.dta"))[, c("id", "mathach", "ses", "sector", "meanses")]; head(hsb)
  
  dat_ALLn <- length(unique(hsb$id))
  which.model <- "normal_cross"; source("./R/190330-hyperparameters.R")
    
  # Define conditions
  n_cond = list("160" = unique(hsb$id),
                  "8"   = c(sample(unique(hsb$id[hsb$sector == 1]), 4),  # sample four chatolic schools
                            sample(unique(hsb$id[hsb$sector == 0]), 4)),
                  "4"   = c(sample(unique(hsb$id[hsb$sector == 1]), 2), 
                            sample(unique(hsb$id[hsb$sector == 0]), 2)))
  conds.index <- matrix(c(160, 45,
                          8, 45,
                          8, 4,
                          8, 3,
                          4, 45,
                          4, 4,
                          4, 3),
                        ncol = 2, byrow = T)
  nconds <- nrow(conds.index)
  
  # Select cases based on these conditons
  dat_list <- vector("list", nconds)
  set.seed(20190503)
  for (outps in 1:nconds) {
    # Select Dataset for condition
    if(conds.index[outps, 2] == 45){
      current_groups <- vector("list", conds.index[outps, 1])
      for (i in 1:length(current_groups)) {
        current_groups[[i]] <- hsb[hsb$id == n_cond[[which(names(n_cond) == conds.index[outps, 1])]][i], ]
      }
      dat_list[[outps]] <- do.call(rbind, current_groups)
    }
    if(conds.index[outps, 2] != 45){
      current_groups <- vector("list", conds.index[outps, 1])
      for (i in 1:length(current_groups)) {
        school_selected <- hsb[hsb$id == n_cond[[which(names(n_cond) == conds.index[outps, 1])]][i], ]
        current_groups[[i]] <- school_selected[sample(seq(1, nrow(school_selected)), conds.index[outps,2]), ]
      }
      dat_list[[outps]] <- do.call(rbind, current_groups)
    }
  }
  
  # MCMC specs
  MCMC_reps   <- 1e4
  MCMC_burnin <- 1/10
  
  # LMEr estimates ####
  lme4_loop_COND <- vector("list", length = nconds)
    names(lme4_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Select Dataset for condition
    dat <- dat_list[[outps]]
    
    # Standard estimation
    lmefit <- lmer(mathach ~ ses + sector*ses + meanses*ses + (1 + ses | id), data = dat, REML = FALSE)
    Psi    <- matrix(VarCorr(lmefit)[[1]][1:4], ncol = 2) # as in MulderPericchi2018 code
    Psi_sd <- matrix(c(attributes(VarCorr(lmefit)[[1]])$stddev[1],         # save sd of intercepts
                       rep(attributes(VarCorr(lmefit)[[1]])$correlation[1,2], 2),  # save correlation
                       attributes(VarCorr(lmefit)[[1]])$stddev[2]),         # save sd of slopes
                     ncol = 2)
    sigma2 <- attr(VarCorr(lmefit), "sc")
    bMat   <- as.matrix(ranef(lmefit)$id)
    theta  <- fixef(lmefit)
    lme4_loop_COND[[outps]] <- list(lme_fixed  = theta,
                                    lme_bMat   = bMat,
                                    lme_Psi    = Psi,
                                    lme_Psi_sd = Psi_sd,
                                    lme_sigma2 = sigma2)
  }
  str(lme4_loop_COND)
  
  # IW SAMPLING ####
  set.seed(19044)
  # Create objects for storing results of all conditions
  IW_loop_COND <- vector("list", length = nconds)
    names(IW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
    
  # Fit the model looping over all conditions
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Create object for storing resutls of condition loop
    output_loop_prior <- vector("list", length = length(IW_PR))
      names(output_loop_prior) <- names(IW_PR)
      
   # Select Dataset for condition
    dat <- dat_list[[outps]]
    
    # Define Data for function
    dat_yvec <- dat$mathach
    dat_Xmat <- matrix(c(rep(1, length(dat$mathach)),
                         dat$ses,
                         dat$sector,
                         dat$meanses,
                         dat$ses*dat$sector,
                         dat$ses*dat$meanses),
                       ncol = 6)
    dat_n  <- length(unique(dat$id)) # <- n
    dat_Zi <- vector("list", dat_n)       # <- Zi
    dat_J  <- rep(NA, dat_n)              # <- J
    cluster <- unique(dat$id)
    for (i in 1:dat_n) {
      dat_Zi[[i]] <- matrix(c(rep(1, sum(dat$id == cluster[i])),
                              dat[dat$id == cluster[i], c("ses")]),
                            ncol = 2)
      dat_J[i] <- sum(dat$id == cluster[i])
    }
    dat_subjects <- dat$id
    
    # Run MCMC for condition "outps" for all priors specified in the IW_PR
    for (PV in 1:length(IW_PR)) {
      quiet(
        try(
          output_loop_prior[[PV]] <- MCMC_cross_invWish(yvec    = dat_yvec, 
                                                        Xmat    = dat_Xmat, 
                                                        Zi      = dat_Zi,
                                                        cluster = dat_subjects,
                                                        J       = dat_J, 
                                                        n       = dat_n, 
                                                        samsize = MCMC_reps, 
                                                        burnin  = MCMC_burnin, 
                                                        B0 = IW_PR[[PV]]), 
          silent = TRUE
        )
      )
    }
    IW_loop_COND[[outps]] <- output_loop_prior
  }
    str(IW_loop_COND)
  
  # MF SAMPLING ####
  set.seed(19044)
  # Create objects for storing results of all conditions
  MF_loop_COND <- vector("list", length = nconds)
    names(MF_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  # Fit the model looping over all conditions
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Create object for storing resutls obtained with different priors in condition "outps"
    output_loop_prior <- vector("list", length = length(MF_PR))
      names(output_loop_prior) <- names(MF_PR)
      
    # Select Dataset for condition
    dat <- dat_list[[outps]]
    
    # Define Data for function
    dat_yvec <- dat$mathach
    dat_Xmat <- matrix(c(rep(1, length(dat$mathach)),
                         dat$ses,
                         dat$sector,
                         dat$meanses,
                         dat$ses*dat$sector,
                         dat$ses*dat$meanses),
                       ncol = 6)
    dat_n  <- length(unique(dat$id)) # <- n
    dat_Zi <- vector("list", dat_n)       # <- Zi
    dat_J  <- rep(NA, dat_n)              # <- J
    cluster <- unique(dat$id)
    for (i in 1:dat_n) {
      dat_Zi[[i]] <- matrix(c(rep(1, sum(dat$id == cluster[i])),
                              dat[dat$id == cluster[i], c("ses")]),
                            ncol = 2)
      dat_J[i] <- sum(dat$id == cluster[i])
    }
    dat_subjects <- dat$id
    
    # Run MCMC for condition "outps" for all priors specified in the MF_PR
    for (PV in 1:length(MF_PR)) {
      quiet(
        try(
          output_loop_prior[[PV]] <- MCMC_cross_matF(yvec    = dat_yvec, 
                                                     Xmat    = dat_Xmat, 
                                                     Zi      = dat_Zi, 
                                                     J       = dat_J, 
                                                     n       = dat_n,
                                                     cluster = dat_subjects,
                                                     samsize = MCMC_reps, 
                                                     burnin  = MCMC_burnin, 
                                                     B0      = MF_PR[[PV]]), 
          silent = TRUE
        )
      )
    }
    MF_loop_COND[[outps]] <- output_loop_prior
  }
    str(MF_loop_COND)
    
  # HW SAMPLING ####
  set.seed(19044)
  # Create objects for storing results of all conditions
  HW_loop_COND <- vector("list", length = nconds)
    names(HW_loop_COND) <- c(paste0("n_", conds.index[, 1], ":J_", conds.index[, 2]))
  # Fit the model looping over all conditions
  for (outps in 1:nconds) {
    print(paste("Condtion", seq(1, nconds)[outps], "Started at:", format(Sys.time())))
    
    # Select Dataset for condition
    dat <- dat_list[[outps]]
    
    # Define Data for function
    dat_yvec <- dat$mathach
    dat_Xmat <- matrix(c(rep(1, length(dat$mathach)),
                         dat$ses,
                         dat$sector,
                         dat$meanses,
                         dat$ses*dat$sector,
                         dat$ses*dat$meanses),
                       ncol = 6)
    dat_n  <- length(unique(dat$id)) # <- n
    dat_Zi <- vector("list", dat_n)       # <- Zi
    dat_J  <- rep(NA, dat_n)              # <- J
    cluster <- unique(dat$id)
    for (i in 1:dat_n) {
      dat_Zi[[i]] <- matrix(c(rep(1, sum(dat$id == cluster[i])),
                              dat[dat$id == cluster[i], c("ses")]),
                            ncol = 2)
      dat_J[i] <- sum(dat$id == cluster[i])
    }
    dat_subjects <- dat$id
    
    # Run MCMC for condition "outps"
    quiet(
      try(
        HW_loop_COND[[outps]] <- MCMC_cross_HWprior(yvec    = dat_yvec, 
                                                    Xmat    = dat_Xmat, 
                                                    Zi      = dat_Zi, 
                                                    J       = dat_J, 
                                                    n       = dat_n, 
                                                    cluster = dat_subjects,
                                                    samsize = MCMC_reps, 
                                                    burnin  = MCMC_burnin), 
        silent = TRUE
      )
    )
  }
    str(HW_loop_COND)

  # Save output
  output <- list(out_lme4 = lme4_loop_COND,
                 out_IW   = IW_loop_COND,
                 out_MF   = MF_loop_COND,
                 out_HW   = HW_loop_COND)
  saveRDS(output, paste0("./output/", "hsb", Sys.Date(), "-ncond_", nconds, "-rep_", MCMC_reps, ".rds"))
  
  
  
  # Define Data for function to use the MCMC Cross function
    # Define Data for function
    yvec <- dat$mathach
    Xmat <- matrix(c(rep(1, length(dat$mathach)),
                         dat$ses,
                         dat$sector,
                         dat$meanses,
                         dat$ses*dat$sector,
                         dat$ses*dat$meanses),
                       ncol = 6)
    n  <- length(unique(dat$id)) # <- n
    Zi <- vector("list", n)       # <- Zi
    J  <- rep(NA, n)              # <- J
    cluster <- unique(dat$id)
    for (i in 1:n) {
      Zi[[i]] <- matrix(c(rep(1, sum(dat$id == cluster[i])),
                              dat[dat$id == cluster[i], c("ses")]),
                            ncol = 2)
      J[i] <- sum(dat$id == cluster[i])
    }
    cluster <- dat$id