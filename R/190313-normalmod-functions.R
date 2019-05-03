### Project:     Master Thesis
### Object:      Posterior Draws functions
### Description: Contains the functions used to draw from the full conditional posterior disitbrutions
###              of the unkown parameters. The full conditonal disitbrutions can be found in the model
###              pdf file in text/ (w/ derivations).
### Date:        2019-03-13

# Normal model:
# y_{ij} = normal(p_{ij})
# y_{ij} = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)

# Packages
  library(lme4)
  library(mvtnorm)
  library(MCMCpack)

# Draws from conditional posterior
draw_theta = function(yvec,Xmat,Zi,bMat,sigma2){
  
  ytilde    <- yvec - c(Zi%*%t(bMat))
  # thetaCovm <- solve(solve(Sigma0) + t(Xmat)%*%Xmat/sigma2)
  # thetaMean <- thetaCovm %*% ( solve(Sigma0)%*%theta0 + t(Xmat)%*%ytilde/sigma2)
  thetaCovm <- solve(t(Xmat)%*%Xmat/sigma2)
  thetaMean <- thetaCovm %*% (t(Xmat)%*%ytilde/sigma2)
  thetaDraw <- rmvnorm(1, thetaMean, thetaCovm)
  return(t(thetaDraw))
  
}

i <- 3
((i-1)*J+1):(i*J)

draw_bMat = function(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J){
  
  ytilde = yvec - c(Xmat%*%theta)
  for(i in 1:n){
    ytilde_i <- ytilde[((i-1)*J+1):(i*J)]
    bi_Covm  <- solve(PsiInv + t(Zi)%*%Zi/sigma2)
    bi_mean  <- bi_Covm %*% (PsiInv %*% bi0 + t(Zi)%*%ytilde_i/sigma2)
    bMat[i,] <- c(rmvnorm(1, bi_mean, bi_Covm))
  }
  return(bMat)
  
}

# sigma^2 (within group variance)
  # sigma^(-2) prior
draw_sigam2_IMPprior = function(yvec,Xmat,Zi,theta,bMat,n,J){
  
  SSR <- t(yvec - Xmat %*% theta - c(Zi %*% t(bMat))) %*% (yvec - Xmat %*% theta - c(Zi %*% t(bMat)))
  sigma2Draw <- rinvgamma(1, (J*n)/2, (SSR)/2)
  return(sigma2Draw)
  
}

  # IG(nu0/2, nu0*S20/2) prior
draw_sigam2_IGprior = function(nu0,sigma20){
  
  SSR <- t(yvec - Xmat %*% theta - c(Zi %*% t(bMat))) %*% (yvec - Xmat %*% theta - c(Zi %*% t(bMat)))
  sigma2Draw <- rinvgamma(1, (nu0 + J*n)/2, (nu0*sigma20 + SSR)/2)
  return(sigma2Draw)
  
}

# Inv-Psi Matrix draws
# Draws w/ Mat-F prior list(nu=2,d=2,e=0,S0=1e3*diag(2))
draw_PsiInv_matF = function(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv,n,d=1,nu=2){
  # Implemented by mixing the scale matrix of an inverse Wishart distribution with a Wishart distribution (eq 3 mulder pericchi)
  # Prior Parameters
    k  <- ncol(bMat)
  # Psi|Omega,.s
      ScaleMatrix = t(bMat)%*%bMat + Omega + 1e-6*diag(2) # current data estimation of random effect covariance matrix + Omega = Previuos draw (or initial guess)
    PsiInvDraw = rwish(v = (n) + (d + k - 1),             # (n) + (d + k - 1) = (posterior) + (prior contribution) [original: n + 2]
                       S = solve(ScaleMatrix)) + 1e-6*diag(2)
    # It's an inverse draw because I would like to get an draw from IW(v, S),
    # but I get a draw from W(v, solve(S)), which becomes a draw from the IW(v, S), when I take its inverse
    # IW(v, S) <=> solve( draw from W(v, solve(S)) )
  # Omega|Psi,.
      ScaleOmega = solve(PsiInvDraw + B0Inv)
    Omega = rwish(v = nu + d + k - 1,           # nu + d + k - 1 = 4 when? (all prior: no information directly available on Omega) [original: 4]
                  S = ScaleOmega)
  return(list(PsiInvDraw,Omega))
}

# Draws w/ HW prior
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

# Draws w/ prior IW(nu0, S0)
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

# # Test Draw Functions ----------------------------------------------------------
# # Test with HW inv-Wish prior
# # Load some data
#   RiesbyDat <- read.table("./data/RiesbyDat.txt")
# 
#   # Make it digestable for the functions
#   yvec     <- RiesbyDat$depr
#   Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
#   J        <- length(unique(RiesbyDat$week))
#   n        <- nrow(RiesbyDat)/J
#   Zi       <- cbind(rep(1,length = J), unique(Xmat[, 2]))
# 
#   subjects <- RiesbyDat$id
# 
# # Define priors
#   bi0      <- rep(0, ncol(Zi))
#   theta0   <- rep(0, ncol(Xmat)) #c(22, -2, 2, -.2) #rep(0, ncol(Xmat))
#   Sigma0   <- t(Xmat) %*% Xmat/(n*3)       # Think about what should this guess be
# 
# # Define initial Values
#   fit <- lmer(depr ~ week + endog + inter + (1 + week | id), data = RiesbyDat, REML = FALSE)
# 
#   theta  <- c(10, 10, 10, 10) #fixef(fit)
#   Psi    <- matrix(c(5,0,0,5), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
#   PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
#   sigma2 <- 1 #attr(VarCorr(fit), "sc")
#   bMat   <- as.matrix(ranef(fit)$id)
#   avec   <- rep(100,2)
# 
# # Run MCMC
#   samsize <- 2000
#   burnin  <- 1/10
#   PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
#   PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
#   PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
#   PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr
#   PD_avec   = matrix(0, nrow = samsize, ncol = 2)
# 
#   #burn-in
#   for(ss in 1:(samsize*burnin)){
#     print(paste("Burn-in:", ss))
# 
#     theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
#     bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv)
#     sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat)
#     outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
#       PsiInv <- outDummy[[1]]
#       avec   <- outDummy[[2]]
# 
#     PD_theta[ss,] = theta
#     PD_bMat[ss,]  = c(bMat)
#     PD_avec[ss,]  = avec
#     PD_Psi[ss,]   = c(solve(PsiInv))
#   }
#   #sampling
#   for(ss in 1:samsize){
#     print(paste("Post-Sample:", ss))
# 
#     theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
#     bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
#     sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
#     outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat)
#       PsiInv <- outDummy[[1]]
#       avec   <- outDummy[[2]]
#     PsiInv <- draw_PsiInv_InvWish(n,bMat,S0=diag(2), e = 1)
# 
# 
#     PD_theta[ss,] = theta
#     PD_bMat[ss,]  = c(bMat)
#     PD_avec[ss,]  = avec
#     PD_Psi[ss,]   = c(solve(PsiInv))
#       sdRI    <- sqrt(PD_Psi[ss,1])       # save the satandard deviation version
#       sdRS    <- sqrt(PD_Psi[ss,4])
#       corRIRS <- PD_Psi[ss,2]/(sdRI*sdRS) # save the correlation
#     PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
#   }
#   out1 <- list(PD_theta, #1
#                PD_bMat,  #2
#                PD_Psi,   #3
#                PD_Psi_sd,#4
#                PD_avec)  #5
# # Results Exploration
#   # Output selection
#   which_o <- 4 # which Object in list
#   which_c <- 1 # which Column in object
#   # Traceplot
#     plot(1:samsize,out1[[which_o]][,which_c],"l")
#   # Estimates
#     mean(out1[[which_o]][,which_c])
#     median(out1[[which_o]][,which_c])
#     # Compare with lme results
#     fit
#   # Posteriors
#     hist(out1[[which_o]][,which_c], breaks = 100,
#          xlab = "Intercept Variance")
#     abline(v = mean(out1[[which_o]][,which_c]), col = "blue", lwd = 2)
#     abline(v = median(out1[[which_o]][,which_c]), col = "red", lwd = 2)