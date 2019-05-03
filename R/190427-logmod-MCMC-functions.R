### Project:     Master Thesis
### Object:      MCMC functions for model fitting
### Description: Contains the functions used to fit all the models of interest. The scripts to actually fit the model is in a different file
###              Comments in the beginning could be used to run this script independently.
### Date:        2019-03-17

# logistic regression mixed model
# y_{ij} = Bernoulli(p_{ij})
# logit(p_{ij}) = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)
#
# Use a Student t distribution to approximate the logistic distribution.
# The Student t distribution has scale=s2tilde=pi**2*(nu-2)/(3*nu) and nu=7.3.
# Use a scale mixture of normals to mimic the Student t distribution.
# This was suggested by Kinney & Dunson (2006).

require(lme4)
require(mvtnorm)
require(MCMCpack)
require(psych)
require(graphics)
require(lattice)
require(logspline)
require(gsl)
require(parallel)
require(doParallel)
require(foreach)

# Functions: Draw from Conditional Posteriors -----------------------------

#posterior draws from conditional posteriors
draw_wvec = function(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec){
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)        # this choice of nu1 and s2tilde allow for nearly exact approximation
                                         # of the logistic model with mixture of normals
  muvec = c(Xmat%*%beta+c(Zi%*%t(bMat))) # values of the log(p(yij)=1) variable
  sigma2vec = s2tilde/phivec             # scale parameter for error variance in the parameter expanded approximated model
  welke1 = which(yvec==1)                # Approximation of binary y_ij as t distributed w_ij
    num1 = length(welke1)
  welke0 = which(yvec==0)
    num0 = length(welke0)
  bounds = pnorm(0,mean=muvec,sd=sqrt(sigma2vec))
  wvecDraw = yvec
  wvecDraw[welke0] = qnorm(runif(num0)*(bounds[welke0]-.00001)+.00001,mean=muvec[welke0],sd=sqrt(sigma2vec[welke0])) #??
  wvecDraw[welke1] = qnorm(runif(num1)*(.99999-bounds[welke1])+bounds[welke1],mean=muvec[welke1],sd=sqrt(sigma2vec[welke1]))
  return(wvecDraw)
}

draw_beta_A = function(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec){
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  sigma2tilde = s2tilde/phivec
  ytilde = wvec - c(Zi%*%t(bMat))
  tXSigma = t(Xmat/sigma2tilde)   # column-wise element-wise division by vector of sigma2tilde
  betaCovm = solve(tXSigma%*%Xmat)
  betaMean = c(betaCovm%*%tXSigma%*%ytilde)
  betaDraw = c(rmvnorm(1,mean=betaMean,sigma=betaCovm))
  return(betaDraw)
}

draw_bMat_A = function(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec){
  ytilde = wvec - c(Xmat%*%beta)
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  for(ii in 1:n){
    phi_i = phivec[((ii-1)*J+1):(ii*J)]
    ytilde_i = ytilde[((ii-1)*J+1):(ii*J)]
    sigma2tilde_i = s2tilde/phi_i
    tZS = t(Zi/sigma2tilde_i)
    tZSZi = tZS%*%Zi
    Sigma_bi = solve(PsiInv+tZSZi)
    mean_bi = Sigma_bi%*%tZS%*%ytilde_i
    bMat[ii,] = c(rmvnorm(1,mean=mean_bi,sigma=Sigma_bi))
  }
  return(bMat)
}

# Draws w/ prior IW(nu0, S0)
draw_PsiInv_InvWish = function(n,bMat,S0=diag(2),e=1){
  # Prior definition
    k <- ncol(bMat)
    nu0 <- k - 1 + e
  # Posterior samples
    ScaleMatrix = t(bMat)%*%bMat + S0
    PsiInvDraw = rwish(v = n + nu0,
                       S = solve(ScaleMatrix))
  return(PsiInvDraw)
}

# General Matrix-F prior (can use to achive, proper neighbour, R*)
draw_PsiInv_matF = function(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv,n,d=1,nu=2){
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

# Posterior of phi (Φ_ij)
  # the vector of precisions in the logistic model approximated as
  # a scale mixture of normals (see eq 8 KinneyDunson2007)
  # The full conditional posterior is reported by KinneyDunson2007 in the supplemental material (checked!)
draw_phivec = function(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec){
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  muvec = c(Xmat%*%beta) + c(Zi%*%t(bMat))
  rate1 = nu1/2+(wvec-muvec)**2/(2*s2tilde)
  shape1 = (nu1+1)/2
  phivecDraw = rgamma(n*J,shape=shape1,rate=rate1)
#  print(c("phivec[1:3]",round(phivecDraw[1:3],3)))
  return(phivecDraw)
}

# Importance weihgts
  # needed because of the apporximation (see KinneyDunson2007 for formula)
compute_weight = function(Xmat,Zi,wvec,beta,bMat){ #?
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  muvec = c(Xmat%*%beta+c(Zi%*%t(bMat)))
  
  weight = exp(sum(dlogis(wvec,location=muvec,log=T) -
                   dt((wvec-muvec)/sqrt(s2tilde),df=nu1,log=T) + .5*log(s2tilde)))
  return(weight)
}

getBounds = function(draws,weights,perc=95){
  probL = (100-perc)/200
  probU = perc/100+probL
  lengte = length(draws)
  draws = draws[order(draws)]
  
  Ftilde = unlist(mclapply(1:length(draws),function(xx){
    mean(weights*(draws<draws[xx]))
  }))/mean(weights)
  
  welkeU = min(which(Ftilde>probU))
  boundU = mean(draws[(welkeU-1):welkeU])
  
  welkeL = max(which(Ftilde<probL))
  boundL = mean(draws[welkeL:(welkeL+1)])
  return(c(boundL,boundU))
}

# Functions: Perform MCMCs ------------------------------------------------

# Help objects
# outps <- 1 # to select the conditon in the other file
# yvec <- dat_yvec     <- dat$yvec
# Xmat <- dat_Xmat     <- cbind(rep(1, nrow(dat)), dat$xvec, dat$cvec, dat$inter)
# J    <- dat_J        <- length(unique(dat$xvec))
# n    <- dat_n        <- nrow(dat)/dat_J
# Zi   <- dat_Zi       <- cbind(rep(1,length = dat_J), unique(dat_Xmat[, 2]))

MCMC_mixlogreg_matF = function(yvec, Xmat, Zi, J, n, iniv = 1, B0 = list(nu=2,d=1,e=0,S0=1e3*diag(2)), samsize = 1e4, burnin = 1/10){
  # Data Prep
  xvec    <- Xmat[, 2]
  covar   <- Xmat[, 3]
  inter   <- Xmat[, 4]
  cluster <- rep(seq(1, n), each = J)
  # Define priors
  bi0   <- rep(0, ncol(Zi)) # prior for bi
  S0    <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    glm2 <- glm(yvec ~ xvec + covar + inter, family=binomial)
    betaHat0 = c(coefficients(glm2))
    pivec = exp(Xmat%*%betaHat0)/(exp(Xmat%*%betaHat0)+1)
    weightVec = pivec*(1-pivec)
    Rstar = solve(t(Zi*weightVec[1:J])%*%Zi)
    S0        <- 2*Rstar # should this be k * Rstar?
    S0Inv     <- solve(S0)
  }
  
  # Define initial Values
  if (iniv == 1){
    fit   <- glmer(yvec ~ xvec + covar + inter + (1 + xvec|cluster), family=binomial)
    #theta  <- fixef(fit)
    Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)# as in MulderPericchi2018 code
    PsiInv <- solve(Psi + 1e-6*diag(2))        
    sigma2 <- attr(VarCorr(fit), "sc")
    bMat   <- as.matrix(ranef(fit)$cluster)
    beta   <- fixef(fit)
  } else {
    #theta  <- t(rep(0, ncol(Xmat)))
    Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
    PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
    sigma2 <- 1 #attr(VarCorr(fit), "sc")
    bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
  }
  Omega = S0
  
  #Student t approximation of logistic regression distribution
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1) # used inside the sampling functions, not defined as input
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2)

  #store draws
  MbetaA   = matrix(0,nrow=samsize,ncol=4)
  MPsiA    = matrix(0,nrow=samsize,ncol=4)
  PD_Psi_sd= matrix(0,nrow=samsize,ncol=4)
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  MOmegaA  = matrix(0,nrow=samsize,ncol=4)
  MbMatA   = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA   = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:(burnin*samsize)){
    print(paste("Burn-in:", ss))
    
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    outDummy = draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega = Omega, B0Inv = S0, n = n, d=B0$d, nu=B0$nu)
      PsiInv = outDummy[[1]]
      Omega = outDummy[[2]]
      
  }
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    outDummy = draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega = Omega, B0Inv = S0, n = n, d=B0$d, nu=B0$nu)
      PsiInv = outDummy[[1]]
      Omega = outDummy[[2]]
    
    #compute weights
    MimpwA[ss,1]  = compute_weight(Xmat,Zi,wvec,beta,bMat)
    MbetaA[ss,]   = beta
    MPsiA[ss,]    = c(solve(PsiInv))
      sdRI    <- sqrt(MPsiA[ss,1])       # save the satandard deviation version
      sdRS    <- sqrt(MPsiA[ss,4])
      corRIRS <- MPsiA[ss,2]/(sdRI*sdRS) # save the correlation
    PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
    MOmegaA[ss,]  = Omega
    MbMatA[ss,]   = c(bMat)
    MphivecA[ss,] = phivec
  }
  
  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean    = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean    = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean     = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds  = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws = MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  
  B0$S0 <- S0 # save the prior again because in R* case this is different at the end
  out <- list(PD_theta   = MbetaA,
              PD_bMat    = MbMatA,
              PD_Psi     = MPsiA,
              PD_Psi_sd  = PD_Psi_sd,
              PD_Omg     = MOmegaA,
              hyperpar   = B0,
              MphivecA   = MphivecA,
              betaMean   = betaMean,
              bMatMean   = bMatMean,
              PsiMean    = PsiMean,
              betaBounds = betaBounds,
              bMatBounds = bMatBounds,
              PsiBounds  = PsiBounds,
              MbetaA     = MbetaA,
              MimpwA     = MimpwA)
  return(out)
}
  
MCMC_mixlogreg_InvWish = function(yvec, Xmat, Zi, J, n, iniv = 1, B0 = list(e=1,S0=diag(2)), samsize = 2000, burnin = 1/10){
 # Data Prep
  xvec    <- Xmat[, 2]
  covar   <- Xmat[, 3]
  inter   <- Xmat[, 4]
  cluster <- rep(seq(1, n), each = J)
  # Define priors
  bi0   <- rep(0, ncol(Zi)) # prior for bi
  e     <- B0$e # small quantity to add to prior nu
  S0    <- B0$S0
  if(is.matrix(S0) == TRUE){
    S0Inv <- solve(S0)
  } else {
    glm2 <- glm(yvec ~ xvec + covar + inter, family=binomial)
    betaHat0 = c(coefficients(glm2))
    pivec = exp(Xmat%*%betaHat0)/(exp(Xmat%*%betaHat0)+1)
    weightVec = pivec*(1-pivec)
    Rstar = solve(t(Zi*weightVec[1:J])%*%Zi)
    S0        <- 2*Rstar # should this be k * Rstar?
    S0Inv     <- solve(S0)
  }
  
  # Define initial Values
  if (iniv == 1){
    fit    <- glmer(yvec ~ xvec + covar + inter + (1 + xvec|cluster), family=binomial)
    Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)# as in MulderPericchi2018 code
    PsiInv <- solve(Psi + 1e-6*diag(2))        
    sigma2 <- attr(VarCorr(fit), "sc")
    bMat   <- as.matrix(ranef(fit)$cluster)
    beta   <- fixef(fit)
  } else {
    #theta  <- t(rep(0, ncol(Xmat)))
    Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
    PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
    sigma2 <- 1 #attr(VarCorr(fit), "sc")
    bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
  }

  #Student t approximation of logistic regression distribution
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1) # used inside the sampling functions, not defined as input
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2)

  #store draws
  MbetaA   = matrix(0,nrow=samsize,ncol=4)
  MPsiA    = matrix(0,nrow=samsize,ncol=4)
  PD_Psi_sd= matrix(0,nrow=samsize,ncol=4)
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  MbMatA   = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA   = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:(burnin*samsize)){
    print(paste("Burn-in:", ss))
    
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    PsiInv = draw_PsiInv_InvWish(n, bMat, S0 = S0, e = e)
  }
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    PsiInv = draw_PsiInv_InvWish(n, bMat, S0 = S0, e = e)
    
    #compute weights
    MimpwA[ss,1]  = compute_weight(Xmat,Zi,wvec,beta,bMat)
    MbetaA[ss,]   = beta
    MPsiA[ss,]    = c(solve(PsiInv))
      sdRI    <- sqrt(MPsiA[ss,1])       # save the satandard deviation version
      sdRS    <- sqrt(MPsiA[ss,4])
      corRIRS <- MPsiA[ss,2]/(sdRI*sdRS) # save the correlation
    PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
    MbMatA[ss,]   = c(bMat)
    MphivecA[ss,] = phivec
  }
  
  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean    = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean    = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean     = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds  = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws = MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  
  B0$S0 <- S0 # save the prior again because in R* case this is different at the end
  out <- list(PD_theta   = MbetaA,
              PD_bMat    = MbMatA,
              PD_Psi     = MPsiA,
              PD_Psi_sd  = PD_Psi_sd,
              hyperpar   = B0,
              MphivecA   = MphivecA,
              betaMean   = betaMean,
              bMatMean   = bMatMean,
              PsiMean    = PsiMean,
              betaBounds = betaBounds,
              bMatBounds = bMatBounds,
              PsiBounds  = PsiBounds,
              MbetaA     = MbetaA,
              MimpwA     = MimpwA)
  return(out)
}

MCMC_mixlogreg_HW = function(yvec, Xmat, Zi, J, n, iniv = 1, samsize = 2000, burnin = 1/10){
  
  # Data Prep
  xvec    <- Xmat[, 2]
  covar   <- Xmat[, 3]
  inter   <- Xmat[, 4]
  cluster <- rep(seq(1, n), each = J)
  # Define priors
  bi0   <- rep(0, ncol(Zi)) # prior for bi
  # Define initial Values
  if (iniv == 1){
    fit   <- glmer(yvec ~ xvec + covar + inter + (1 + xvec|cluster), family=binomial)
    Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)# as in MulderPericchi2018 code
    PsiInv <- solve(Psi + 1e-6*diag(2))        
    sigma2 <- attr(VarCorr(fit), "sc")
    bMat   <- as.matrix(ranef(fit)$cluster)
    beta   <- fixef(fit)
  } else {
    Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
    PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
    sigma2 <- 1 #attr(VarCorr(fit), "sc")
    bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
  }
  avec   <- rep(100,2)

  #Student t approximation of logistic regression distribution
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1) # used inside the sampling functions, not defined as input
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2)

  #store draws
  MbetaA   = matrix(0,nrow=samsize,ncol=4)
  MPsiA    = matrix(0,nrow=samsize,ncol=4)
  PD_Psi_sd= matrix(0,nrow=samsize,ncol=4)
    colnames(PD_Psi_sd) <- c("RI_sd", "cor", "cor", "RS_sd")
  MbMatA   = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA   = matrix(0,nrow=samsize,ncol=1)
  PD_avec   = matrix(0, nrow = samsize, ncol = 2)
  
  #burn-in
  for(ss in 1:(burnin*samsize)){
    print(paste("Burn-in:", ss))
    
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
      PsiInv <- outDummy[[1]]
      avec   <- outDummy[[2]]
      
  }
  
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,n,J,wvec,beta,bMat,PsiInv,phivec)
    outDummy <- draw_PsiInv_HW(PsiInv,avec,bMat,n)
      PsiInv <- outDummy[[1]]
      avec   <- outDummy[[2]]
    
    #compute weights
    MimpwA[ss,1]  = compute_weight(Xmat,Zi,wvec,beta,bMat)
    #store results
    MbetaA[ss,]   = beta
    MPsiA[ss,]    = c(solve(PsiInv))
      sdRI    <- sqrt(MPsiA[ss,1])       # save the satandard deviation version
      sdRS    <- sqrt(MPsiA[ss,4])
      corRIRS <- MPsiA[ss,2]/(sdRI*sdRS) # save the correlation
    PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
    PD_avec[ss,]   = avec
    MbMatA[ss,]    = c(bMat)
    MphivecA[ss,]  = phivec
    
  }
  
  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean    = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean    = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean     = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds  = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws = MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  
  out <- list(PD_theta   = MbetaA,
              PD_bMat    = MbMatA,
              PD_Psi     = MPsiA,
              PD_Psi_sd  = PD_Psi_sd,
              MphivecA   = MphivecA,
              betaMean   = betaMean,
              bMatMean   = bMatMean,
              PsiMean    = PsiMean,
              betaBounds = betaBounds,
              bMatBounds = bMatBounds,
              PsiBounds  = PsiBounds,
              MbetaA     = MbetaA,
              MimpwA     = MimpwA)
  return(out)
}


# # Check results (posterior and traceplots)
#   out$PD_Psi   
#   par <- 1
#   # Posterior
#   plot(density(out$PD_Psi_sd[,par]))
#   median(out$PD_Psi_sd[,par])
#   hist(out$PD_Psi_sd[,par],
#        xlim = c(0, 25))
#   # Traceplot
#   x <- out$PD_Psi_sd[,par]
#   plot(1:samsize, x,"l",
#        xlim = c(0, samsize))
#   # Point estimates
#   glme4_loop_COND$`n_400:J_7`$glme_Psi_sd
#   colMeans(out$PD_Psi_sd)