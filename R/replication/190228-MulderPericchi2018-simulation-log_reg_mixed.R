
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

detectCores()
getDoParWorkers()

# Functions: Generate Data ------------------------------------------------

#generate data
generate_yvec_logreg = function(n = 30,                          # number of clusters, ie individuals
                                beta = c(-.625,.25,-.25,.125),   # true fixed effects (intercept, time, covariate, interaction)
                                PsiMat = matrix(c(1,0,0,0),      # var-cov matrix of random effects
                                                ncol = 2)){
  J = 7                                                          # cluster size (nj), ie n of repeated measures
  xvec = c(rep(1, length = J*n/2), rep(0, length = J*n/2))       # Covariate (eg treatment)
  tvec = rep(-3:3, length = n*J)                                 # Time varibale
  Xmat = matrix(c(rep(1,length=J*n),tvec,xvec,tvec*xvec),ncol=4) # Design matrix
  Zi = matrix(c(rep(1,length=J),-3:3),ncol=2)                    # Predictors with random effects (intercept and random slope for time)
  bMat = rmvnorm(n, mean = c(0, 0), sigma = PsiMat)              # Random effects b, one per predictor w/ random eff, per individual
  muvec = c(Xmat %*% beta) + c(Zi %*% t(bMat))                   # Outcome variable (in log form) for each measurement (n*J)
                                                                 # is the sum of two components:
                                                                 # - the fixed effect contribution c(Xmat %*% beta)
                                                                 #   (a vector of contributions constant between clusters/individuals)
                                                                 # - the random effect contribution c(Zi %*% t(bMat))
                                                                 #   (a vector of contributions of constant within a cluster/individual)
  probvec = exp(muvec)/(exp(muvec)+1)                            # Transform: log -> probability
  yvec = 1*(runif(n*J)<probvec)                                  # Transform: ?
  return(list(yvec,Xmat,bMat))                                   # Return Dataset [outcome], [design matrix], [random effects]
}

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
#  print(c("wvec",wvecDraw[1:5]))
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
#  print(c("beta",round(betaDraw,3)))
  return(betaDraw)
}

draw_bMat_A = function(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec){
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

#  print(c("bMati",round(bMat[1,],3)))
  return(bMat)
}

# Using prior invW(nu0, solve(S0)), where S0 is a prior guess for the random effects var-cov matrix
draw_PsiInv_InvWish = function(n,bMat,B0=diag(2)){
  ScaleMatrix = t(bMat)%*%bMat + B0
  PsiInvDraw = rwish(v = n + 2, # nu = 2
                     S = solve(ScaleMatrix))
#  print(c("Psi[1,]",round(PsiInvDraw[1,],3)))
  return(PsiInvDraw)
}

draw_PsiInv_JB = function(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec){
  ScaleMatrix = t(bMat)%*%bMat
  PsiInvDraw = rwish(v = n - 2,
                     S = solve(ScaleMatrix))
#  print(c("Psi[1,]",round(PsiInvDraw[1,],3)))
  return(PsiInvDraw)
}

# Proper neighbour of Sigma**-(1/2) ?
draw_PsiInv_MSBeta2 = function(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec,Omega,B0Inv){
  ScaleMatrix = t(bMat)%*%bMat + Omega # current data estimation of random effect covariance matrix
                                       # + Omega = Previuos draw (or initial guess)
  PsiInvDraw = rwish(v = n + 2,
                     S = solve(ScaleMatrix))
  
  ScaleOmega = solve(PsiInvDraw + B0Inv)
  Omega = rwish(v = 4,                 # k = 4 (4 predicotrs: intercept, time, covariate, time*covariate)
                S = ScaleOmega)
  
  return(list(PsiInvDraw,Omega))
}

# Posterior Using HW2013 inv-Wishart prior
  # Look at section 4 of the paper for the full conditional of Sigma, and a_k
  # How is it = F(nu = 1, δ = 2, B = 10**5)?
draw_PsiInv_HW = function(PsiInv,avec,bMat){
  #prior of Huang & Wand with scale A_k=10^5 and nu=2.
  Ak = 10**5
  ScaleMatrix = t(bMat)%*%bMat +               # sum of ui %*% ui'
                4*diag(1/avec)                 # nu = 2
  PsiInvDraw = rwish(v = n + 3,                # nu + N + q - 1, where: nu = 2, q = 2, n = 30 (number of clusters), 
                                               # result is m + 3
                     S = solve(ScaleMatrix))
  
  scale_avec = 2*diag(PsiInvDraw) + 1/Ak**2    # a_k scale parameter in Haung Wnad 2013 (section 4 full conditionals)
  avec = rinvgamma(2,                          # q = 2: 1 random intercept, 1 random slope
                   shape = 2,                  # (nu + q) /2, where q = 2, nu = 2
                   scale = scale_avec)
  
  return(list(PsiInvDraw,avec))
}

# Posterior of phi (Φ_ij)
  # the vector of precisions in the logistic model approximated as
  # a scale mixture of normals (see eq 8 KinneyDunson2007)
  # The full conditional posterior is reported by KinneyDunson2007 in the supplemental material (checked!)
draw_phivec = function(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec){
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

MCMC_mixlogreg_JB_A = function(yvec,n,J,samsize=1e3){
  xvec = c(rep(1,length=J*n/2),rep(0,length=J*n/2))
  tvec = rep(-3:3,length=n*J)
  Xmat = matrix(c(rep(1,length=J*n),tvec,xvec,tvec*xvec),ncol=4)
  Zi = matrix(c(rep(1,length=J),-3:3),ncol=2)
  subject = rep(1:n,each=J)#nvec
  
  #Student t approximation of logistic regression distribution
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  #scale1

  #estimate model with lme
  data1 = data.frame(y=yvec,x=xvec,t=tvec,inact=tvec*xvec,subject=subject)
  glm1 = glmer(y ~ 1 + t + x + inact + (1+t|subject), family=binomial, data=data1)
  
  #glm1 = glmer(y ~ 1 + t + x + inact + (1|subject), family=binomial, data=data1)

  #get estimates
  PsiHat = matrix(VarCorr(glm1)[[1]][1:4],ncol=2)
  betaHat = fixef(glm1)
  bHat = ranef(glm1)[[1]]

  #derive Kass & Natarajan (2006) hyperparameter for prior of Psi
  glm2 = glm(y ~ 1 + t + x + inact, family=binomial, data=data1)
  betaHat0 = c(coefficients(glm2))
  pivec = exp(Xmat%*%betaHat0)/(exp(Xmat%*%betaHat0)+1) # predicted probabilities with fixed effects only
  weightVec = pivec*(1-pivec)                           # <- the diagonal of the weight matrix in Kass Natarajan!
  #note W_i and Z_i are equal for all i in this setup (no random effects)
  Rstar = solve(t(Zi*weightVec[1:J])%*%Zi)
#  Rmat = matrix(0,ncol=2,nrow=2)
#  Rmat = 0
#  for(rr in 1:n){
#    Rmat = Rmat + t(Zi*weightVec[((rr-1)*J+1):(rr*J)])%*%Zi
#    #Rmat = Rmat + t(rep(1,J)*weightVec[((rr-1)*J+1):(rr*J)])%*%rep(1,J)
#  }
#  Rstar = solve(Rmat/n)
  
  #initial values
  beta = betaHat
  Psi = PsiHat + .01*diag(2) # ? why adding .01 to the vairances?
  PsiInv = solve(Psi)        # var-covar matrix of random effects estimated with glmer
  bMat = as.matrix(bHat)     # fixed effects estimated with glm
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2) # ?

  #store draws
  MbetaA = matrix(0,nrow=samsize,ncol=4)
  MPsiA = matrix(0,nrow=samsize,ncol=4)
  MbMatA = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:1e2){
    Sys.time()
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    Sys.time()
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    Sys.time()
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    Sys.time()
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    Sys.time()
    PsiInv = draw_PsiInv_JB(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    Sys.time()
  }
  
  for(ss in 1:samsize){
    wvec   = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta   = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat   = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    PsiInv = draw_PsiInv_JB(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    #PsiInv = draw_PsiInv_InvWish(yvec,wvec,beta,bMat,PsiInv,phivec,Omega0 = diag(2)*.5)
    
    #compute weights
    MimpwA[ss,1] = compute_weight(Xmat,Zi,wvec,beta,bMat)
    
    MbetaA[ss,] = beta
    MPsiA[ss,] = c(solve(PsiInv))
    MbMatA[ss,] = c(bMat)
    MphivecA[ss,] = phivec
  }

  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws=MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  return(list(betaMean,bMatMean,PsiMean,betaBounds,bMatBounds,PsiBounds,MbetaA,MbMatA,MPsiA,MphivecA,MimpwA,betaHat,bHat,PsiHat))
}

# Simulation using F(nu = k, δ = 1, B = R* or R_hat)

MCMC_mixlogreg_SBeta2_A = function(yvec,n,J,B0=1,samsize=1e3){
  xvec = c(rep(1,length=J*n/2),rep(0,length=J*n/2))
  tvec = rep(-3:3,length=n*J)
  Xmat = matrix(c(rep(1,length=J*n),tvec,xvec,tvec*xvec),ncol=4)
  Zi = matrix(c(rep(1,length=J),-3:3),ncol=2)
  subject = rep(1:n,each=J)#nvec
  
  #Student t approximation of logistic regression distribution
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  #scale1
  
  #estimate model with lme
  data1 = data.frame(y=yvec,x=xvec,t=tvec,inact=tvec*xvec,subject=subject)
  glm1 = glmer(y ~ 1 + t + x + inact + (1+t|subject), family=binomial, data=data1)
  
  #glm1 = glmer(y ~ 1 + t + x + inact + (1|subject), family=binomial, data=data1)
  
  #get estimates
  PsiHat = matrix(VarCorr(glm1)[[1]][1:4],ncol=2)
  betaHat = fixef(glm1)
  bHat = ranef(glm1)[[1]]
  #else the input B0 is the prior scale matrix
  
  #derive Kass & Natarajan (2006) hyperparameter for prior of Psi
  glm2 = glm(y ~ 1 + t + x + inact, family=binomial, data=data1)
  betaHat0 = c(coefficients(glm2))
  pivec = exp(Xmat%*%betaHat0)/(exp(Xmat%*%betaHat0)+1)
  weightVec = pivec*(1-pivec)
  #note W_i and Z_i are equal for all i in this setup
  Rstar = solve(t(Zi*weightVec[1:J])%*%Zi)
  if(is.matrix(B0)==F){ # B0==1, then use default prior scale of Kass&Nat'06
                        # B0==2, then set B0 equal to the ML estimate
    if(B0==1){
      B0 = Rstar
      if(min(eigen(Rstar)$values)<.0001){
        B0 = Rstar+diag(2)*.0001
      }
    }
    else{
      B0 = PsiHat
      if(min(eigen(Rstar)$values)<.0001){
        B0 = PsiHat+diag(2)*.0001
      }
    }
  }
  #  Rmat = matrix(0,ncol=2,nrow=2)
  #  Rmat = 0
  #  for(rr in 1:n){
  #    Rmat = Rmat + t(Zi*weightVec[((rr-1)*J+1):(rr*J)])%*%Zi
  #    #Rmat = Rmat + t(rep(1,J)*weightVec[((rr-1)*J+1):(rr*J)])%*%rep(1,J)
  #  }
  #  Rstar = solve(Rmat/n)
  
  #initial values
  beta = betaHat
  Psi = PsiHat + .01*diag(2)
  PsiInv = solve(Psi)
  bMat = as.matrix(bHat)
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2)
  Omega = B0
  B0Inv = solve(B0)
  
  
  #store draws
  MbetaA = matrix(0,nrow=samsize,ncol=4)
  MPsiA = matrix(0,nrow=samsize,ncol=4)
  MOmegaA = matrix(0,nrow=samsize,ncol=4)
  MbMatA = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:1e2){
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    #PsiInv = draw_PsiInv_JB(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    outDummy = draw_PsiInv_MSBeta2(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec,Omega,B0Inv)
    PsiInv = outDummy[[1]]
    Omega = outDummy[[2]]
  }
  
  for(ss in 1:samsize){
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    #PsiInv = draw_PsiInv_JB(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    outDummy = draw_PsiInv_MSBeta2(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec,Omega,B0Inv)
    PsiInv = outDummy[[1]]
    Omega = outDummy[[2]]
    
    #compute weights
    MimpwA[ss,1] = compute_weight(Xmat,Zi,wvec,beta,bMat)
    
    MbetaA[ss,] = beta
    MPsiA[ss,] = c(solve(PsiInv))
    MOmegaA[ss,] = Omega
    MbMatA[ss,] = c(bMat)
    MphivecA[ss,] = phivec
  }
  
  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws=MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  return(list(betaMean,bMatMean,PsiMean,betaBounds,bMatBounds,PsiBounds,MbetaA,MbMatA,MPsiA,MphivecA,MimpwA,betaHat,bHat,PsiHat,MOmegaA))
}

MCMC_mixlogreg_InvWish_A = function(yvec,n,J,B0=1,samsize=1e3){
  xvec = c(rep(1,length=J*n/2),rep(0,length=J*n/2))
  tvec = rep(-3:3,length=n*J)
  Xmat = matrix(c(rep(1,length=J*n),tvec,xvec,tvec*xvec),ncol=4)
  Zi = matrix(c(rep(1,length=J),-3:3),ncol=2)
  subject = rep(1:n,each=J)#nvec
  
  #Student t approximation of logistic regression distribution
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)
  #scale1
  
  #estimate model with lme
  data1 = data.frame(y=yvec,x=xvec,t=tvec,inact=tvec*xvec,subject=subject)
  glm1 = glmer(y ~ 1 + t + x + inact + (1+t|subject), family=binomial, data=data1)
  
  #glm1 = glmer(y ~ 1 + t + x + inact + (1|subject), family=binomial, data=data1)
  
  #get estimates
  PsiHat = matrix(VarCorr(glm1)[[1]][1:4],ncol=2)
  betaHat = fixef(glm1)
  bHat = ranef(glm1)[[1]]
  #else the input B0 is the prior scale matrix
  
  #derive Kass & Natarajan (2006) hyperparameter for prior of Psi
  glm2 = glm(y ~ 1 + t + x + inact, family=binomial, data=data1)
  betaHat0 = c(coefficients(glm2))
  pivec = exp(Xmat%*%betaHat0)/(exp(Xmat%*%betaHat0)+1)
  weightVec = pivec*(1-pivec)
  #note W_i and Z_i are equal for all i in this setup
  Rstar = solve(t(Zi*weightVec[1:J])%*%Zi)
  #  Rmat = matrix(0,ncol=2,nrow=2)
  #  Rmat = 0
  #  for(rr in 1:n){
  #    Rmat = Rmat + t(Zi*weightVec[((rr-1)*J+1):(rr*J)])%*%Zi
  #    #Rmat = Rmat + t(rep(1,J)*weightVec[((rr-1)*J+1):(rr*J)])%*%rep(1,J)
  #  }
  #  Rstar = solve(Rmat/n)
  if(is.matrix(B0)==F){#B0==1, then use default prior of Kass&Nat'06, B0==2, then set B0 equal to the ML estimate
    if(B0==1){
      B0 = Rstar
      if(min(eigen(Rstar)$values)<.0001){
        B0 = Rstar+diag(2)*.0001
      }
    }
    else{
      B0 = PsiHat
      if(min(eigen(Rstar)$values)<.0001){
        B0 = PsiHat+diag(2)*.0001
      }
    }
  }
  
  #initial values
  beta = betaHat
  Psi = PsiHat + .01*diag(2)
  PsiInv = solve(Psi)
  bMat = as.matrix(bHat)
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2)
  
  #store draws
  MbetaA = matrix(0,nrow=samsize,ncol=4)
  MPsiA = matrix(0,nrow=samsize,ncol=4)
  MbMatA = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:1e2){
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    #PsiInv = draw_PsiInv_JB(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    PsiInv = draw_PsiInv_InvWish(n,bMat,B0=B0)
  }
  
  for(ss in 1:samsize){
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    #PsiInv = draw_PsiInv_JB(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    PsiInv = draw_PsiInv_InvWish(n,bMat,B0=B0)
    
    #compute weights
    MimpwA[ss,1] = compute_weight(Xmat,Zi,wvec,beta,bMat)
    
    MbetaA[ss,] = beta
    MPsiA[ss,] = c(solve(PsiInv))
    MbMatA[ss,] = c(bMat)
    MphivecA[ss,] = phivec
  }
  
  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws=MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  return(list(betaMean,bMatMean,PsiMean,betaBounds,bMatBounds,PsiBounds,MbetaA,MbMatA,MPsiA,MphivecA,MimpwA,betaHat,bHat,PsiHat))
}

# Simulation using F(nu = k, δ = 2, B = 10**10*Ik)

MCMC_mixlogreg_HW_A = function(yvec,n,J,samsize=1e3){
  #prior of Huang and Wand (2013) with uniform priors for the correlations
  #and half-t(nu=2,Ak=10^5) for the standard deviations.
  xvec = c(rep(1,length=J*n/2),rep(0,length=J*n/2))
  tvec = rep(-3:3,length=n*J)
  Xmat = matrix(c(rep(1,length=J*n),tvec,xvec,tvec*xvec),ncol=4)
  Zi = matrix(c(rep(1,length=J),-3:3),ncol=2)
  subject = rep(1:n,each=J)#nvec
  
  #Student t approximation of logistic regression distribution (Kinney & Dunson, 2006)
  nu1 = 7.3
  s2tilde = pi**2*(nu1-2)/(3*nu1)

  #estimate model with lme
  data1 = data.frame(y=yvec,x=xvec,t=tvec,inact=tvec*xvec,subject=subject)
  glm1 = glmer(y ~ 1 + t + x + inact + (1+t|subject), family=binomial, data=data1)
  
  #get estimates
  PsiHat = matrix(VarCorr(glm1)[[1]][1:4],ncol=2)
  betaHat = fixef(glm1)
  bHat = ranef(glm1)[[1]]

  #derive Kass & Natarajan (2006) hyperparameter for prior of Psi
  glm2 = glm(y ~ 1 + t + x + inact, family=binomial, data=data1)
  betaHat0 = c(coefficients(glm2))
  pivec = exp(Xmat%*%betaHat0)/(exp(Xmat%*%betaHat0)+1)
  weightVec = pivec*(1-pivec)
  #note W_i and Z_i are equal for all i in this setup
  Rstar = solve(t(Zi*weightVec[1:J])%*%Zi)
  
  #initial values
  beta = betaHat
  Psi = PsiHat + .01*diag(2)
  PsiInv = solve(Psi)
  bMat = as.matrix(bHat)
  phivec = rgamma(n*J,shape=nu1/2,rate=nu1/2)
  avec = rep(100,2)
  
  #store draws
  MbetaA = matrix(0,nrow=samsize,ncol=4)
  MPsiA = matrix(0,nrow=samsize,ncol=4)
  MavecA = matrix(0,nrow=samsize,ncol=2)
  MbMatA = matrix(0,nrow=samsize,ncol=n*2)
  MphivecA = matrix(0,nrow=samsize,ncol=n*J)
  MimpwA = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:1e2){
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    outDummy = draw_PsiInv_HW(PsiInv,avec,bMat)
    PsiInv = outDummy[[1]]
    avec = outDummy[[2]]
  }
  
  for(ss in 1:samsize){
    wvec = draw_wvec(yvec,Xmat,Zi,beta,bMat,PsiInv,phivec)
    beta = draw_beta_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    phivec = draw_phivec(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    bMat = draw_bMat_A(yvec,Xmat,Zi,wvec,beta,bMat,PsiInv,phivec)
    outDummy = draw_PsiInv_HW(PsiInv,avec,bMat)
    PsiInv = outDummy[[1]]
    avec = outDummy[[2]]
    
    #compute weights
    MimpwA[ss,1] = compute_weight(Xmat,Zi,wvec,beta,bMat)
    
    MbetaA[ss,] = beta
    MPsiA[ss,] = c(solve(PsiInv))
    MbMatA[ss,] = c(bMat)
    MphivecA[ss,] = phivec
    MavecA[ss,] = avec
  }
  
  #compute estimates and credibility bounds
  meanWeights = mean(MimpwA)
  betaMean = c(t(MimpwA)%*%MbetaA/meanWeights)/samsize
  bMatMean = c(t(MimpwA)%*%MbMatA/meanWeights)/samsize
  PsiMean = c(t(MimpwA)%*%MPsiA/meanWeights)/samsize
  betaBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    betaBounds[bb,] = getBounds(draws=MbetaA[,bb],weights=MimpwA[,1],perc=95)
  }
  bMatBounds = matrix(0,ncol=2,nrow=2*n)
  for(bb in 1:(2*n)){
    bMatBounds[bb,] = getBounds(draws=MbMatA[,bb],weights=MimpwA[,1],perc=95)
  }
  PsiBounds = matrix(0,ncol=2,nrow=4)
  for(bb in 1:4){
    PsiBounds[bb,] = getBounds(draws=MPsiA[,bb],weights=MimpwA[,1],perc=95)
  }
  return(list(betaMean,bMatMean,PsiMean,betaBounds,bMatBounds,PsiBounds,MbetaA,MbMatA,MPsiA,MphivecA,MimpwA,betaHat,bHat,PsiHat,MavecA))
}

# Perform Simulations -----------------------------------------------------
#> Basic simulation ####
n = 30
J = 7
datagen = generate_yvec_logreg(n=n,beta = c(-.625,.25,-.25,.125),PsiMat = matrix(c(1,0,0,1),ncol=2))
yvec = datagen[[1]]
Xmat = datagen[[2]]
bMat = datagen[[3]]

nu1=7.3
s2tilde=pi**2*(nu1-2)/(3*nu1)

samsize = 1e4
out1 = MCMC_mixlogreg_JB_A(yvec,n,J,samsize=samsize)                      # uniform prior on Sigma2 (with Jacobian Transform) / or improper?
out2 = MCMC_mixlogreg_SBeta2_A(yvec,n,J,B0 = 1e3*diag(2),samsize=samsize) # proper neighbour
out3 = MCMC_mixlogreg_InvWish_A(yvec,n,J,B0 = diag(2),samsize=samsize)    #
out4 = MCMC_mixlogreg_SBeta2_A(yvec,n,J,B0 = 1,samsize=samsize)           # mat-F w/ R*
out5 = MCMC_mixlogreg_InvWish_A(yvec,n,J,B0 = 1,samsize=samsize)          #
out6 = MCMC_mixlogreg_SBeta2_A(yvec,n,J,B0 = 2,samsize=samsize)           # mat-F w/ Psi_Hat (ML)
out7 = MCMC_mixlogreg_HW_A(yvec,n,J,samsize=samsize)                      # inv-Wish à la Huang Wand

out1[[12]]
out1[[1]]
out2[[1]]
out3[[1]]
out4[[1]]
out5[[1]]
out6[[1]]
out7[[2]]

#check ene beta
ene = 1
plot(1:1e4,out1[[7]][,ene],"l")
lines(1:1e4,out2[[7]][,ene],col=2)
lines(1:1e4,out3[[7]][,ene],col=3)
lines(1:1e4,out4[[7]][,ene],col=4)
lines(1:1e4,out5[[7]][,ene],col=5)
lines(1:1e4,out6[[7]][,ene],col=6)
lines(1:1e4,out7[[7]][,ene],col=7)

ymax = 1.2
plot(density(out1[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax))
par(new=T)
plot(density(out2[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax),col=2)
par(new=T)
plot(density(out3[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax),col=3)
par(new=T)
plot(density(out4[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax),col=4)
par(new=T)
plot(density(out5[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax),col=5)
par(new=T)
plot(density(out6[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax),col=6)
par(new=T)
plot(density(out7[[7]][,ene]),xlim=c(-4,1),ylim=c(0,ymax),col=7)

#check ene random effect b_i
ene1 = 1
plot(1:1e4,out1[[8]][,ene],"l")
lines(1:1e4,out2[[8]][,ene],col=2)
lines(1:1e4,out3[[8]][,ene],col=3)
lines(1:1e4,out4[[8]][,ene],col=4)
lines(1:1e4,out5[[8]][,ene],col=5)
lines(1:1e4,out6[[8]][,ene],col=6)
lines(1:1e4,out7[[8]][,ene],col=7)

plot(density(out1[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3))
par(new=T)
plot(density(out2[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3),col=2)
par(new=T)
plot(density(out3[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3),col=3)
par(new=T)
plot(density(out4[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3),col=4)
par(new=T)
plot(density(out5[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3),col=5)
par(new=T)
plot(density(out6[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3),col=6)
par(new=T)
plot(density(out7[[8]][,ene]),xlim=c(-3,2),ylim=c(0,1.3),col=7)

#check ene Psi
ene = 2
plot(1:1e4,out4[[9]][,ene],"l")
lines(1:1e4,out2[[9]][,ene],col=2)
lines(1:1e4,out3[[9]][,ene],col=3)
lines(1:1e4,out4[[9]][,ene],col=4)
lines(1:1e4,out5[[9]][,ene],col=5)
lines(1:1e4,out6[[9]][,ene],col=6)
lines(1:1e4,out7[[9]][,ene],col=7)

plot(density(out1[[9]][,ene]),xlim=c(0,8),ylim=c(0,2))
par(new=T)
plot(density(out2[[9]][,ene]),xlim=c(0,8),ylim=c(0,2),col=2)
par(new=T)
plot(density(out3[[9]][,ene]),xlim=c(0,8),ylim=c(0,2),col=3)
par(new=T)
plot(density(out4[[9]][,ene]),xlim=c(0,8),ylim=c(0,2),col=4)
par(new=T)
plot(density(out5[[9]][,ene]),xlim=c(0,8),ylim=c(0,2),col=5)
par(new=T)
plot(density(out6[[9]][,ene]),xlim=c(0,8),ylim=c(0,2),col=6)
par(new=T)
plot(density(out7[[9]][,ene]),xlim=c(0,8),ylim=c(0,2),col=7)

#> Natarajan & Kass (1999) in JASA set up ####
#simulatie based on setup of Natarajan & Kass (1999) in JASA
ddsets = 1e3
samsize = 2e4
n = 30
J = 7
beta = c(-.625,.25,-.25,.125)
Psi = matrix(c(.5,0,0,.25),ncol=2)
Psivec = c(Psi)

numPriors = 3
betaMeanMat = array(0,dim=c(ddsets,4,numPriors))
betaHits = matrix(0,nrow=4,ncol=numPriors)
betaWidthMat = array(0,dim=c(ddsets,4,numPriors))
PsiMeanMat = array(0,dim=c(ddsets,4,numPriors))
PsiHits = matrix(0,nrow=4,ncol=numPriors)
PsiWidthMat = array(0,dim=c(ddsets,4,numPriors))
bMatErrorMat = array(0,dim=c(ddsets,2*n,numPriors))
bMatHits = matrix(0,nrow=n*2,ncol=numPriors)
bMatWidthMat = array(0,dim=c(ddsets,2*n,numPriors))

for(dd in 1:ddsets){
  #generate data
  datagen = generate_yvec_logreg(n=n, beta = beta, PsiMat = Psi)
  yvec = datagen[[1]]
  Xmat = datagen[[2]]
  bMat = datagen[[3]]
  bVec = c(bMat)
  
#  out1 = MCMC_mixlogreg_JB_A(yvec,n,J,samsize=samsize)
#  pr = 1
#  out1 = MCMC_mixlogreg_SBeta2_A(yvec,n,J,B0 = 1e3*diag(2),samsize=samsize)
#  pr = 2
  out1 = MCMC_mixlogreg_SBeta2_A(yvec,n,J,B0 = 1,samsize=samsize) #B0 = 1, uses R* for guess
  pr = 3
  
  #store output
  betaMeanMat[dd,,pr] = out1[[1]]
  bMatErrorMat[dd,,pr] = out1[[2]] - bVec
  PsiMeanMat[dd,,pr] = out1[[3]]
  
  betaHits[,pr] = betaHits[,pr] + (out1[[4]][,2]>beta)*(out1[[4]][,1]<beta)
  bMatHits[,pr] = bMatHits[,pr] + (out1[[5]][,2]>bVec)*(out1[[5]][,1]<bVec)
  PsiHits[,pr] = PsiHits[,pr] + (out1[[6]][,2]>c(Psi))*(out1[[6]][,1]<c(Psi))
  
  betaWidthMat[dd,,pr] = out1[[4]][,2]-out1[[4]][,1]
  bMatWidthMat[dd,,pr] = out1[[5]][,2]-out1[[5]][,1]
  PsiWidthMat[dd,,pr] =  out1[[6]][,2]-out1[[6]][,1]

  print(list(dd,Sys.time()))
}

#> Huang & Wand (2013) vs Matrix-F R* prior ####
#Extra simulation for revision using the prior by Huang & Wand (2013)
#use multi cores for simulation
library(doParallel)
# Find out how many cores are available (if you don't already know)
detectCores()
# Create cluster with desired number of cores
cl <- makeCluster(3) # Register cluster
registerDoParallel(cl)
# Find out how many cores are being used
getDoParWorkers()

#simulatie based on setup of Natarajan & Kass (1999) in JASA
ddsets = 1e3
samsize = 2e4
n = 30
J = 7
beta = c(-.625,.25,-.25,.125)
Psi = matrix(c(.5,0,0,.25),ncol=2)
Psivec = c(Psi)

#NOT IN PARALLEL
betaMeanMat = array(0,dim=c(ddsets,4))
betaHits = matrix(0,nrow=4,ncol=1)
betaWidthMat = array(0,dim=c(ddsets,4))
PsiMeanMat = array(0,dim=c(ddsets,4))
PsiHits = matrix(0,nrow=4,ncol=1)
PsiWidthMat = array(0,dim=c(ddsets,4))
bMatErrorMat = array(0,dim=c(ddsets,2*n))
bMatHits = matrix(0,nrow=n*2,ncol=1)
bMatWidthMat = array(0,dim=c(ddsets,2*n))

storeAll = matrix(0,nrow=ddsets,ncol=3*4+3*60+3*4)
for(dd in 1:ddsets){
  #generate data
  datagen = generate_yvec_logreg(n=n, beta=beta, PsiMat=Psi)
  yvec = datagen[[1]]
  Xmat = datagen[[2]]
  bMat = datagen[[3]]
  bVec = c(bMat)
  
  out1 = MCMC_mixlogreg_HW_A(yvec,n,J,samsize=samsize)

  #store output
  storeAll[dd,1:4] = out1[[1]]
  storeAll[dd,4+(1:60)] = out1[[2]] - bVec
  storeAll[dd,64+(1:4)] = out1[[3]]
  
  storeAll[dd,68+(1:4)] = (out1[[4]][,2]>beta)*(out1[[4]][,1]<beta)
  storeAll[dd,72+(1:60)] = (out1[[5]][,2]>bVec)*(out1[[5]][,1]<bVec)
  storeAll[dd,132+(1:4)] = (out1[[6]][,2]>c(Psi))*(out1[[6]][,1]<c(Psi))
  
  storeAll[dd,136+(1:4)] = out1[[4]][,2]-out1[[4]][,1]
  storeAll[dd,140+(1:60)] = out1[[5]][,2]-out1[[5]][,1]
  storeAll[dd,200+(1:4)] =  out1[[6]][,2]-out1[[6]][,1]
  
  print(list(dd,Sys.time()))
}

#IN PARALLEL
genfit_logreg = function(samsize = 2e4){
  n = 30
  J = 7
  beta = c(-.625,.25,-.25,.125)
  Psi = matrix(c(.5,0,0,.25),ncol=2)
  Psivec = c(Psi)
  
  datagen = generate_yvec_logreg(n=n, beta=beta, PsiMat=Psi)
  yvec = datagen[[1]]
  Xmat = datagen[[2]]
  bMat = datagen[[3]]
  bVec = c(bMat)
  
  #out1 = MCMC_mixlogreg_HW_A(yvec,n,J,samsize=samsize)         #inv-Wish following HW
  out1 = MCMC_mixlogreg_SBeta2_A(yvec,n,J,B0=1,samsize=samsize) #matrix-F with B = R*

  #store output
  store1 = rep(0,204)
  store1[1:4] = out1[[1]]
  store1[4+(1:60)] = out1[[2]] - bVec
  store1[64+(1:4)] = out1[[3]]
  
  store1[68+(1:4)] = (out1[[4]][,2]>beta)*(out1[[4]][,1]<beta)
  store1[72+(1:60)] = (out1[[5]][,2]>bVec)*(out1[[5]][,1]<bVec)
  store1[132+(1:4)] = (out1[[6]][,2]>c(Psi))*(out1[[6]][,1]<c(Psi))
  
  store1[136+(1:4)] = out1[[4]][,2]-out1[[4]][,1]
  store1[140+(1:60)] = out1[[5]][,2]-out1[[5]][,1]
  store1[200+(1:4)] =  out1[[6]][,2]-out1[[6]][,1]
  return(store1)
}
print(Sys.time())
resultsF_EB = t(matrix(unlist(mclapply(1:1000,function(ss){
  genfit_logreg(samsize = 2e4)
})),nrow=204))
print(Sys.time())
#resultsHWbackup = resultsHW
#resultsF_EBbackup = resultsF_EB

# Results -----------------------------------------------------------------

results = resultsF_EB

#risk beta
diffs = results[,1:4] - rep(1,1000)%*%t(beta)
mean(apply(diffs**2,1,sum))
sd(apply(diffs**2,1,sum))/sqrt(1000)
#noncoverage beta
1-apply(results[,68+(1:4)],2,mean)
#interval width beta
apply(results[,136+(1:4)],2,mean)

#risk random effects b
mean(apply(results[,4+(1:30)]**2,1,sum))
sd(apply(results[,4+(1:30)]**2,1,sum))/sqrt(1000)
mean(apply(results[,4+(31:60)]**2,1,sum))
sd(apply(results[,4+(31:60)]**2,1,sum))/sqrt(1000)
#noncoverage b
mean(1-apply(results[,72+(1:30)],2,mean))
mean(1-apply(results[,72+(31:60)],2,mean))
#interval width b
mean(results[,140+(1:30)])
mean(results[,140+(31:60)])

#risk Psi
errorPsi = rep(0,1000)
for(ii in 1:1000){
  invPsi = solve(Psi)
  diffMat = matrix(results[ii,64+(1:4)],ncol=2)%*%invPsi-diag(2)
  errorPsi[ii] = sum(diag(diffMat%*%diffMat))
}
mean(errorPsi)
sd(errorPsi)/sqrt(1000)
#noncoverage Psi
1-apply(results[,132+(1:4)],2,mean)
#interval width Psi
apply(results[,200+(1:4)],2,mean)



