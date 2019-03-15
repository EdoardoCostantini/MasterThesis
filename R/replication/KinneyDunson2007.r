#FIXED AND RANDOM EFFECTS SELECTION IN LINEAR AND LOGISTIC MODELS (Kinney and Dunson, Biometrics, forthcoming)
#R code for sampling functions for full conditional posterior distributions given in Appendix
#this code is not black-box ready-to-run, or supported, or optimized for computational efficiency,
#it is provided simply to help other researchers with their own Gibbs samplers.
#uses libraries mvtnorm and MASS

# Definitions
# N = numbers obs
# n = number groups
# group = grouping variable
# XA = potential fixed effects matrix
# X = fixed effects matrix for current model
# Z = random effects matrix (I set Z=XA)
# p=ncol(X); q=ncol(Z); r = q*(q-1)/2
# fixc = X%*%beta #fixed effects component
# randc = random effects component

#parameters updated
# w = latent variable
# phi = vector of precisions, see eq 8
# beta = fixed effects coefficients
# J = vector of indicators for fixed effects in current model
# lambda = diagonal elements of diagonal matrix Lambda, if parameter expansion not used. In context here lambda is actually diagonal elements of A.
# Gamma = covariance component Gamma (matrix)
# b = random effects vector, in context here with parameter expanded model it is called xi
# DA = matrix D
# c = Zellner prior parameter
# sig2 = sigma^2 (normal case)

## RECOMMENDED PRIOR PARAMETER VALUES ##
ac = .5; bc = N/2 # prior pars for c
lam0 = rep(0,q); lV0 = rep(1,q) # prior pars lambda
gam0 = rep(0,r); gV0 = diag(.5,r) #prior gamma
pb0 = rep(.5,p) #prior prob beta
pl0 = rep(.5,q) #prior prob lam
ad = .5; bd = N/2 # prior pars d

nu = 7.3; sig2f = pi^2*(nu-2)/(3*nu) #fixed values for t approximation

#HELPER FUNCTIONS
bind = numeric(r);s=1;e=1;l=1;  #used for updating Gamma, only need to run once
  zind = numeric(r)
  for(m in 1:(q-1)){
    bind[s:e]= 1:m; 
    zind[s:e]= rep(m+1,m)
    s=e+1; e=s+l; l=l+1
  }

#calculate random effects component vector
rand = function(N,n,group,Z,lamdba,Gamma,b){
  randc = numeric(N)
  for(i in 1:n){
    randc[group==i]= Z[group==i,]%*%diag(lambda)%*%Gamma%*%b[i,]}
  randc
}

#used for lambda and gamma
mbij = function(q,N,group,b){
  bij = matrix(0,q,N)
  for(i in 1:n){bij[,group==i]= b[i,]}
  bij = t(bij)
}

#PARAMETER UPDATE FUNCTIONS
# update w
wsam=function(fixc,randc,y,sig2){
eta = fixc + randc
    w[y==0] <- qnorm(runif(sum(1-y),0, pnorm(0,eta[y==0],sqrt(sig2[y==0]))),eta[y==0],sqrt(sig2[y==0]))
    w[y==1] <- qnorm(runif(sum(y),  pnorm(0,eta[y==1],sqrt(sig2[y==1])),1),eta[y==1],sqrt(sig2[y==1]))
       w}

# update phi
phisam= function(nu,w,fixc,randc,sig2f){
   phia = (1+nu)/2
   phib = (w-randc-fixc)^2/(2*sig2f) + nu/2
   phi = rgamma(N,phia,phib)

   sig2=sig2f/phi
   sig2}

## sample beta #zellner prior version
  #first need to change dim dep on J
  Js = (J*1:p)[J>0]
  X = as.matrix(XA[,Js])
   xst = X*sqrt(phi/sig2f + 1/c)
   betav = solve(t(xst)%*%xst)
   betam = apply(phi*as.vector(w-randc)*X,2,sum)%*%betav/sig2f
   beta = as.matrix(mvrnorm(1,betam,betav))

  # expand to full dimension, for storage purposes, otherwise not needed
  betaf = numeric(p)
  betaf[Js] = beta  
  # update fixc
  fixc = X%*%beta
  

## sample lambda
  bij = mbij(q,N,group,b)
  lambda = lamsamlog(p,q,N,n,group,b,Gamma,fixc,sig2f,lam0,lV0,pl0,bij,phi)
  lamsamlog = function(p,q,N,n,group,b,Gamma,fixc,sig2f,lam0,lV0,pl0,bij,phi){
  t=matrix(0,N,q)
  for(k in 1:N){
      t[k,1] = Z[k,1]*bij[k,1]
    for(l in 2:q){
      t[k,l] = Z[k,l]*(bij[k,l] + bij[k,1:(l-1)]%*%Gamma[l,1:(l-1)])
    }
  }
  for(l in 1:q){
    v<- 1/(1/lV0[l] + sum(phi*t[,l]^2)/sig2f)
    Tij = w - fixc - t[,-l]%*%lambda[-l]
    E = v*(sum(phi*t[,l]*Tij)/sig2f + lam0[l]/lV0[l])
    bkn = dnorm(0,lam0[l],sqrt(lV0[l]))*(1-pnorm(0,E,sqrt(v)))
    bkd = dnorm(0,E,sqrt(v))*(1-pnorm(0,lam0[l],sqrt(lV0[l])))
    phat = pl0[l]/(pl0[l] + (1-pl0[l])*bkn/bkd)
 
    m<- 1-rbinom(1,1,phat)
    lambda[l]<- m*qnorm(runif(1,pnorm(0,E,sqrt(v)),1),E,sqrt(v))
  }
  lambda}

# sample Gamma
  Gamma = gamsamlog(p,q,r,n,N,bind,zind,lambda,Z,fixc,gV0,sig2f,w,bij,phi)
gamsamlog = function(p,q,r,n,N,bind,zind,lambda,Z,fixc,gV0,sig2f,w,bij,phi){
  u = matrix(0,N,r)               
  for(k in 1:N){
    for(t in 1:r){
      u[k,t]=bij[k,bind[t]]*lambda[zind[t]]*Z[k,zind[t]]
    }}
  us = sqrt(phi)*u
  Vgam = solve(t(us)%*%us/sig2f + solve(gV0))
  Egam = apply(phi*as.vector(w-fixc)*u,2,sum)%*%Vgam/sig2f

  gamma = mvrnorm(1,Egam,Vgam)
  Gamma=diag(q)
  Gamma[upper.tri(Gamma)]=gamma
  Gamma=t(Gamma)
  for(k in 1:q){  #apply 1(Rlam)
    if(lambda[k]==0){
      Gamma[k,1:(k-1)]= 0
      if(k<q){
      Gamma[(k+1):q,k]= 0}
     }
   }
  Gamma[1,1]=1
  Gamma}

# sample d_l             ##updated only in parameter expanded model
 DA=dsam(q,n,ad,bd,b) 
dsam = function(q,n,ad,bd,b){
  for(l in 1:q){DA[l,l]=1/rgamma(1,ad + n/2, bd + t(b[,l])%*%b[,l]/2)}
  DA}

# sample b
  b = bsamlog(lambda,n,phi,group,Gamma,sig2f,DA,fixc,Z)
  randc = rand(N,n,group,Z,lamdba,Gamma,b) #update random component
bsamlog=function(lambda,n,phi,group,Gamma,sig2f,DA=diag(q),fixc,Z){
  K = ((I(lambda!=0))*1:q)[lambda>0]
  qs = length(K)
  if(qs>0){
  bs = matrix(0,n,qs)
    for(i in 1:n){
      vb = sqrt(phi[group==i])*as.matrix(Z[group==i,K])%*%as.matrix(diag(lambda)[K,K])%*%Gamm
a[K,K]
      Vb = solve(t(vb)%*%vb/sig2f + solve(DA[K,K]))
      Eb = t(as.matrix(apply(sqrt(phi[group==i])*as.matrix(w-fixc)[group==i]*vb,2,sum)))%*%Vb
/sig2f
      bs[i,] = mvrnorm(1,Eb,Vb)
    }
    b[,K] = bs
  }
  b}

# sample c
   c1 = (sum(J) + 1)/2
   c2 = t(beta)%*%t(X)%*%X%*%beta/2 + N/2
   c = 1/rgamma(1,c1,c2)

# sample J
  J = Jsamlog(p,J,c,pb0,w,randc,phi,XA,sig2f)

Jsamlog = function(p,J,c,pb0,w,randc,phi,XA,sig2f){
  for(k in 2:p){
    J0=J; J0[k]=0
    J1 = J; J1[k]=1
    hk=0
    if(sum(J)>1){
      hk = sqrt(c)*exp(-.5*(SJL1(J0,w,randc,phi,XA,sig2f,c)-SJL1(J1,w,randc,phi,XA,sig2f,c)))
*(SJL2(J0,XA,phi,sig2f,c)/SJL2(J1,XA,phi,sig2f,c))*(1-pb0[k])/pb0[k]}
    pk = 1/(1+hk)
    J[k] = rbinom(1,1,pk)
  }
  J}

SJL1=function(R,w,randc,phi,XA,sig2f,c){
   Js = (R*1:p)[R>0]
   X = as.matrix(XA[,Js])
   xst = X*sqrt(phi/sig2f + 1/c)
   betav = solve(t(xst)%*%xst)
   betam = apply(phi*as.vector(w-randc)*X,2,sum)%*%betav/sig2f
   S = sum(phi*(w-randc)^2) - betam%*%solve(betav)%*%t(betam)
   return(S)}
SJL2=function(R,XA,phi,sig2f,c){
   Js = (R*1:p)[R>0]
   X = as.matrix(XA[,Js])
   xst = X*sqrt(phi/sig2f + 1/c)
   betav = solve(t(xst)%*%xst)
   S=sqrt(det(t(X)%*%X%*%betav))
   S}

#for normal case
sigsam = function(J,N,fixc,randc,w,beta,X,c){
  siga = sum(J)/2 + N/2
  sigb = (t(w-fixc-randc)%*%(w-fixc-randc))/2 + t(beta)%*%t(X)%*%X%*%beta/(2*c)
  sig2 = 1/rgamma(1,siga,sigb)
  sig2}

csam = function(J,ac,bc,sig2,X,beta){
  ca = sum(J)/2 +ac
  cb = t(beta)%*%t(X)%*%X%*%beta/(2*sig2) + bc 
  c = 1/rgamma(1,ca,cb)
  c}
Jsam = function(p,J,c,pb0,w,randc,XA){
  for(k in 2:p){
    J0 = J; J0[k]=0
    J1 = J; J1[k]=1
    hk=0
    if(sum(J)>1){
    hk = sqrt(c+1)*(SJ(J0,w,randc,XA,c)/SJ(J1,w,randc,XA,c))^(-N/2)*(1-pb0[k])/pb0[k]}
    pk = 1/(1+hk)
    J[k] = rbinom(1,1,pk)
  }
  J}
SJ=function(R,w,randc,XA,c){
    s = (R*1:p)[R>0]
    Xs = XA[,s]
    S = t(w-randc)%*%(w-randc) - t(w-randc)%*%Xs%*%solve(t(Xs)%*%Xs)%*%t(Xs)%*%(w-randc)*c/(c
+1)
    return(S)}
bsam = function(n,q,lambda,Gamma,sig2,X,beta,w,group,DA=diag(q)){

  K = ((I(lambda!=0))*1:q)[lambda>0]
  qs = length(K)
  if(qs>0){
  bs = matrix(0,n,qs)
    for(i in 1:n){
      vb = Z[group==i,K]%*%as.matrix(diag(lambda)[K,K]%*%Gamma[K,K])
      Vb = solve(t(vb)%*%vb/sig2 + solve(DA[K,K]))
      Eb = Vb%*%(t(vb)%*%(w[group==i]-X[group==i,]%*%beta))/sig2
      bs[i,] = mvrnorm(1,Eb,Vb)
    }
    b[,K] = bs
  }
  b}
