### Project:     Master Thesis
### Object:      MCMC simulation (fitting the model)
### Description: Contains the functions and specification to fit the bayesian model as defined in
###              model.pdf file in text/.
### Date:        2019-03-15

# Load some data
  RiesbyDat <- read.table("./data/RiesbyDat.txt")
  
  # Make it digestable for the functions
  yvec     <- RiesbyDat$depr
  Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
  J        <- length(unique(RiesbyDat$week))
  n        <- nrow(RiesbyDat)/J
  Zi       <- cbind(rep(1,length = J), unique(Xmat[, 2]))
  subjects <- RiesbyDat$id

# Define priors
  B0       <- 10**3*diag(2)
  B0Inv    <- solve(B0)
  nu0      <- .001
  sigma20  <- .001
  bi0      <- rep(0, ncol(Zi))
  theta0   <- rep(0, ncol(Xmat)) #c(22, -2, 2, -.2) #rep(0, ncol(Xmat))
  Sigma0   <- t(Xmat) %*% Xmat/(n*3)       # Think about what should this guess be

# Define initial Values
  fit <- lmer(depr ~ week + endog + inter + (1 + week | id), data = RiesbyDat, REML = FALSE)
  
  theta  <- c(10, 10, 10, 10) #fixef(fit)
  Psi    <- matrix(c(5,0,0,5), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
  PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
  sigma2 <- 1 #attr(VarCorr(fit), "sc")
  bMat   <- as.matrix(ranef(fit)$id)
  Omega  <- B0
  avec = rep(100,2)
  
# Run MCMC
  samsize <- 2000
  burnin  <- 1/10
  MbetaA   = matrix(0,nrow=samsize,ncol=ncol(Xmat))
  MPsiA    = matrix(0,nrow=samsize,ncol=ncol(Zi)*2)
  MOmegaA  = matrix(0,nrow=samsize,ncol=ncol(Zi)*2)
  MbMatA   = matrix(0,nrow=samsize,ncol=n*2)
  MavecA = matrix(0,nrow=samsize,ncol=2)
  
  #MimpwA   = matrix(0,nrow=samsize,ncol=1)
  
  #burn-in
  for(ss in 1:(samsize*burnin)){
    print(paste("Burn-in:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,theta0,Sigma0,bMat)
    bMat     <- draw_bMat_A(yvec,Xmat,Zi,theta,bi0,bMat)
    sigma2   <- draw_sigam2(nu0, sigma20)
    # outDummy <- draw_PsiInv_MSBeta2(yvec,Xmat,Zi,PsiInv,Omega,B0Inv)
    #   PsiInv <- outDummy[[1]]
    #   Omega  <- outDummy[[2]]
    outDummy = draw_PsiInv_HW(PsiInv,avec,bMat)
      PsiInv = outDummy[[1]]
      avec = outDummy[[2]]
    
    MbetaA[ss,]   = theta
    MbMatA[ss,]   = c(bMat)
    # MPsiA[ss,]    = c(solve(PsiInv))
    # MOmegaA[ss,]  = Omega
    MPsiA[ss,] = c(solve(PsiInv))
    MavecA[ss,] = avec
  }
  for(ss in 1:samsize){
    print(paste("Post-Sample:", ss))
    
    theta    <- draw_theta(yvec,Xmat,Zi,theta0,Sigma0,bMat)
    bMat     <- draw_bMat_A(yvec,Xmat,Zi,theta,bi0,bMat)
    sigma2   <- draw_sigam2(nu0, sigma20)
    # outDummy <- draw_PsiInv_MSBeta2(yvec,Xmat,Zi,PsiInv,Omega,B0Inv)
    #   PsiInv <- outDummy[[1]]
    #   Omega  <- outDummy[[2]]
    outDummy = draw_PsiInv_HW(PsiInv,avec,bMat)
      PsiInv = outDummy[[1]]
      avec = outDummy[[2]]
    
    MbetaA[ss,]   = theta
    MbMatA[ss,]   = c(bMat)
    # MPsiA[ss,]    = c(solve(PsiInv))
    # MOmegaA[ss,]  = Omega
    MPsiA[ss,] = c(solve(PsiInv))
    MavecA[ss,] = avec
  }
  #out1 <- list(MbetaA,MbMatA,MPsiA,MOmegaA)
  out1 <- list(MbetaA,MbMatA,MPsiA,MavecA)
  
  which_o <- 1
  which_c <- 1
  plot(1:samsize,out1[[which_o]][,which_c],"l")
  mean(out1[[which_o]][,which_c])
  median(out1[[which_o]][,which_c])
  # Compare with lme results
  fit
  hist(out1[[which_o]][,which_c], breaks = 100,
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_o]][,which_c]), col = "blue", lwd = 2)
      abline(v = median(out1[[which_o]][,which_c]), col = "red", lwd = 2)
  
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  
      
      
  hist(out2[[which_n]][[9]][,1], breaks = 100,
       main = "mat-F w/ 1e3*I",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 2)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 2)
  hist(out3[[which_n]][[9]][,1], breaks = 100,
       main = "inv-Wish HW2013",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 2)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 2)
      
      
      
      # Scrop notes for starting values
      
        theta  <- t(rep(0, ncol(Xmat))) #fixef(fit)
        Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
        PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
        sigma2 <- 1 #attr(VarCorr(fit), "sc")
        bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
        