### Project:     Master Thesis
### Object:      MCMC functions for model fitting
### Description: Contains the functions used to fit all the models of interest
### Date:        2019-03-17

# Normal model:
# y_{ij} = normal(p_{ij})
# y_{ij} = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)

# # Set up
#   post_draws_fun_filename <- "./R/190313-normalmod-functions.R"
#   source(post_draws_fun_filename)
# 
# # Data
# # Depression Scale (n = 46, J = 6)
#   RiesbyDat <- read.table("./data/RiesbyDat.txt")
# 
#   # Make it digestable for the functions
#   yvec     <- RiesbyDat$depr
#   Xmat     <- cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog)
#   J        <- length(unique(RiesbyDat$week))
#   n        <- nrow(RiesbyDat)/J
#   Zi       <- cbind(rep(1, length = J), unique(Xmat[, 2]))
# 
# # Test Functions
#   out1 <- MCMC_HWprior(yvec = RiesbyDat$depr,
#                        Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
#                        J    = length(unique(RiesbyDat$week)),
#                        n    = nrow(RiesbyDat)/J,
#                        Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
#                        iniv = 0,
#                        samsize = 2000,
#                        burnin = 1/10)
#   out2 <- MCMC_invWish(yvec = RiesbyDat$depr,
#                        Xmat = cbind(rep(1, nrow(RiesbyDat)), RiesbyDat$week, RiesbyDat$endog, RiesbyDat$week*RiesbyDat$endog),
#                        J    = length(unique(RiesbyDat$week)),
#                        n    = nrow(RiesbyDat)/J,
#                        Zi   = cbind(rep(1, length = J), unique(Xmat[, 2])),
#                        iniv = 0,
#                        B0   = diag(2),
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
#   
# MCMC Functions ----------------------------------------------------------

# Model fitting with HW prior
    
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
          Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
          PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
          sigma2 <- attr(VarCorr(fit), "sc")
          bMat   <- as.matrix(ranef(fit)$cluster)
        } else {
          #theta  <- t(rep(0, ncol(Xmat)))
          Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
          PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
          sigma2 <- 1 #attr(VarCorr(fit), "sc")
          bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
        }
        avec   <- rep(100,2)
        
      # Run MCMC
        PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
        PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
        PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
        PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr
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
      
          PD_theta[ss,] = theta
          PD_bMat[ss,]  = c(bMat)
          PD_avec[ss,]  = avec
          PD_Psi[ss,]   = c(solve(PsiInv))
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
        return(list(PD_theta,  #1
                    PD_bMat,   #2
                    PD_Psi,    #3
                    PD_Psi_sd, #4
                    PD_avec))  #5
    }
    
# Model fitting with inv-Wish standard prior
    
    MCMC_invWish = function(yvec, Xmat, Zi, J, n, iniv = 1, B0 = diag(2), samsize = 2000, burnin = 1/10){
      # Data Prep
        xvec    <- Xmat[, 2]
        covar   <- Xmat[, 3]
        inter   <- Xmat[, 4]
        cluster <- rep(seq(1, n), each = J)
      # Define priors
        bi0   <- rep(0, ncol(Zi)) # prior for bi
        S0    <- B0               # prior guess for Psi
      # Define initial Values
        if (iniv == 1){
          fit   <- lmer(yvec ~ xvec + covar + inter + (1 + xvec | cluster), REML = FALSE)
          #theta  <- fixef(fit)
          Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
          PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
          sigma2 <- attr(VarCorr(fit), "sc")
          bMat   <- as.matrix(ranef(fit)$cluster)
        } else {
          #theta  <- t(rep(0, ncol(Xmat)))
          Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
          PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
          sigma2 <- 1 #attr(VarCorr(fit), "sc")
          bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
        }
        
      # Run MCMC
        PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
        PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
        PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
        PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr

        #burn-in
        for(ss in 1:(samsize*burnin)){
          print(paste("Burn-in:", ss))
          
          theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
          bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
          sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
          PsiInv   <- draw_PsiInv_InvWish(n,bMat,S0 = S0)
      
          # PD_theta[ss,] = theta
          # PD_bMat[ss,]  = c(bMat)
          # PD_Psi[ss,]   = c(solve(PsiInv))
        }
        #sampling
        for(ss in 1:samsize){
          print(paste("Post-Sample:", ss))
          
          theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
          bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
          sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
          PsiInv   <- draw_PsiInv_InvWish(n,bMat,S0 = S0)
          
          PD_theta[ss,] = theta
          PD_bMat[ss,]  = c(bMat)
          PD_Psi[ss,]   = c(solve(PsiInv))
            sdRI    <- sqrt(PD_Psi[ss,1])       # save the satandard deviation version
            sdRS    <- sqrt(PD_Psi[ss,4])
            corRIRS <- PD_Psi[ss,2]/(sdRI*sdRS) # save the correlation
          PD_Psi_sd[ss,] = c(sdRI, corRIRS, corRIRS, sdRS)
        }
        return(list(PD_theta,  #1
                    PD_bMat,   #2
                    PD_Psi,    #3
                    PD_Psi_sd))#4
    }
    
# Model fitting with inv-Wish standard prior
    
    MCMC_matF = function(yvec, Xmat, Zi, J, n, iniv = 1, B0 = 1e3*diag(2), samsize = 2000, burnin = 1/10){
      # Data Prep
        xvec    <- Xmat[, 2]
        covar   <- Xmat[, 3]
        inter   <- Xmat[, 4]
        cluster <- rep(seq(1, n), each = J)
      # Define priors
        bi0   <- rep(0, ncol(Zi)) # prior for bi
        B0Inv <- solve(B0) 
      # Define initial Values
        if (iniv == 1){
          fit   <- lmer(yvec ~ xvec + covar + inter + (1 + xvec | cluster), REML = FALSE)
          #theta  <- fixef(fit)
          Psi    <- matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
          PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
          sigma2 <- attr(VarCorr(fit), "sc")
          bMat   <- as.matrix(ranef(fit)$cluster)
        } else {
          #theta  <- t(rep(0, ncol(Xmat)))
          Psi    <- matrix(c(1,0,0,1), ncol = 2) #matrix(VarCorr(fit)[[1]][1:4], ncol = 2)
          PsiInv <- solve(Psi)        # var-covar matrix of random effects estimated with glmer
          sigma2 <- 1 #attr(VarCorr(fit), "sc")
          bMat   <- matrix(rep(0, n*2), ncol = 2) #as.matrix(ranef(fit)$cluster)
        }
          Omega = B0

      # Run MCMC
        PD_theta  = matrix(0, nrow = samsize, ncol = ncol(Xmat))
        PD_bMat   = matrix(0, nrow = samsize, ncol = n*2)
        PD_Omg    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2)
        PD_Psi    = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS var and cov
        PD_Psi_sd = matrix(0, nrow = samsize, ncol = ncol(Zi)*2) # store RI, RS sd and corr

        #burn-in
        for(ss in 1:(samsize*burnin)){
          print(paste("Burn-in:", ss))
          
          theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
          bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
          sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
          outDummy <- draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv,n)
              PsiInv <- outDummy[[1]]
              Omega  <- outDummy[[2]]
      
          # PD_theta[ss,] = theta
          # PD_bMat[ss,]  = c(bMat)
          # PD_Psi[ss,]   = c(solve(PsiInv))
        }
        #sampling
        for(ss in 1:samsize){
          print(paste("Post-Sample:", ss))
          
          theta    <- draw_theta(yvec,Xmat,Zi,bMat,sigma2)
          bMat     <- draw_bMat(yvec,Xmat,Zi,theta,bi0,bMat,sigma2,PsiInv,n,J)
          sigma2   <- draw_sigam2_IMPprior(yvec,Xmat,Zi,theta,bMat,n,J)
          outDummy <- draw_PsiInv_matF(yvec,Xmat,Zi,bMat,PsiInv,Omega,B0Inv,n)
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
        return(list(PD_theta,  #1
                    PD_bMat,   #2
                    PD_Psi,    #3
                    PD_Psi_sd, #4
                    PD_Omg))   #5 
    }
    