# Title:       Replicate Hoff 2009 Chapter 7
# Description: GLM Mixed effects with inverse-Wishart distribution for the var-cov matrix

library(mvtnorm)
library(MCMCpack)

# 11.1 A hierarchical regression model ------------------------------------
# p199 (pdf), p.195 (book)
#> Data Preparation ####
  odat <- dget("./Data/nels_2002.txt")
  colnames(odat)<-c("sch_id","sch_enroll","sch_freelunch","sch_cnrtl",
                    "sch_urban","mteach_deg","eteach_deg","mteach_years","eteach_years" , 
                    "stu_sex","stu_lang","stu_pared","stu_income","stu_mathscore",
                    "stu_readscore","stu_mhw","stu_ehw","stu_readhours","stu_ses")
  odat <- as.data.frame(odat)
  t(head(odat, 3))
  
  # Perform Some Form of selection (read book again for insights)
  ids <- dget("./Data/ids_selectschools.txt")
  group <- odat$sch_id
  indset <- apply( odat[, 1, drop = FALSE],
                   1,
                   is.element,
                   ids)
  dat <- odat[indset,]
  mathdat <- dat[,c(1, 3, 19, 14)]
  mathdat[, 3] <- (mathdat[,3] - mean(mathdat[,3])) / sd(mathdat[,3])
  #head(mathdat)
  #groups <- ids
  m <- length(ids)
  Y <- list() ; X <- list() ; N <- NULL
  for(j in 1:m) {
    Y[[j]] <- mathdat[mathdat[, 1] == ids[j], 4] # students' scores in group j
    N[j]   <- sum(dat$sch_id == ids[j])          # number of students in group j
    xj     <- mathdat[mathdat[, 1] == ids[j], 3] # students' SES in group j
    xj     <- (xj - mean(xj))                    # (centered)
    X[[j]] <- cbind( rep(1, N[j]), xj  )         # design matrix (intercept, slope 1, ..., slope p)
                                                 # for group j!
  }

  # Create the OLS fits
  BETA.LS <- NULL # OLS coefficients estiamtes (intercept and slopes)
  S2.LS   <- NULL # OLS with-in group residual variance
  for(j in 1:m) {
    fit     <- lm(Y[[j]] ~ -1 + X[[j]] ) # we get rid of the intercept because there is a 
                                         # column of 1s in our desing matrix that already estiamtes
                                         # the intercept
    BETA.LS <- rbind(BETA.LS, c(fit$coef)) 
    S2.LS   <- c(S2.LS, summary(fit)$sigma^2)
  }
  
#> Plot OLS ####
  # Plot Least squared regression lines for the math score data,
  # and plots of estimates versus group sample size
  par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
  par(mfrow=c(1,3))
  
  # Fig1 - math score vs SES
    plot( range(mathdat[,3]), range(mathdat[,4]), # create the space, no dots
          type = "n", xlab = "SES", ylab = "math score")
    for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  } # add regression lines
    BETA.MLS <- apply(BETA.LS,2,mean)
    abline(BETA.MLS[1], BETA.MLS[2], lwd=2)
  # Fig2 - intercept vs sample size
    plot(N,BETA.LS[, 1],
         xlab = "sample size", ylab = "intercept")
    abline(h = BETA.MLS[1], col = "black", lwd = 2)
  # Fig3 - slope vs sample size
    plot(N,BETA.LS[, 2],
         xlab = "sample size", ylab = "slope")
    abline(h = BETA.MLS[2], col = "black", lwd = 2)

  #dev.off()

# Hierarchical Model 
  p <- dim(X[[1]])[2] # number of coefficients (intercept + number of predictors)
#> Priors ####
  # theta_vec ~ multivariate-Norm(mu0_vec, L0)
    mu0    <- colMeans(BETA.LS) # use intercept and coefficients OLS estimates for prior guess
    L0     <- cov(BETA.LS)
    iL0    <- solve(L0)
  # sigma2 ~ invGamma(nu0/2, nu0*s20/2)
    nu0    <- 1
    s20    <- mean(S2.LS) # prior guess for sigma2 (error variance)
  # Sigma ~ invWish(eta0, S0)
    eta0   <- p + 2
    S0     <- cov(BETA.LS) # prior guess for Sigma
    
# Initial values and other objs
  theta  <- colMeans(BETA.LS) # working theta for full conditional computations (initial value)
  s2     <- mean(S2.LS)  # working sigma2 for full conditional computations (initial value)
  Sigma  <- cov(BETA.LS) # working Sigma for full conditional computations (initial values)
    iSigma <- solve(Sigma)

# Prepare OBJs for MCMC
  REPS <- 1000 # how many iterations
  SAVES <- 10 # save the result of which iterations? 1 in 10? then this is equal 10
  BETA <- BETA.LS # copy OLS beta in a BETA object for the dimensionality
                  # this object will keep temporary updates of betas
  # Create empty objects to store posterior draws
  
  # Original
  # S2.b     <- NULL # posterior draws of sigma2
  # THETA.b  <- NULL # posterior draws of theta
  # Sigma.ps <- matrix(0, p, p) # posterior draws of Sigma
  # SIGMA.PS <- NULL # posterior draws of SIGMA?
  # BETA.ps  <- BETA*0 # posterior draws (samples) of betas (posterior distirbution)
  
  # My method
  THETA.b       <- matrix(rep(NA, p*REPS/SAVES), ncol = p)
  S2.b          <- rep(NA, REPS/SAVES)
  SIGMA.PS      <- matrix(rep(NA, 4*REPS/SAVES), ncol = 4)
  BETA.int.ps   <- matrix(rep(NA, length(N)*1*REPS/SAVES), ncol = length(N))
  BETA.coefs.ps <- matrix(rep(NA, length(N)*(p-1)*REPS/SAVES), ncol = length(N))

  head(t(BETA.int.ps))
  BETA.pp  <- NULL # posterior predictive distribution of betas for a to-be-sampled school
  set.seed(1)
  
  start <- Sys.time()
  
# MCMC
for(s in 1:REPS) {
  ##update beta_j (group specific effects)
  for(j in 1:m) {   # for each group
    Vj        <- solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
    Ej        <- Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[[j]]/s2 )
    BETA[j, ] <- rmvnorm(1, Ej, Vj)
  }

  ##update theta (fixed effects: constant across groups, one for each predictor, that's why is a vector)
  Lm    <- solve( iL0 +  m*iSigma )
  mum   <- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA, 2, sum)) # why there is no m ? I have it the full conditonal for theta
  theta <- t(rmvnorm(1, mum, Lm))

  ##update Sigma <- INVERSE WISHART! We are going to work here:)
  mtheta <- matrix(theta, m, p, byrow = TRUE)
  iSigma <- rwish(eta0 + m, solve( S0 + t(BETA-mtheta)%*%(BETA-mtheta) ))

  ##update s2
  RSS <- 0
  for(j in 1:m) { 
    RSS <- RSS + sum( (Y[[j]] - X[[j]]%*%BETA[j,] )^2 ) 
  }
  s2 <-1/rgamma(1,
                (nu0 + sum(N))/2, (nu0*s20 + RSS)/2 )

  ##store results
  if(s%%10 == 0) 
  { 
    cat(s, s2, "\n")
    # My method
    S2.b[s/SAVES]            <- s2
    THETA.b[s/SAVES, ]       <- t(theta)  #fixed effects
    BETA.int.ps[s/SAVES, ]   <- BETA[, 1] 
    BETA.coefs.ps[s/SAVES, ] <- BETA[, 2]
    SIGMA.PS[s/SAVES, ] <- c(solve(iSigma))
    
    # Original
    # S2.b  <- c(S2.b, s2)
    # THETA.b  <- rbind(THETA.b, t(theta)) #fexied effects
    # Sigma.ps <- Sigma.ps + solve(iSigma) #
    # BETA.ps  <- BETA.ps + BETA                    # what is this sum?
    # SIGMA.PS                 <- rbind(SIGMA.PS, c(solve(iSigma)))
    # BETA.pp  <- rbind(BETA.pp, rmvnorm(1, theta, solve(iSigma)) )
  }
}
  
  end <- Sys.time()
  end - start
  
  # Results
  # Point Estimates and Credibility Intervals
    quantile(THETA.b[, 2], prob = c(.025,.5,.975))
    mean(BETA.pp[, 2] < 0)
  
  # Plots
    par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
    par(mfrow=c(1,2))
    
    # Posterior Density of theta.2 of a randomly sampled school
    # with posterior predictive distribution of a randomly sampled slope
    plot(density(THETA.b[, 2], adj = 2),
         xlim = range(BETA.pp[,2]), xlab = "slope parameter", ylab = "posterior density",
         main = "", lwd = 2)
    lines(density(BETA.pp[, 2], adj = 2), col = "gray", lwd = 2)
    legend(-3, 1.0,
           legend = c( expression(theta[2]),
                       expression(tilde(beta)[2])), 
           lwd = c(2,2), col = c("black", "gray"), bty = "n")
    
    # Posterior expectations of the 100 school-specific regression lines,
    # average line given in black
    BETA.PM <- BETA.ps/1000
    plot( range(mathdat[,3]), range(mathdat[,4]),
          type = "n", xlab = "SES", ylab = "math score" )
    for(j in 1:m) { abline(BETA.PM[j, 1], BETA.PM[j, 2], col = "gray") }
    abline( mean(THETA.b[, 1]), mean(THETA.b[,2]), lwd = 2)
    
    dev.off()
    
# 11.4 Generalized linear mixed effects models ----------------------------
# p205 (pdf), p.201 (book)
#> Data (tumor example) ####
  # Load
    XY.tumor <- dget("http://www.stat.washington.edu/~hoff/Book/Data/data/XY.tumor")
  # Structure
    # Mices are units (J = 21)
    nrow(XY.tumor$Y)
    # Observations (nj = 20) are counts of tumor in 20 diff segments of mouse intestine
      # Each mouse (cluster) contains 20 observations (nj = 20)
    ncol(XY.tumor$Y)
    # Total sample size (N = 420)
    ncol(XY.tumor$Y) * nrow(XY.tumor$Y)
    # The desing matrix of each cluster (mouse) contains:
    attributes(XY.tumor$X)
    # - 20 rows (observations)
    # - 5 columns (1 for intercept, 4 for predicotrs,
    #   which is actually 1 predictor in 4 different powers)
    #   the predicotr x, has 20 values between 0 and 1 (.05, .10,..., .95)
    
  # Prep for analysis
  Y <- XY.tumor$Y
  X <- XY.tumor$X
  m <- dim(Y)[1]
  p <- dim(X)[2]
  
#> Priors ####
  # Theta ~ multivaraite-N(mu0, L0)
  BETA <- NULL
  colnames(BETA) <- c("int", "b_x1","b_x2","b_x3","b_x4")
  for(j in 1:m) {
    BETA <- rbind(BETA, lm(log(Y[j, ]+1/20) ~ -1+X[,,j])$coef) # Serves as initial values for beta as well
  }
  mu0 <- apply(BETA, 2, mean) # "prior" guesses (initial values) are mean of the 21 OLS
                              # coefficients estiamted for each group separately
  L0  <- cov(BETA); iL0 <- solve(L0)
  
  # Sigma ~ inverse-Wishart(eta0, S0**-1)
  eta0   <- p + 2             # prior "sample size" (usual rule for the first invwish parameter)
  S0     <- cov(BETA)         # prior guess for Sigma (sample covairance matrix as prior guess)
    # How do I use the MLE estiamte of this? mle4 model? investigate the vcov matrix there 
  
  # Initial Values
  iSigma <- solve(S0)         # initial value for iSigma
  BETA                        # initial values for regression coefiients
  
fit  <- lm(log(Y[1, ]+1/20) ~ -1+X[,,1])
vcov(fit)

#> MCMC ####
set.seed(1)
THETA.post <- NULL
for(s in 1:50000) {
  ## update theta (a mean fixed effect for each predictor)
  Lm    <- solve( iL0 + m*iSigma )
  mum   <- Lm %*% ( iL0%*%mu0 + iSigma%*%apply(BETA, 2, sum) ) # maybe multiply second part by m?
  theta <- t(rmvnorm(1, mum, Lm))
  
  mtheta <- matrix(theta, m, p, byrow = TRUE)
  
  ## update Sigma (use of inve-Wish)
  Stheta <- t(BETA - mtheta)%*%(BETA - mtheta)
  iSigma <- rwish(eta0 + m,
                  solve( S0 + Stheta ))
  
  ## update beta (Metropolis algorithm)
  Sigma  <- solve(iSigma)
  Vj <-.5*Sigma # a scaled version of Sigma^(s) for the proposal dist of Bj is a safe bet,
                # often need trial and error to produce a good mixing MC
  dSigma <- det(Sigma)
  for(j in 1:m) {
    # Sample of beta.j from proposal distirbution
    beta.p <- t( rmvnorm(1,
                         mean  = BETA[j ,],
                         sigma = Vj) )
    
    # Acceptance Ratio
    lr <- sum( dpois(Y[j,], exp(X[,,j]%*%beta.p),
                     log = TRUE ) -
                 dpois(Y[j,], exp(X[,,j]%*%BETA[j,]),
                       log = TRUE ) ) +
      dmvnorm( t(beta.p),
               theta,
               Sigma,
               log = TRUE) -
      dmvnorm( t(BETA[j,]),
               theta,
               Sigma,
               log = TRUE )
    
    if( log(runif(1)) < lr ) { BETA[j,] <- beta.p }
  }
  ##store some output (only 1 every 10 samples)
  if(s %% 10 == 0){THETA.post <- rbind(THETA.post, t( theta ))}
}
