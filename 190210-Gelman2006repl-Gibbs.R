# Gelman 2006
# Replicated with gibbs sapler

# Package
library(R2OpenBUGS)
library(MCMCpack)
# 1. Data ####
  data("schools")
  J <- nrow(schools)
  nj <- 1 # I couldnt find it, but it's written more than 30 per school
  N <- J*nj
  y <- schools$estimate
  sigma <- schools$sd # known variance within group varaince (different among groups)
  # Parameters definition
  parameters <- c ("theta", "mu", "tau")
  data <- list ("J", "y", "sigma")
  
# 2. Updating functions ####
# > Usuing Priors from the book: Gelman et al 2014 #### 
# For details see book on page 289 and appendix C

  # Full Conditional Posterior of theta.J (or mu.J, the gorup means)
  theta_update <- function (){
      theta_hat <- (mu/tau^2 + y/sigma^2)/(1/tau^2 + 1/sigma^2)
       # vector containing theta_hat for each group
       # Neede: 
       # > y, is the vector of observed group means
       # > sigma, is the vector of group sd (or just sd if equal among all groups)
       # > mu is current value of the grand mean
       # > tau, current value of btw group variance
      V_theta <- 1/(1/tau^2 + 1/sigma^2)
       # Current vector (or single value) of the posterior variance of the group means
      rnorm (J, theta_hat, sqrt(V_theta))
       # Actual sample of the posterior group means for an interation
  }
  
  # Full conditional Posterior of mu (given uniform prior on mu)
    mu_update <- function (){
      rnorm (1, mean(theta), tau/sqrt(J))
    }
    
  # Full conditional Posterior of tau (group level variance)
    tau_update <- function (){
      sqrt(sum((theta-mu)^2)/rchisq(1,J-1))
  }

# > Usuing Priors from: Gelman 2006 #### 
# For details see book on page 289 and appendix C

  # Full Conditional Posterior of theta.J (or mu.J, the gorup means)
  theta_update <- function (){
      theta_hat <- (mu/tau^2 + y*nj/sigma^2)/(1/tau^2 + nj/sigma^2)
       # vector containing theta_hat for each group
       # Neede: 
       # > y, is the vector of observed group means
       # > sigma, is the vector of group sd (or just sd if equal among all groups)
       # > mu is current value of the grand mean
       # > tau, current value of btw group variance
      V_theta <- 1/(1/tau^2 + nj/sigma^2)
       # Current vector (or single value) of the posterior variance of the group means
      rnorm (J, theta_hat, sqrt(V_theta))
       # Actual sample of the posterior group means for an interation
  }
  
  # Full conditional Posterior of mu (given uniform prior on mu)
    # Normal prior on mu
    mu_update <- function (){
      vmu <- 1/(J/tau**2+1/g20**2)
      emu <- vmu*(J*mean(theta)/tau**2 + mu0/g20**2)
      mu <- rnorm(1,
                  emu,
                  sqrt(vmu))
    }
      #sample a new value of mu

  # Full conditional Posterior of tau (group level variance)
    # Uniform or inv-G prior depending on how we specify the parameters
    tau_update <- function (){
      etam <- eta0 + J
      ss   <- eta0*t20 + sum( (theta-mu)^2 )
      sqrt(rinvgamma(1,etam/2,ss/2))
  }
    
# 3. gibss Sampler ####
  # Prior specification
    # Following the indication of Browne and Draper,
    # the U(0, 1/e) and invG(e,e) can be specified using an
    # inverse gamma disitbrution with parameters nu0 and tau02
    # -2 and 0, and 2*e and 1, respectively (where e small number
    # like .001)
    # Prior for tau
    eta0_vec <- c(-2, 2*1, 2*.001)
    t20_vec  <- c(0, 1, 1)
    # Prior for mu
    mu0<- 0 ; g20 <- 1e-6 # match OpenBUGS model
    
priortype <- length(eta0_vec)
iter <- 6000
sims <- array (NA, c(iter, priortype, J+2))
dimnames (sims) <- list (NULL, 
                         c("UNIF", "INV.1.1", "INV.0.0"),
                         c (paste ("theta[", 1:J, "]", sep=""), "mu", "tau"))
for (m in 1:priortype){
  mu   <- rnorm (1, mean(y), sd(y)) # starting value for mu
  tau  <- runif (1, 0, sd(y))       # starting value for tau
  eta0 <- eta0_vec[m]
  t20  <- t20_vec[m]
  for (t in 1:iter){
    theta      <- theta_update ()
    mu         <- mu_update ()
    tau        <- tau_update ()
    sims[t,m,] <- c (theta, mu, tau)
  }
}
# 4. Results ####
  # Point Estimates
    data.frame(mean = colMeans(cbind(sims[-c(1:burnin),1,10],
                                     sims[-c(1:burnin),2,10],
                                     sims[-c(1:burnin),3,10])),
               median = apply((cbind(sims[-c(1:burnin),1,10],
                                     sims[-c(1:burnin),2,10],
                                     sims[-c(1:burnin),3,10])), 2, median)
               )

  # Posterior Plots
    attributes(sims)
    par(mfcol=c(3,1))
    burnin <- 1000
    # Uniform
    hist(sims[,1,10][-c(1:burnin)],
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "tau",
       main = "8 Schools: tau posterior (uniform prior on tau)",
       breaks = 100)
    
    # IG(1, 1)
    hist(sims[,2,10][-c(1:burnin)],
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "tau",
       main = "8 Schools: tau posterior (IG(1,1) prior on tau)",
       breaks = 30)
    
    # IG(.001,.001)
    hist(sims[,3,10][-c(1:burnin)],
       yaxt = "n", ylab = ".", xlim = c(0, 30), xlab = "tau",
       main = "8 Schools: tau posterior (IG(.001,.001) prior on tau)",
       breaks = 100)
    
    