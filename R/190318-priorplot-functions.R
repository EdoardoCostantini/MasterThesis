### Project:     Master Thesis
### Object:      Function to compute density for univ-F distribution with free scale parameter choice.
### Description: Contains the R-function for:
###              1) density for univariate-F (free scale parameter) defined in Mulder Pericchi 2018 in eq 1 and 2
###              2) density for noncentral-t
###              3) density for folded-t
###              You can call this functions to plot the priors for the parameters of interest (variances)
### Date:        2019-03-19


# R-functions for univariate-F --------------------------------------------
  
# F disitrbution with free scale parameter (b) for variance
  df_freeb <- function(x, nu, d, b){
    ( gamma( (d+nu)/2 ) / ( gamma(nu/2)*gamma(d/2)*b**(nu/2) ) ) * (x)**((nu/2) - 1) * (1 + (x)/b)**(-(nu+d)/2)
  }
# Resulting disitribution of sd (from F dist w/ free scale parameter (b) on var)
  df_freeb_SD <- function(x, nu, d, b){
    ( 2*gamma( (d+nu)/2 ) / ( gamma(nu/2)*gamma(d/2)*b**(nu/2) ) ) * (x)**(nu - 1) * (1 + (x**2)/b)**(-(nu+d)/2)
  }
  
# Non centralised t density function example
# Source: Defining an arbitrary density function (Bayesian Statisitcs Course: exercise set 6.3c) 
# pdf: https://en.wikipedia.org/wiki/Student%27s_t-distribution#Generalized_Student's_t-distribution
  #  logarithmic version
  dt_nc_sigma2_log <- function(x,nu,mu,sigma2){
    exp( lgamma( (nu+1)/2 ) - lgamma(nu/2) - .5*log(pi*nu*sigma2) - (nu+1)/2 * log( 1+(x-mu)**2/(nu*sigma2) ) )
  }
  
  dt_nc_sigma <- function(x, nu, mu, sigma){
    (gamma((nu + 1)/2))/(gamma(nu/2)*sqrt(pi*nu)*sigma) * (1+1/nu*((x-mu)**2)/sigma)**(-(nu+1)/2)
  }
  
  # # How does it look like?
  # sdseq <- sqrt(seq(0, 3000, length = 100000))
  # plot(sdseq, dt_nc_sigma2_log(sdseq, nu = 2, mu = 0, sigma2 = 100), type = "l")
  # plot(sdseq, dt_nc_sigma(sdseq, nu = 2, mu = 0, sigma = 100), type = "l")
  
# Folded-t (Half-Cauchy when mu = 0, nu = 1)
  # pdf: https://en.wikipedia.org/wiki/Student%27s_t-distribution#Generalized_Student's_t-distribution
  # w/ taking absolute value of x (folded), setting mu = 0 (half-t)
  # In terms of sigma**2
  dt_folded_sigma2 <- function(x,nu,mu,sigma2){
    exp( lgamma( (nu+1)/2 ) - lgamma(nu/2) - .5*log(pi*nu*sigma2) - (nu+1)/2 * log( 1+(x-mu)**2/(nu*sigma2) ) ) + 
      exp( lgamma( (nu+1)/2 ) - lgamma(nu/2) - .5*log(pi*nu*sigma2) - (nu+1)/2 * log( 1+(x+mu)**2/(nu*sigma2) ) )
  }
  # In terms of sigma
  dt_folded_sigma <- function(x, nu, mu, sigma){
    (gamma((nu + 1)/2))/(gamma(nu/2)*sqrt(pi*nu)*sigma) * (1+1/nu*((x-mu)**2)/sigma)**(-(nu+1)/2) + 
      (gamma((nu + 1)/2))/(gamma(nu/2)*sqrt(pi*nu)*sigma) * (1+1/nu*((x+mu)**2)/sigma)**(-(nu+1)/2)
  }
  
  # # How deos it look like
  # sdseq <- sqrt(seq(0, 3000, length = 100000)) # important to plot the piror
  # plot(sdseq, dt_folded_sigma2(sdseq, nu = 2, mu = 0, sigma2 = 100), type = "l")
  # plot(sdseq, dt_folded_sigma(sdseq, nu = 2, mu = 0, sigma = 100), type = "l")
  
# Inverse Wishart
  draw_InvWish = function(n,bMat,S0=diag(2)){
  ScaleMatrix = S0
  PsiInvDraw = rwish(v = n + 1,               # nu0 = 1
                                              # usually nu0 = 2, v = n + 2; you changed to 1 for reasons specified in the manuscript
                     S = solve(ScaleMatrix))  # distirbution is just Wishart, not inverse! Therefore, the guess priovaded is invterted!
  return(PsiInvDraw)
}

  
  
# Scrap  

  # x <- abs(seq(-30, 30))
  # plot(x, 
  #      dt_nc(x,nu = 1, mu = 0, sigma2 = 1), # these values of nu and mu make this dist. an half cauchy
  #      type = "l", ylim = c(0, 1))
  # plot(x,
  #      dt_folded(x,nu = 1, mu = 0, sigma2 = 1), # these values of nu and mu make this dist. an half cauchy
  #      type = "l", ylim = c(0, 1))
  # #Large but ﬁnite values of sigma2 represent prior distributions which we call “weakly informative” because, 
  # #even in the tail, they have a gentle slope

# # Try the functions out to plot priors
#   varseq <- seq(0, 100, length = 1000)
#   plot(varseq, df_freeb(varseq, nu = 2, d = 1, b = 20), type = "l")
#   
#   sdseq <- sqrt(varseq)
#   plot(sdseq, 
#        df_freeb_SD(sdseq, nu = 2, d = 1, b = sqrt(20)), type = "l")
#   
#   # Standard F
#   plot(x, 
#        df_freeb(x, nu = 2, d = 1, b = 1), type = "l",
#        ylim = c(0, 1))
#   lines(x, df(x, 2, 1), type = "l", # standard F, df_freeb should give same plot when b = 1
#        ylim = c(0, 1), col = 3)
#   plot(x, 
#        df_freeb_SD(x, nu = 2, d = 1, b = 1), type = "l",
#        ylim = c(0, 1))
#   lines(x, df(x, 2, 1), type = "l", # standard F, df_freeb should give same plot when b = 1
#        ylim = c(0, 1), col = 3)
# # Plotting Priors
# # (from Bayesian Statistics Course ex5_1b)
#   # Data
#     n <- 22
#     x <- c(3,0,0,2,0,1,1,2,4,2,0,7,2,0,1,3,3,1,2,2,0,0)
#     sigmax <- sum(x)
#   # Prior: gamma(a, b)
#     a <- 2
#     b <- 1
#   # Plots
#     lambdaseq <- seq(0,10,length=1e4)
#     # Prior distribution
#       # Actual
#         plot(lambdaseq, dgamma(lambdaseq, a, b), type = "l", col = 2)
#       # Simulated
#         numdraws <- 1e7
#         draws <- rgamma(numdraws, shape = a, rate = b)
#         plot(density(draws), xlim = c(0, 15), ylim = c(0, 1.5), col = 2)
#     # Legend
#       legend(4, 1.5,
#              legend = c("prior", "norm. likelihood", "posterior"),
#              lwd = rep(1, 3), col = c("red", "blue", "green"))
#     # Postirior distribution
#       # Actual
#         # plot(lambdaseq,
#         #       dgamma(lambdaseq, shape = a + sigmax, rate = b + n),
#         #       lty = 1, col = "green")
#       # Simulate (monte carlo estimate)
#         postdraws <- rgamma(numdraws,
#                         shape = a + sigmax,
#                         rate = b + n)
#         #plot(density(postdraws), xlim = c(0, 15), col = "red")
#         lines(density(postdraws), xlim = c(0, 15), col = "green")
#     # Likelihood
#       par(new=T)
#       plot(lambdaseq,
#            exp( - n * lambdaseq ) * lambdaseq**sum(x), #poissan likelihood function (kernel)
#            type = "l", col = 4, xaxt = "n", yaxt = "n", ylab = "", xlab = "", xlim = c(0, 15))
  
