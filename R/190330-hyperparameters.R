### Project:     Master Thesis
### Object:      Prior Definition
### Description: Contains the code I used to define the prior guesses. You can source this code to get the lists cointaining
###              the priors of interest in your environment. The call to the data is kept in the comments in case one wants
###              to run this code independently. However, the data will usually already be in the environment.
### Date:        2019-03-30

# Normal model:
# y_{ij} = normal(p_{ij})
# y_{ij} = beta0 + beta1*t_j + beta2*x_i + beta_3*t_j*x_i + b_{i0} + b_{i1}*t_j
# (b_{i0},b_{i1})' ~ N(c(0,0),Psi)

# # Riesbt dat
#   RiesbyDat <- read.table("./data/RiesbyDat.txt")
#   dat_ALL <- data.frame(cluster = RiesbyDat$id,    #constant dataset w/ naming that fits the loop
#                         yvec    = RiesbyDat$depr,  #you only need to change the variables included here
#                         xvec    = RiesbyDat$week,  #and everything will be adapted in the loop
#                         cvec    = RiesbyDat$endog,
#                         inter   = RiesbyDat$inter)
#   dat_Zi   <- cbind(rep(1, length = length(unique(dat_ALL$xvec))), unique(dat_ALL$xvec))
#   dat_ALLn <- nrow(dat_ALL)/length(unique(dat_ALL$xvec)) # need it for conditons definition

  dimens <- ncol(dat_Zi)   # k
  epsilon <- c(1, 1/2, .1) # choose the e values to plug in "nu = k - 1 + e" (and "d = e" when needed)
  
# Educated Guess ----------------------------------------------------------
  # Define the educated guess by fitting a regression line for each group and getting
  # the variance of the parameters (intercept and slope)
  intercept <- rep(NA, dat_ALLn)
  slope <- rep(NA, dat_ALLn)
  i <- 1
  for (id in 1:dat_ALLn) {
    #print(paste("Before:", i))
    fit <- lm(dat_ALL$yvec[i:(i+5)] ~ scale(dat_ALL$xvec[i:(i+5)]))
    intercept[[id]] <- coef(fit)[1]
    slope[[id]]     <- coef(fit)[2]
    
    i <- i+6
    #print(paste("After:", i))
  }
  vi <- var(intercept) # 20 # these guesses are based on the entire sample
  vs <- var(slope)     # 9
    
  B0_ed  <- matrix(c(vi, 0, 0, vs), ncol = 2) # guess based on data exploration and knowledge
  
# Empirical Bayes ---------------------------------------------------------
  # mock version: need to work on this definition
  Rstar <- solve(t(dat_Zi)%*%dat_Zi) #this will be defined in the loop because it depedns on Zi
  
# Proper Neighbor ---------------------------------------------------------
  # see Mulder Pericchi 2018 (or equivalent) for proper neightbur of |Psi|^(-1/2) prior
  B0_pn  <- 1e3*diag(2) # mat-f proper neighbour guess
  
# Inverse Wishart Priors --------------------------------------------------
# IW uninformative
  # Start by creating the priors with the k - 1 + e scheeme for different values of e
  IW_PR <- sapply(epsilon, 
                    FUN = function(k, e){
                      nu = k - 1 + e
                      B0_IWunin = (k - 1 + e)*diag(2)
                      return(list(e = e, S0 = B0_IWunin))},
                    k = dimens, simplify = F)
  # Then add the educated prior guess (computed above) with a given e
  IW_PR[[length(IW_PR)+1]] <- list(e=1,S0=B0_ed)
  # Name all the priors
  names(IW_PR) <- c("IW_e2","IW_e05","IW_e01","IW_eg")
  # Check the structure of the list
  str(IW_PR)


# Matrix F ----------------------------------------------------------------
  # Create the different versions of the uninformative prior 5
  MF_PR <- sapply(epsilon, 
                  FUN = function(k, e){
                    nu = k - 1 + e
                    d = e
                    B0_IWunin = B0_ed
                    return(list(nu = nu, d = e, e = e, S0 = B0_IWunin))},
                  k = dimens, simplify = F)
  # Add improper prior, educated guess informative, R* uninformative?
  MF_PR[[length(MF_PR)+1]] <- list(nu=2,d=2,e=0,S0=10**3*diag(2))
  MF_PR[[length(MF_PR)+1]] <- list(nu=2,d=2,e=0,S0=B0_ed) # this adds the educated prior guess
  MF_PR[[length(MF_PR)+1]] <- list(nu=2,d=2,e=0,S0=Rstar) # this adds the educated prior guess
  # Provide names
  names(MF_PR) <- c("MF_e1","MF_e05","MF_e01","MF_pn","MF_eg","MF_R*")
  # Check structure
  str(MF_PR)
