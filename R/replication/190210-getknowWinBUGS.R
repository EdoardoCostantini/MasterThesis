# Understanding WinBUGs

# Trying to replicate what was done in Gelman2006
# while learning how to use WinBUGS in R.
# Important: Bugs requires proper prior distirbutions.

# Sources: Gelman2006, ManualR2WinBUGS.pdf
# Manual for OpenBUGS: https://oliviergimenez.github.io/post/run_openbugs_on_mac/

#Packages
library(R2OpenBUGS)

data(schools)
schools

#Data
y       <- schools$estimate # dependent variable
J       <- nrow(schools)    # number of clusters
sigma.y <- schools$sd       # observed within gorups sd
data    <- list ("J", "y", "sigma.y")

# Function to specify the initial values
inits <- function(){
  list(theta = rnorm(J, 0, 100),
       mu.theta = rnorm(1, 0, 100), 
       sigma.theta = runif(1, 0, 100))
}

# Set the Wine working directory and the directory to OpenBUGS,
# and change the OpenBUGS.exe location as necessary:
WINE = "/usr/local/bin/wine"
WINEPATH = "/usr/local/bin/winepath"
OpenBUGS.pgm = "/Applications/OpenBUGS323/OpenBUGS.exe"

# Define Parameters
parameters = c("theta", "mu.theta", "sigma.theta")

# Define Model File Location
model.file1 <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/BugsModels/ManualR2OpenBUGS_model1.txt"

# Start Simulation
schools.sim <- bugs(data, inits, model.file = model.file1,
                    parameters = parameters,
                    n.chains = 3, n.iter = 1000, 
                    OpenBUGS.pgm = OpenBUGS.pgm, WINE = WINE, WINEPATH = WINEPATH, useWINE = T) # very imp on mac
  # ERRORS: - You might get a weird message starting by err:ole, just ignore it.
  #         - It seems that you have to run it twice to get it correct, why?

# Results
  print(schools.sim)
  plot(schools.sim)
# Access posterior simulations in R
  attributes(schools.sim$sims.list)
# Plot posterior
  par(mfcol=c(3,1))
  hist(schools.sim$sims.list$sigma.theta,
       yaxt = "n", ylab = ".",
       xlim = c(0, 30),
       main = "8 Schools: posterior on sigma_a given uniform prior on sigma_a",
       breaks = 20)
# Posterior predictive simulation
  theta <- schools.sim$sims.list$theta
  y.rep <- array (NA, c(1000, J))
  for (sim in 1:1000)
    y.rep[sim,] <- rnorm(J, theta[sim, ], sigma.y)

