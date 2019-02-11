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
<<<<<<< HEAD
model.file1 <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/BugsModels/ManualR2OpenBUGS_model1.txt"
=======
model.file1 <- "/Users/Edoardo/DriveUni/MasterThesis/BayesianModeling/ManualR2OpenBUGS_model1.txt"
>>>>>>> 76f3faa1e28b9862b98956cf2b0407cd28eb0e6f

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

