# Title:       Replicate findings in Haung Wand 2013
# Description: Trying to replicate to understand the prior for Sigma proposed by them
#              thoguh replciation of the analysis of the real longitudinal data.

# Ste up
library(mvtnorm)
library(MCMCpack)

# Data
  pigsdata <- read.table("./Data/pigsdata.txt")
  # The weight of 16 pigs is recorded (48 in original data)  
  # every week for 9 weeks.
  
# Hierarchical Model
  
