# GRAPHS (posterior distribution for Sigma components, histograms)

# Contains reading functions and graphs for the mix logistic regression 
# approximated models fitted to GPA data and schizophrenia data (+ data desc)

# GPA Data ----------------------------------------------------------------
  # Source: Hox 2010 (p. 92, section 5.2)
  #         [original: unkown]
  # J: 6 (semesters)
  # tot N: 200 individuals
  # Y: GPA (continuous) dichotomized (0/1, above/below sample average)
  # X: Sex (0 = male, 1 = female)
  # T: measurement occasion (0-5 in original, 
  #    centered to -1.46, -0.88, -0.29, 0.29, 0.88, 1.46)

# Load Results
  out1 <- readRDS("./output/GPAdat_matF_R.rds")
  out2 <- readRDS("./output/GPAdat_matF_1e3I.rds")
  out3 <- readRDS("./output/GPAdat_invW_HW.rds")
# N conditions
conds <- c("100", "50", "30", "8", "3") # How many clusters to consider

#> N = 100 ####
which_n <- 1 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 25,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 20),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 25,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 20),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(0, 20),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(0, 10),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 10),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(0, 10),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(-2, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-2, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(-2, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 50 ####
which_n <- 2 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 25,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 20),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 25,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 20),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(0, 20),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(0, 10),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 10),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(0, 10),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(-2, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-2, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(-2, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 30 ####
which_n <- 3 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 50,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 30),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 50,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 30),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 50,
       main = "inv-Wish HW2013",
       xlim = c(0, 30),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 50, 
       main = "mat-F w/ R*",
       xlim = c(0, 15),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 50, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 15),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 50, 
       main = "inv-Wish HW2013",
       xlim = c(0, 15),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 50, 
       main = "mat-F w/ R*",
       xlim = c(-2, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 50, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-2, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 50, 
       main = "inv-Wish HW2013",
       xlim = c(-2, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 8 ####
which_n <- 4 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 50,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 4000),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 150,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 4000),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 250,
       main = "inv-Wish HW2013",
       xlim = c(0, 4000),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 100, 
       main = "mat-F w/ R*",
       xlim = c(0, 500),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 350, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 500),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 150, 
       main = "inv-Wish HW2013",
       xlim = c(0, 500),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 50, 
       main = "mat-F w/ R*",
       xlim = c(-250, 1000),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 500, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-250, 1000),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 250, 
       main = "inv-Wish HW2013",
       xlim = c(-250, 1000),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)
      
  dev.off()
  
#> N = 3 ####
which_n <- 5 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 250,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 100),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 500,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 2000),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 50000,
       main = "inv-Wish HW2013",
       xlim = c(0, 2000),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 1000, 
       main = "mat-F w/ R*",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 10000, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 500),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 100000, 
       main = "inv-Wish HW2013",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 1000, 
       main = "mat-F w/ R*",
       xlim = c(-500, 500),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 1000, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-500, 500),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 3000, 
       main = "inv-Wish HW2013",
       xlim = c(-500, 500),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  

# Schiz 3 Priors --------------------------------------------------------------
  # Source: Hedeker Gibbons 2006, Longitudinal Data Analysis (p.183)
  #         [original: National Institute of Mental Health Schizophrenia Collaborative Study]
  # J = 4 repeated measures (weeks: 0, 1, 3, 6)
  # tot N: 312 individuals
  # Y: Severity of Illness (Schizophrenia),
  #    ordinal (0-7, not_ill-extremely_ill)
  #    dichotomized (< 4 -> 0 ; >= 4 -> 1) (book used only 3 and 4?)
  # X: 0 placebo, 1 drug
  # t: 0, 1, 3, 6 (unit: weeks)
  
# Load Results
  out1 <- readRDS("./output/schizdata_matF_R.rds")
  out2 <- readRDS("./output/schizdata_matF_1e3I.rds")
  out3 <- readRDS("./output/schizdata_invW_HW.rds")
  
# n conditions (randomly sampled; 8 and 4 were hand picked for determination)
conds <- c("100", "50", "30", "20", "8", "4")

#> N = 100 ####
which_n <- 1 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 25,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 300),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 25,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 300),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(0, 300),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(-60, 20),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-60, 20),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(-60, 20),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  
#> N = 50 ####
which_n <- 2 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 25,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 25,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  
#> N = 30 ####
which_n <- 3 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 1000,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 500,
       main = "mat-F w/ 1e3*I",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 20 ####
which_n <- 4 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 500,
       main = "inv-Wish HW2013",
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 200,
       main = "inv-Wish HW2013",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 500,
       main = "inv-Wish HW2013",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 8 ####
which_n <- 5 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 300,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 50,
       main = "mat-F w/ R*",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 400,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 400,
       main = "mat-F w/ 1e3*I",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  
#> N = 4 ####
which_n <- 6 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out1[[which_n]][[9]][,1], breaks = 1000,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 250),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 300,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 5000),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 300,
       main = "inv-Wish HW2013",
       xlim = c(0, 5000),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 500,
       main = "mat-F w/ R*",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 250,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 5000),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 5000,
       main = "inv-Wish HW2013",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 500,
       main = "mat-F w/ R*",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 1000,
       main = "mat-F w/ 1e3*I",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 1000,
       main = "inv-Wish HW2013",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  

# Schiz 4 Priors ----------------------------------------------------------
# Load Results
  out1 <- readRDS("./output/schizdata_matF_R_Ncond50302084.rds")
  out2 <- readRDS("./output/schizdata_matF_1e3I_Ncond50302084.rds")
  out3 <- readRDS("./output/schizdata_invW_HW_Ncond50302084.rds")
  out4 <- readRDS("./output/schizdata_impropSig_Ncond50302084.rds")
# n conditions (randomly sampled; 8 and 4 were hand picked for determination)
conds <- c("50", "30", "20", "8", "4")

#> N = 50 ####
which_n <- 1 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(4,3))
  hist(out1[[which_n]][[9]][,1], breaks = 25,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 25,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,1], breaks = 25,
       main = "uniform",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out4[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 25,
       main = "inv-Wish HW2013",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,4], breaks = 25,
       main = "uniform",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out4[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(0, 20),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)
  
  hist(out1[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ R*",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 25, 
       main = "mat-F w/ 1e3*I",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,2], breaks = 25,
       main = "uniform",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out4[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 25, 
       main = "inv-Wish HW2013",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  
#> N = 30 ####
which_n <- 2 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(4,3))
  hist(out1[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 1000,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,1], breaks = 100,
       main = "uniform",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out4[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 150),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,4], breaks = 100,
       main = "uniform",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out4[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 100),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 500,
       main = "mat-F w/ 1e3*I",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,2], breaks = 500,
       main = "uniform",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out4[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(-15, 15),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 20 ####
which_n <- 3 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(4,3))
  hist(out1[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,1], breaks = 100,
       main = "uniform",
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out4[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 500,
       main = "inv-Wish HW2013",
       xlim = c(0, 600),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,4], breaks = 100,
       main = "uniform",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out4[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 200,
       main = "inv-Wish HW2013",
       xlim = c(0, 300),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ 1e3*I",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,2], breaks = 100,
       main = "uniform",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out4[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 500,
       main = "inv-Wish HW2013",
       xlim = c(-100, 100),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()

#> N = 8 ####
which_n <- 4 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(4,3))
  hist(out1[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 300,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,1], breaks = 300,
       main = "uniform",
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out4[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 500),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 50,
       main = "mat-F w/ R*",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 400,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,4], breaks = 400,
       main = "uniform",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out4[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 100,
       main = "mat-F w/ R*",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 400,
       main = "mat-F w/ 1e3*I",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,2], breaks = 400,
       main = "uniform",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out4[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 100,
       main = "inv-Wish HW2013",
       xlim = c(-200, 10),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()
  
#> N = 4 ####
which_n <- 5 # which n condition
# Posterior Distributions
# With Histograms
  pdf(paste0("./output/schizdata_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(4,3))
  hist(out1[[which_n]][[9]][,1], breaks = 1000,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       xlim = c(0, 250),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,1], breaks = 300,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 5000),
       xlab = "Intercept Variance")
      abline(v = mean(out2[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,1], breaks = 300,
       main = "uniform",
       xlim = c(0, 5000),
       xlab = "Intercept Variance")
      abline(v = mean(out4[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,1]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,1], breaks = 300,
       main = "inv-Wish HW2013",
       xlim = c(0, 5000),
       xlab = "Intercept Variance")
      abline(v = mean(out3[[which_n]][[9]][,1]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,1]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,4], breaks = 500,
       main = "mat-F w/ R*",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,4], breaks = 250,
       main = "mat-F w/ 1e3*I",
       xlim = c(0, 5000),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,4], breaks = 250,
       main = "uniform",
       xlim = c(0, 5000),
       xlab = "Slope Variance")
      abline(v = mean(out4[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,4]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,4], breaks = 5000,
       main = "inv-Wish HW2013",
       xlim = c(0, 250),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 1)

  hist(out1[[which_n]][[9]][,2], breaks = 500,
       main = "mat-F w/ R*",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out2[[which_n]][[9]][,2], breaks = 1000,
       main = "mat-F w/ 1e3*I",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out4[[which_n]][[9]][,2], breaks = 1000,
       main = "uniform",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out4[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out4[[which_n]][[9]][,2]), col = "red", lwd = 1)
  hist(out3[[which_n]][[9]][,2], breaks = 1000,
       main = "inv-Wish HW2013",
       xlim = c(-1000, 200),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 1)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 1)

  dev.off()