# Analysis of real data with Logistic Model

# Contains script to fit mixed logistic regression approximated
# for longitudianl data with different priors (3 at the moment)

# Load Functions
source("./190304-log_reg_mixed-functions.R")
library(ddpcr) # for quiet function

# Store results
  conds <- c("50", "30", "20", "8", "4") # How many clusters to consider
  out1 <- out2 <- out3 <- out4 <- vector("list", length = length(conds))
  names(out1) <- names(out2) <- names(out3) <- names(out4) <- conds
  set.seed(20190308)
  go <- Sys.time()
  
# Perform MCMC ####
for (outps in 1:length(conds)) {
  
  print(paste("n condition:", conds[outps]))
  # Define Data
  dat <- data.frame(cluster = schizdata.noNA$id, #= gpadata$student,
                    yvec    = schizdata.noNA$SevIll,     #= gpadata$gpa,
                    tvec    = schizdata.noNA$week,  #= gpadata$occas,
                    xvec    = schizdata.noNA$drug)    #= gpadata$sex)
  # Reduce number of clusters
  n_goal <- as.numeric(conds[outps])                                         # how many clusters do you want to work with?
  clusters_goal <- sample(unique(dat$cluster), n_goal) # sample that many from full dataset
  # For 8 and 4 cases we need special attention:
  # If by random choice all elements in the covariate or outcome are the same, then
  # you can't go on. Therefore we select specific cases instaed of oging at random
  # For the 8 cases scenario we select these ids: 1129, 1306, 2302, 3106 (no drug), 1103, 1113, 1109, 2316 (drug)
  # For the 4 cases scenario we select these ids: 1129, 1306 (no drug), 2316, 1103 (drug)
  if(n_goal == 8){
    clusters_goal <- c(1129, 1306, 2302, 3106, 1103, 1113, 1109, 2316)
  }
  if(n_goal == 4){
    clusters_goal <- c(1129, 1306, 2316, 1103)
  }
  dat <- dat[dat$cluster %in% clusters_goal, ]         # keep only those clusters
  
  # Prepare Inputs for functions
  yvec     <- dat$yvec
  Xmat     <- cbind(rep(1, nrow(dat)), dat$tvec, dat$xvec, dat$tvec*dat$xvec)
  # dat      <- list(yvec, Xmat)
  # 
  # yvec <- dat[[1]]
  # Xmat <- dat[[2]]
  J    <- length(unique(Xmat[,2]))
  n    <- length(yvec)/J # must be even!
  
  # Define Priors
  nu1=7.3
  s2tilde=pi**2*(nu1-2)/(3*nu1)
  
  # Define iterations
  samsize = 1e4
  
  # Perform MCMC
  print("> prior 1: Matrix-F prior with R*")
  quiet(
    out1[[outps]] <- MCMC_mixlogreg_SBeta2_A(yvec,Xmat,n,J,
                                             B0 = 1,
                                             samsize=samsize)                        # Matrix-F prior with R*
  )
  print("> prior 2: Neighbour of Sigma**-(1/2)")
  quiet(
    out2[[outps]] <- MCMC_mixlogreg_SBeta2_A(yvec,Xmat,n,J,samsize=samsize,
                                             B0 = 1e3*diag(2))                       # Neighbour of Sigma**-(1/2)
  )
  print("> prior 3: Huang Wand")
  quiet(
    out3[[outps]] <- MCMC_mixlogreg_HW_A(yvec, Xmat, n, J, samsize=samsize) # Haung Wand Prior
  )
  print("> prior 4: Improper solve(sigma)")
  quiet(
    out4[[outps]] <- MCMC_mixlogreg_JB_A(yvec, Xmat, n, J, samsize=samsize) # Improper inverse
  )
}
  
stop <- Sys.time()
stop - go

# Opitonal
out4 <- MCMC_mixlogreg_SBeta2_A(yvec,Xmat,n,J,
                                B0 = 0,
                                samsize=samsize) # Matrix-F prior with B_hat

# Save Resulsts
  saveRDS(out1, paste0("./output/schizdata_", "matF_R", "_Ncond", paste(conds,collapse=""), ".rds"))
  saveRDS(out2, paste0("./output/schizdata_", "matF_1e3I", "_Ncond", paste(conds,collapse=""), ".rds"))
  saveRDS(out3, paste0("./output/schizdata_", "invW_HW", "_Ncond", paste(conds,collapse=""), ".rds"))
  saveRDS(out4, paste0("./output/schizdata_", "impropSig", "_Ncond", paste(conds,collapse=""), ".rds"))
# Load Results
  out1 <- readRDS("./output/GPAdat_matF_R.rds")
  out2 <- readRDS("./output/GPAdat_matF_1e3I.rds")
  out3 <- readRDS("./output/GPAdat_invW_HW.rds")

# Compare Results ####
which_n <- 3 # which n condition
# For GPA data
conds <- c("100", "50", "30", "8", "3") # How many clusters to consider
# Trace Plots
plot(1:samsize,out1[[which_n]][[9]][,1],"l")
  lines(1:samsize,out2[[which_n]][[9]][,1],col=2)
  lines(1:samsize,out3[[which_n]][[9]][,1],col=2)

# Posterior Distributions
# With Histograms
  pdf(paste0("./output/GPAdat_", "n", conds[which_n], ".pdf"))
  par(mfcol = c(3,3))
  hist(out4[[which_n]][[9]][,1], breaks = 100,
       main = paste("n =", conds[which_n], ", mat-F w/ R*"),
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Intercept Variance")
      abline(v = mean(out1[[which_n]][[9]][,1]), col = "blue", lwd = 2)
      abline(v = median(out1[[which_n]][[9]][,1]), col = "red", lwd = 2)
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
  
  hist(out1[[which_n]][[9]][,4], breaks = 100, 
       main = "mat-F w/ R*",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Slope Variance")
      abline(v = mean(out1[[which_n]][[9]][,4]), col = "blue", lwd = 2)
      abline(v = median(out1[[which_n]][[9]][,4]), col = "red", lwd = 2)
  hist(out2[[which_n]][[9]][,4], breaks = 100, 
       main = "mat-F w/ 1e3*I",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Slope Variance")
      abline(v = mean(out2[[which_n]][[9]][,4]), col = "blue", lwd = 2)
      abline(v = median(out2[[which_n]][[9]][,4]), col = "red", lwd = 2)
  hist(out3[[which_n]][[9]][,4], breaks = 100, 
       main = "inv-Wish HW2013",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Slope Variance")
      abline(v = mean(out3[[which_n]][[9]][,4]), col = "blue", lwd = 2)
      abline(v = median(out3[[which_n]][[9]][,4]), col = "red", lwd = 2)
  
  hist(out1[[which_n]][[9]][,2], breaks = 100, 
       main = "mat-F w/ R*",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out1[[which_n]][[9]][,2]), col = "blue", lwd = 2)
      abline(v = median(out1[[which_n]][[9]][,2]), col = "red", lwd = 2)
  hist(out2[[which_n]][[9]][,2], breaks = 100, 
       main = "mat-F w/ 1e3*I",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out2[[which_n]][[9]][,2]), col = "blue", lwd = 2)
      abline(v = median(out2[[which_n]][[9]][,2]), col = "red", lwd = 2)
  hist(out3[[which_n]][[9]][,2], breaks = 100, 
       main = "inv-Wish HW2013",
       #xlim = c(0, 20), ylim = c(0, 1500),
       xlab = "Slope-Int Covariance")
      abline(v = mean(out3[[which_n]][[9]][,2]), col = "blue", lwd = 2)
      abline(v = median(out3[[which_n]][[9]][,2]), col = "red", lwd = 2)
  
  dev.off()
  
  median(out1[[which_n]][[9]][,2])
  
# With density plots
  jpeg("./output/hist_gpa_sat.jpg")
  
  plot(density(out1[[which_n]][[9]][,1]),xlim=c(0,20),ylim=c(0,.5), lty = 1,
       main = "Posterior Density Intercepts Variance",
       sub = "Matrix F R*")
  #par(new=T)
  lines(density(out2[[which_n]][[9]][,1]),xlim=c(0,20),ylim=c(0,.5), lty = 2)
  lines(density(out3[[which_n]][[9]][,1]),xlim=c(0,20),ylim=c(0,.5), lty = 3)
  
  plot(density(out1[[which_n]][[9]][,4]),xlim=c(0,8),ylim=c(0,2), lty = 1)
  par(new=T)
  plot(density(out2[[which_n]][[9]][,4]),xlim=c(0,8),ylim=c(0,2), lty = 2)
  par(new=T)
  plot(density(out3[[which_n]][[9]][,4]),xlim=c(0,8),ylim=c(0,2),col=2, lty = 3)
  par(new=T)

  plot(density(out1[[which_n]][[1]][[9]][,2]),ylim=c(0,2), lty = 1)
  par(new=T)
  plot(density(out2[[which_n]][[9]][,2]),ylim=c(0,2), lty = 2)
  par(new=T)
  plot(density(out3[[which_n]][[9]][,2]),ylim=c(0,2), lty = 3)
  par(new=T)