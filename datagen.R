# Title:       Data generation script
# Description: Generation of data based on random intercept/slope model
#              default: repeated measures in individuals ####

library(lme4)
library(lattice)
library(mvtnorm)

# values are based on HedekerGibbons2006 chapter 4 data
# to get a sense

# MODEL:
# LVL1: y_ij = b0_i + b1_i * t_ij + eij   where eij ~ N(0, sigma2e)
# LVL2: b0_i = gamma00 + u0_i             where u_vec ~ multiv-N(0_vec, SIGMA)
#       b1_i = gamma10 + u1_i             w/ SIGMA has: - diagonal sigma2u0, sigma2u1; 
#                                                       - off-diagonal sigmau01
set.seed(201902)
I <- 65 # number of individuals original: 65
nj <- 6 # number of students (start with balanced desing/no missings)
N <- nj*I # number of observations, original: 6*65

# True values
  sigma2e  <- 12    # between observations variance
  sigma2u0 <- 12    # between individuals (intercepts) variance
  sigma2u1 <- 2     # between individuals (slopes) variance
  sigmau01 <- -1.5  # random slope intercept covariance
                    # positive -> higer initial values = greater positve slopes
  SIGMA <- matrix(c(sigma2u0,sigmau01,
                    sigmau01,sigma2u1),
                  nrow = 2, byrow = T,
                  dimnames = list(c("u0", "u1"), c("u0", "u1")))
  
  gamma00 <- 24 # grand mean fixed
  gamma10 <- -3 # fixed effect of time
  
  eij <- rnorm(N, 0, sqrt(sigma2e))
  uj  <- rmvnorm(I,
                 mean = c(0,0),
                 sigma = SIGMA)
    colnames(uj) <- c("u0", "u1")

# Random effects for data
  uj_data <- cbind(rep(uj[,1], each = nj), rep(uj[,2], each = nj))

# Predictors
  tij <- rep(seq(0, nj-1), I) # time (this inerval should cover all
                           # the repeated observations for each individual)

# Data
  yij <- gamma00 + uj_data[,1] + gamma10 * tij + uj_data[,2]*tij + eij
  ind_id <- rep(seq(1:I), each = nj)
  dataset <- data.frame(ind_id = ind_id, 
                        yij = yij, 
                        tij = tij)
# Simple OLS
  summary(lm(yij ~ tij))
# Random intercept and slope model
  head(dataset)
  RS_model <- lmer(yij ~ tij + (1 + tij | ind_id), data = dataset, REML = FALSE)
  summary(RS_model)
  as.data.frame(VarCorr(RS_model))
  # everything very close to what you expected!
  # To do:
  # - make it a function and think about how
  # - make it flexible for more than 1 predictors
  
# Plot
  # OLS trajectories superimposed on the empirical growth plots
  xyplot(yij ~ tij | ind_id, data = dataset, 
         panel = function(x, y){
           panel.xyplot(x, y)
           panel.lmline(x, y)
         }, as.table=TRUE)

  # Plot each line
  predscore <- fitted(RS_model)
  datapred <- data.frame(predscore = predscore,
                    tij = dataset$tij,
                    ind_id = dataset$ind_id)
  xyplot(predscore ~ tij, data = datapred,
         groups = ind_id, 
         type = c("l"), col = "brown"
         # panel = function(tij, predscore,...) {
         #   panel.xyplot(tij, predscore,...)
         #   panel.text(1,65,labels=letters[panel.number()])
         #   })
         # panel.groups = function(tij, predscore, ind_id){
         #   panel.xyplot(tij, predscore, t="b",...) 
         #   panel.grid(v=-1, h=-1, lty=3)
         #   xt <- tij[tij == max(tij)]
         #   yt <- predscore[tij == max(tij)]
         #   panel.text(xt, yt, labels = unique(ind_id))
         # }
         )
  # Selected individuals plot
  ind_plot1 <- datapred[, 3] == 4
  ind_plot2 <- datapred[, 3] == 24
  ind_plot <- datapred[, 3] == 4 | datapred[, 3] == 24 # which individuals do you want to plot
  dataplot <- datapred[ind_plot, ]
  selec_plot <- xyplot(predscore ~ tij, data = dataplot,
                       groups = ind_id, 
                       type = c("l"), col = "brown")
  selec_obs1 <- latticeExtra::layer(panel.points(0:5, # fix every time the upper bound!
                                                 dataset[ind_plot1, 2], pch = 1, col = "blue"))
  selec_obs2 <- latticeExtra::layer(panel.points(0:5, 
                                                 dataset[ind_plot2, 2], pch = 1, col = "red"))
  selec_plot + selec_obs1 + selec_obs2

# Predictors
  # if you wanted to add a predictor you could do something like this
EXTRmeans <- runif(J, 8, 12) 
EXTR <- vector("list", J)
for (j in 1:J) {
  EXTR[[j]] <- rnorm(nj,
                   EXTRmeans[j],
                   1)
}
mean(EXTR[[1]])



# SCRAP ####

# Students in School
# Sample
J <- 48 # number of schools
N <- 4800 # number of students
nj <- N/J # number of students (start with balanced desing)

#ture values
sigma2e  <- 10 # between students variance
sigma2u0 <- 40 # between schools (intercepts) vairance
sigma2u1 <- 15 # between schools (slopes) vairance
sigmau01 <- 2  # random slope intercept covairance
               # positive -> higer initial values = higergreater positve slopes
SIGMA <- matrix(c(sigma2u0,sigmau01,
                  sigmau01,sigma2u1), nrow = 2, byrow = T)

gamma00 <- 30 # grand mean fixed
gamma01 <- 10 # fixed effect of teacher experience
gamma10 <- .5 # fixed effect of pupil's sex
gamma20 <- 5 # fixed effect of pupil extroversion

eij <- rnorm(N, 0, sigma2e)
uj <- rmvnorm(J,
              mean = c(0,0),
              sigma = SIGMA)

# Predictors
EXTRmeans <- runif(J, 8, 12) 
EXTR <- vector("list", J)
for (j in 1:J) {
  EXTR[[j]] <- rnorm(nj,
                   EXTRmeans[j],
                   1)
}
mean(EXTR[[1]])

group_id <- rep(seq(1:J), each = nj)
colMeans(do.call(cbind, EXTR))
