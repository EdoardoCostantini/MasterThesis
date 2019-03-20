### Project:     Master Thesis
### Object:      Loading and exploring data
### Description: Contains the loading functions and exploration of datasets found online and in books.
###              The goal of the script is to find a dataset suitable for the purposes of the project
### Date:        2019-03

# Loading datasets found online to see which one could be used
  library(tidyverse)
  library(haven)
  
  # Data Restrictions:
  # - Balanced Desing
  # - No missing values
  # - Social Sciences
  # - Repeated measures
  # - Time constant covariate
  
# NIMH Schizophrenia ####
  # Used by Hedeker Gibbons 2006, Longitudinal Data Analysis
  # Ordinal Outcome, dichotomized
  schizdata <- as.data.frame(read.table("./data/SCHIZX1.DAT.txt")[-c(3060),])
  colnames(schizdata) <- c("id", "SevIll", ".", ".", ".", "drug", "week")
  
  # Study Dataset
  head(schizdata, 12)
  summary(schizdata[,6][which(schizdata[,2]!=c(-9))])
  length(unique(schizdata[which(schizdata[,6] == 1), 1])) # how many in place/drug conditon
  
  store <- NA # correct number of observations per id
  for (i in 1:length(unique(schizdata$id))) {
    store[i] <- sum(schizdata$id == unique(schizdata$id)[i])
  }
  sum(schizdata[schizdata[, 1] == unique(schizdata[which(schizdata[,6] == 1), 1]), ])
  
  # Select Columns
  schizdata <- schizdata[, c(1, 2, 6, 7)]
  
  # Missings
  schizdata$SevIll[which(schizdata$SevIll == c(-9))] <- NA
  
  # Recode Time
  schizdata$week <- rep(seq(-3, 3, 1), nrow(schizdata)/7)
  
  #Dichotomize y
  schizdata$SevIll[which(schizdata$SevIll < 4)] <- 0
  schizdata$SevIll[which(schizdata$SevIll >= 4)] <- 1
  
  # Create list for functions
  yvec <- schizdata$SevIll[-c(1:7)]
  Xmat <- cbind(rep(1, nrow(schizdata)), schizdata$week, schizdata$drug, schizdata$drug*schizdata$week)[-c(1:7), ]
  dat <- list(yvec, Xmat)
  
  #Convert to wide format (no need)
  outdata <- data.frame(id = schizdata$id,
                        y = schizdata$SevIll,
                        week = schizdata$week)
  which(schizdata$id == "\032")
  which(store != 7)
  schizdata_wide <- spread(outdata, week, y)
  # Missings
  for (i in 2:ncol(schizdata_wide)) {
    schizdata_wide[, i][which(schizdata_wide[, i] == -9)] <- NA
  }
  colSums(is.na(schizdata_wide) == FALSE) 
  # Add treatment type (0 = placebo, 1 = drugs)
  schizdata_wide$drug1 <- schizdata[, 6][seq(1, nrow(schizdata), by = 7)]
  
  head(schizdata_wide) # you can use this data. but missings are problem.
  
  #> No missing version ####
  # Source: Hedeker Gibbons 2006, Longitudinal Data Analysis (p.183)
  #         [original: National Institute of Mental Health Schizophrenia Collaborative Study]
  # J = 4 repeated measures (weeks: 0, 1, 3, 6)
  # tot N: 312 individuals (final)
  # Y: Severity of Illness (Schizophrenia),
  #    ordinal (0-7, not_ill-extremely_ill)
  #    dichotomized (< 4 -> 0 ; >= 4 -> 1) (book used only 3 and 4?)
  # X: 0 placebo, 1 drug
  # t: 0, 1, 3, 6 weeks (final)
  schizdata.noNA <- as.data.frame(read.table("./Data/SCHIZREP.DAT_nomiss.txt"))
  head(schizdata, 12)
  head(schizdata.noNA, 12)
  
  colnames(schizdata.noNA) <- c("id", "SevIll", "week", "drug")
  
  #get rid of 4 and 5 in week
  schizdata.noNA$week[which(schizdata.noNA$week == 2)] <- NA
  schizdata.noNA$week[which(schizdata.noNA$week == 4)] <- NA
  schizdata.noNA$week[which(schizdata.noNA$week == 5)] <- NA
  schizdata.noNA <- na.omit(schizdata.noNA)
  
  
  index <- NULL #index those that have 4 observations (0, 1, 3, 6)
  for (i in 1:nrow(schizdata.noNA)) {
    if(sum(unique(schizdata.noNA$id)[i] == schizdata.noNA$id) == 4){
      index[i] <- unique(schizdata.noNA$id)[i]
    }
  }

  new_index <- na.omit(index) # leave out NAs
  schizdata.noNA <- schizdata.noNA[schizdata.noNA$id %in% new_index, ] #select cases with 0, 1, 3, 6 observations
  
  unique(schizdata.noNA$week)
  unique(schizdata.noNA$id)
  
  #Dichotomize y
  # schizdata.noNA$SevIll[which(schizdata.noNA$SevIll < 4)] <- 0
  # schizdata.noNA$SevIll[which(schizdata.noNA$SevIll >= 4)] <- 1
  
  schizdata.noNA <- schizdata.noNA[, c(1:4)]
  schizdata.noNA$inter <- schizdata.noNA$week*schizdata.noNA$drug
  # For 8 and 3 cases you need special attention:
  # If by random choice all elements in the covariate are the same, then
  # You can't go on. Therefore we select specific cases instaed of oging at random
  # For the 8 case scenario we select these ids:
  
  schizdata.noNA[order(schizdata.noNA$drug, decreasing = TRUE), ]
  nrow(schizdata.noNA)/4
  head(schizdata.noNA)
  write.table(schizdata.noNA, "./data/schizdata_noNA.txt")
  read.table("./data/schizdata_noNA.txt")
  
# Riesby data (depression scale) ####
  # Used by Hedeker Gibbons 2006, Longitudinal Data Analysis (p.68)
  # Can look at slides for more info on data
  RiesbyDat <- read.table("./data/RIESBY.DAT.txt")[, -3]
  colnames(RiesbyDat) <- c("id", "depr", "week", "endog", "inter")
  head(RiesbyDat)
  which(is.na(RiesbyDat))
  RiesbyDat <- RiesbyDat[-397, ]
  RiesbyDat <- apply(RiesbyDat, 2, as.numeric)
  RiesbyDat <- as.data.frame(RiesbyDat)
  RiesbyDat <- RiesbyDat[!(RiesbyDat$id %in% RiesbyDat[which(is.na(RiesbyDat$depr)), "id"]) , ]
  
  write.table(RiesbyDat, "./data/RiesbyDat.txt")
  read.table("./data/RiesbyDat.txt")
  
# Poling data ####
  # Used by Gelman et al 2014 (p.437, section 16.5)
  polls <- as.data.frame(read_dta(file = "./data/polls.dta"))
  census88 <- as.data.frame(read_dta(file = "./data/census88.dta"))
  
  head(polls)
  head(census88)
  summary(ffc.stata$state)
  
# Thailand education data ####
  # Used by Hox 2010 (p.132, section 6.3)
  # and Raudenbush Bryk 2002 (p.320)
  thaidat <- read.table("./Data/UTHAI1.dat")[,1:5] #long format
  colnames(thaidat) <- c("schid", "male", "pped", "rep1", "msesc") # pped = preschool pupil education)
  summary(thaidat)
  head(thaidat) 
  unique(thaidat$msesc)

  nj <- NULL
  for (i in 1:length(unique(thaidat$schid))) {
    nj[i] <- sum(thaidat$schid == unique(thaidat$schid)[i])
  }
  # Which school size more common
  sort(table(nj),decreasing=TRUE)[1:3]
  J <- 26
  sum(nj == J) # there are 24 schools with 26 students
  school_index <- NULL
  for (i in 1:length(unique(thaidat$schid))) {
    cond <- sum(thaidat$schid == unique(thaidat$schid)[i])
    if(cond == J){
      school_index[i] <- unique(thaidat$schid)[i]
    }
  }
  school_index <- school_index[is.na(school_index) == FALSE]
  thaidat <- thaidat[which(thaidat$schid %in% school_index), ] # select only schools with same nj
  
  # Dichotomize level 2 variable
  thaidat$msesc[which(thaidat$msesc < 0)] <- 0
  thaidat$msesc[which(thaidat$msesc > 0)] <- 1
  
  # Create an "occasion of measurement" vairable
  occ <- rep(c(1:J), nrow(thaidat)/J) - J/2
  occ <- rep(scale(c(1:J)), nrow(thaidat)/J)
  
  yvec <- thaidat$rep1
  Xmat <- cbind(rep(1, nrow(thaidat)), occ, thaidat$msesc, occ*thaidat$msesc)
  dat <- list(yvec, Xmat)

# GPA longitudinal data college students (GPA dichotmoized as below/above average) ####
  # Continuous Outcome Dichotomized
  # Used by Hox 2010 (p. 92, section 5.2)
  gpadata <- as.data.frame(read_sav("./Data/gpa2long.sav"))
  
  # Select Columns
  gpadata <- gpadata[, c(1, 2, 3, 5)]

  # Recode Time
  #schizdata$week <- rep(seq(-3, 3, 1), nrow(schizdata)/7)
  
  #Dichotomize y
  gpadata$gpa[which(gpadata$gpa < mean(gpadata$gpa))] <- 0
  gpadata$gpa[which(gpadata$gpa >= mean(gpadata$gpa))] <- 1
  
  # Occasion centered
  gpadata$occas <- scale(gpadata$occas)
  round(as.vector(unique(scale(gpadata$occas))), 2)
  # Create for functions
  dat <- data.frame(cluster = gpadata$student,
                    yvec    = gpadata$gpa,
                    tvec    = gpadata$occas,gpadata$occas,
                    xvec    = gpadata$sex)
  
  yvec <- gpadata$gpa
  clusters <- gpadata$student
  Xmat <- cbind(rep(1, nrow(gpadata)), gpadata$occas, gpadata$sex, gpadata$occas*gpadata$sex)
  J <- length(unique(gpadata$occas))
  n <- nrow(dat)/J
  dat <- list(yvec, Xmat)
  
# NELS data (students within schools math score) ####
  load("./Data/nelsSES.RData")
  
# Work and Schooling data from Liss panel ####
  #Little within variance in outcome
# Source: https://www.dataarchive.lissdata.nl/study_units/view/17
  WSdat1 <- as.data.frame(read_dta(file = "./Data/cw08a_1.1p_EN.dta"))
  WSdat2 <- as.data.frame(read_dta(file = "./Data/cw09b_EN_3.0p.dta"))
  WSdat3 <- as.data.frame(read_dta(file = "./Data/cw10c_EN_1.0p.dta"))
  WSdat4 <- as.data.frame(read_dta(file = "./Data/cw11d_EN_1.0p.dta"))
  WSdat5 <- as.data.frame(read_dta(file = "./Data/cw12e_EN_1.0p.dta"))
  WSdat6 <- as.data.frame(read_dta(file = "./Data/cw13f_EN_1.0p.dta"))
  
# Select Varibales of interest
  # Codes:
  # - 001 paid work
  # - 002 respondent year birth
  # - 005 highest level educ
  # - 035 followd corse last 12 months
  
  WSdat1$cw08a005[WSdat1$cw08a005 > 26] <- NA
  WSdat2$cw09b005[WSdat2$cw09b005 > 26] <- NA
  WSdat3$cw10c005[WSdat3$cw10c005 > 26] <- NA
  WSdat4$cw11d005[WSdat4$cw11d005 > 26] <- NA
  WSdat5$cw12e005[WSdat5$cw12e005 > 26] <- NA
  WSdat6$cw13f005[WSdat6$cw13f005 > 26] <- NA
  
  WSdat1 <- na.omit(WSdat1[, c(which(grepl("nomem", colnames(WSdat1))),
                       which(grepl("002", colnames(WSdat1))),
                       which(grepl("001", colnames(WSdat1))),
                       which(grepl("005", colnames(WSdat1)))
                       )])
  WSdat2 <- na.omit(WSdat2[, c(which(grepl("nomem", colnames(WSdat2))),
                       which(grepl("002", colnames(WSdat2))),
                       which(grepl("001", colnames(WSdat2))),
                       which(grepl("005", colnames(WSdat2)))
                       )])
  WSdat3 <- na.omit(WSdat3[, c(which(grepl("nomem", colnames(WSdat3))),
                       which(grepl("002", colnames(WSdat3))),
                       which(grepl("001", colnames(WSdat3))),
                       which(grepl("005", colnames(WSdat3)))
                       )])
  WSdat4 <- na.omit(WSdat4[, c(which(grepl("nomem", colnames(WSdat4))),
                       which(grepl("002", colnames(WSdat4))),
                       which(grepl("001", colnames(WSdat4))),
                       which(grepl("005", colnames(WSdat4)))
                       )])
  WSdat5 <- na.omit(WSdat5[, c(which(grepl("nomem", colnames(WSdat5))),
                       which(grepl("002", colnames(WSdat5))),
                       which(grepl("001", colnames(WSdat5))),
                       which(grepl("005", colnames(WSdat5)))
                       )])
  WSdat6 <- na.omit(WSdat6[, c(which(grepl("nomem", colnames(WSdat6))),
                       which(grepl("002", colnames(WSdat6))),
                       which(grepl("001", colnames(WSdat6))),
                       which(grepl("005", colnames(WSdat6)))
                       )])
# Select Recidive Cases
  index <- intersect(
    intersect(
      intersect(
        intersect(
          intersect(WSdat1$nomem_encr, WSdat2$nomem_encr),
                                WSdat3$nomem_encr), 
        WSdat4$nomem_encr),
      WSdat5$nomem_encr),
    WSdat6$nomem_encr)
  WSdat1_1 <- WSdat1[WSdat1$nomem_encr %in% index, ]
  WSdat2_1 <- WSdat2[WSdat2$nomem_encr %in% index, ]
  WSdat3_1 <- WSdat3[WSdat3$nomem_encr %in% index, ]
  WSdat4_1 <- WSdat4[WSdat4$nomem_encr %in% index, ]
  WSdat5_1 <- WSdat5[WSdat5$nomem_encr %in% index, ]
  WSdat6_1 <- WSdat6[WSdat6$nomem_encr %in% index, ]
  
    sum(WSdat1_1$nomem_encr != WSdat6_1$nomem_encr)
  head(WSdat1_1)
  head(WSdat2_1)
  head(WSdat3_1)
  head(WSdat4_1)
  head(WSdat5_1)
  head(WSdat6_1)

  colnames(WSdat1_1) <- colnames(WSdat2_1) <- colnames(WSdat3_1) <- c("id", "yearB", "emp", "edu")
  colnames(WSdat4_1) <- colnames(WSdat5_1) <- colnames(WSdat6_1) <- c("id", "yearB", "emp", "edu")
  dat <- rbind(WSdat1_1, WSdat2_1, WSdat3_1,
               WSdat4_1, WSdat5_1, WSdat6_1)
  dat <- dat[order(dat$id),]
  
  # Check you have same measurment for each id
  nj <- NULL
  for (i in 1:length(unique(dat$id))) {
    nj[i] <- sum(dat$id == unique(dat$id)[i])
  }
  sum(nj != 6)
  # Select only young cases
  dat <- dat[dat$yearB > 1984 & dat$yearB < 1989, ]
  
  edu <- as.numeric(dat$edu)
  
  # edu[edu <= 11] <- 1 # low
  # edu[edu >= 12 & edu <= 21] <- 2 #non-acc
  edu[edu <= 21] <- 1 #non-acc
  edu[edu >= 22 & edu <= 26] <- 2 #acc (including others beacuse they were mistakes)
  
  dat$eduRec <- edu # factor(edu, levels = c(1, 2, 3), labels = c("low", "int", "acc"))
  dat$tvec <- scale(rep(1:6, nrow(dat)/6))
  WSliss <- dat
  
  J <- 6
  n <- nrow(WSliss)/J
  
  ###
  
  # Check you have same measurment for each id
  nj <- NULL
  for (i in 1:length(unique(dat$id))) {
    nj[i] <- sum(dat$id == unique(dat$id)[i])
  }
  sum(nj != 6)
  
  index_educ <- NULL
  ss <- 1
  nrow(dat)/6
  for (i in 1:85) {
    print(ss)
    index_educ[ss:(ss+5)] <- dat$edu[ss] == dat$edu[ss+5]
    print(dat[ss:(ss+5), ])
    print(ss)
    ss <- ss + 6
  }
  dat <- dat[index_educ, ]
  sum(na.omit(index_educ)) # 82 cases in total
  index_educ[which(is.na(index_educ))] <- FALSE
  sum(index_educ)
  
  unique(dat$id)[index_educ]
  
# Simulated data ####
  simdata <- generate_yvec_logreg(n = 30, beta = c(-.625,.25,-.25,.125),
                                  PsiMat = matrix(c(1,-.25,-.25,.5),
                                                  ncol = 2))
  simdata.df <- data.frame(cluster = rep(1:30, each = 7),
                    yvec = simdata[[1]],
                    tvec = simdata[[2]][,2],
                    xvec = simdata[[2]][,3])
  saveRDS(simdata.df, "./Data/simdata_n30_j7_Psi1-neg0_25-0_5.rds")
  simdata.df <- readRDS("./Data/simdata_n30_j7_Psi1-neg0_25-0_5.rds")
  