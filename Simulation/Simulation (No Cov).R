### Load the UK Data----
setwd("~/desktop/Dissertation/Data")
UK6 <- read.csv(file = "UK6.csv") # With no NA cell
UK5 <- read.csv(file = "UK5.csv") # With the NA cell
UK5_complete <- read.csv(file = "UK5_complete.csv") # Used for dga package, replaced the NA count to 0
# Add the NA column to data
UK5_withNA <- rbind(c(0, 0, 0, 0, 0, NA), UK5)
UK5_withNA


library(mefa)
## Load the data
setwd("~/desktop/Dissertation/Data")
RM <- read.csv(file = "Romania 2015-2016.csv")[, -1]

## Data Handling for Romania data -----
# Convert all categorical variables to factors to fit Poisson/NB models
RM_3fac <- RM
RM_3fac[, 6:10] <- data.frame(lapply(RM_3fac[, 6:10], as.factor))
str(RM_3fac)

# Separate data by each covariate
# rm(sum)
aggregate(Freq ~ S, RM_3fac, FUN = sum) 
aggregate(Freq ~ A, RM_3fac, FUN = sum) 
aggregate(Freq ~ A + Y, RM_3fac, FUN = sum) 
aggregate(Freq ~ E, RM_3fac, FUN = sum) 
aggregate(Freq ~ D, RM_3fac, FUN = sum) 
aggregate(Freq ~ D + Y, RM_3fac, FUN = sum) 

# Separate data by year 
aggregate(Freq ~ Y, RM_3fac, FUN = sum) 

RM_3fac_15 <- RM_3fac[RM_3fac$Y == "2015", ]
length(which(RM_3fac_15$Freq == 0)) 
length(which(RM_3fac_15$Freq <= 10 & RM_3fac_15$Freq > 0)) 

RM_3fac_16 <- RM_3fac[RM_3fac$Y == "2016", ]
length(which(RM_3fac_16$Freq == 0)) 
length(which(RM_3fac_16$Freq <= 10 & RM_3fac_16$Freq > 0)) 


# For two-factor treatment, replace all table entries to 0 and 1 
for (i in c("male", "adult", "other", "internal")){
  RM[RM == i] <- 0
}
for (i in c("beggary")){
  RM[RM == i] <- 0 
  # Treat 'beggary' & 'other' together as 'non-sexual' exploitation
}
for (i in c("female", "minor", "sexual", "external")){
  RM[RM == i] <- 1
}
str(RM)
# Convert all factors to numeric to fit other models
RM <- data.frame(lapply(RM, as.numeric))


# For comparison across models, only allow 0 and 1 for exploitation type
# As packages for Bayes threshold approach/Bayes log-linear model only work with tables of 0 and 1 
RM_2fac <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D + Y, RM, FUN = sum)
# Replace all unobserved cells with NA 
NArows <- which(RM_2fac$R1==0 & RM_2fac$R2==0 & RM_2fac$R3==0 & RM_2fac$R4==0 & RM_2fac$R5==0)
RM_withNA <- RM_2fac
RM_withNA$Freq[NArows] <- NA
# Remove all NA cells
RM_2fac <- RM_withNA[-which(is.na(RM_withNA$Freq)), ]
length(which(RM_2fac$Freq == 0)) 
length(which(RM_2fac$Freq <= 5 & RM_2fac$Freq > 0)) 
# Separate data with NA by year
RM_withNA_15 <- RM_withNA[RM_withNA$Y == "2015", -10]
RM_withNA_16 <- RM_withNA[RM_withNA$Y == "2016", -10]
RM_withNA_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D, 
                              RM_withNA, FUN = sum, na.action = na.pass)

# Separate data by year
RM_2fac_15 <- RM_2fac[RM_2fac$Y == "2015", -10]
length(which(RM_2fac_15$Freq == 0)) 
length(which(RM_2fac_15$Freq <= 10 & RM_2fac_15$Freq > 0)) 

RM_2fac_16 <- RM_2fac[RM_2fac$Y == "2016", -10]
length(which(RM_2fac_16$Freq == 0)) 
length(which(RM_2fac_16$Freq <= 10 & RM_2fac_16$Freq > 0)) 

# Sum over years - still many zero counts
RM_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D, RM_2fac, FUN = sum)
length(which(RM_noyear$Freq == 0)) 
length(which(RM_noyear$Freq <= 5 & RM_noyear$Freq > 0)) 

# Sum over all covariates - less zero counts but still quite sparse
RM_nocov <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac, FUN = sum)
length(which(RM_nocov$Freq == 0)) 
length(which(RM_nocov$Freq <= 5 & RM_nocov$Freq > 0)) 

# Prepare data with no covariates for each year
RM_nocov_15 <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac_15, FUN = sum)
RM_nocov_16 <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac_16, FUN = sum)




## Data Handling for Serbia data -----
# Replace all unobserved cells with NA 
# Convert all categorical variables to factors to fit Poisson/NB models
SB_3fac <- SB
SB_3fac[, 5:8] <- data.frame(lapply(SB_3fac[, 5:8], as.factor))
# Change covariate names 
levels(SB_3fac$S) <- c("female", "male")
levels(SB_3fac$A) <- c("adult", "minor")
levels(SB_3fac$N) <- c("non-Serbian", "Serbian")
levels(SB_3fac$E) <- c("beggary", "other", "sexual")
str(SB_3fac)

# Separate data by each covariate
aggregate(Freq ~ S, SB_3fac, FUN = sum) 
# 365 female, 509 male
aggregate(Freq ~ A, SB_3fac, FUN = sum) 
# 626 adult, 248 minor
aggregate(Freq ~ A + Y, SB_3fac, FUN = sum) 
aggregate(Freq ~ E, SB_3fac, FUN = sum) 
# 230 sexual, 589 other, 55 beggary
aggregate(Freq ~ N, SB_3fac, FUN = sum) 
# 776 Serbian, 98 non-Serbian
aggregate(Freq ~ N + Y, SB_3fac, FUN = sum) 
# Separate data by year 
aggregate(Freq ~ Y, SB_3fac, FUN = sum) 
# 156 obs for 2013, 480 obs for 2014
# 96 obs for 2015, 142 obs for 2016

SB_3fac_13 <- SB_3fac[SB_3fac$Y == "2013", ]
length(which(SB_3fac_13$Freq == 0)) 
# 346 out of 384 cells with zero counts
length(which(SB_3fac_13$Freq <= 10 & SB_3fac_13$Freq > 0)) 
# 34 out of 384 cells have counts less than 10

SB_3fac_14 <- SB_3fac[SB_3fac$Y == "2014", ]
length(which(SB_3fac_14$Freq == 0)) 
# 351 out of 384 cells with zero counts
length(which(SB_3fac_14$Freq <= 10 & SB_3fac_14$Freq > 0)) 
# 28 out of 384 cells have counts less than 10

SB_3fac_15 <- SB_3fac[SB_3fac$Y == "2015", ]
length(which(SB_3fac_15$Freq == 0)) 
# 352 out of 384 cells with zero counts
length(which(SB_3fac_15$Freq <= 10 & SB_3fac_15$Freq > 0)) 
# 31 out of 384 cells have counts less than 10

SB_3fac_16 <- SB_3fac[SB_3fac$Y == "2016", ]
length(which(SB_3fac_16$Freq == 0)) 
# 346 out of 384 cells with zero counts
length(which(SB_3fac_16$Freq <= 10 & SB_3fac_16$Freq > 0)) 
# 36 out of 384 cells have counts less than 10

# For two-factor treatment, replace all table entries to 0 and 1 
SB <- SB_3fac
SB[, 5:8] <- data.frame(lapply(SB[, 5:8], as.character))
for (i in c("male", "adult", "other", "Serbian")){
  SB[SB == i] <- 0
}
for (i in c("beggary")){
  SB[SB == i] <- 0 
  # Treat 'beggary' & 'other' together as 'non-sexual' exploitation
}
for (i in c("female", "minor", "sexual", "non-Serbian")){
  SB[SB == i] <- 1
}
str(SB)
# Convert all factors to numeric to fit other models
SB <- data.frame(lapply(SB, as.numeric))


# For comparison across models, only allow 0 and 1 for exploitation type
# As packages for Bayes threshold approach/Bayes log-linear model only work with tables of 0 and 1 
SB_2fac <- aggregate(Freq ~ R1 + R2 + R3 + R4 + S + A + N + E + Y, SB, FUN = sum)
# Remove year 2014 data
SB_2fac <- SB_2fac[-which(SB_2fac$Y == '2014'), ]
# Replace all unobserved cells with NA 
NArows <- which(SB_2fac$R1==0 & SB_2fac$R2==0 & SB_2fac$R3==0 & SB_2fac$R4==0)
SB_withNA <- SB_2fac
SB_withNA$Freq[NArows] <- NA
# Remove all NA cells
SB_2fac <- SB_withNA[-which(is.na(SB_withNA$Freq)), ]
length(which(SB_2fac$Freq == 0)) 
# 625 out of 720 cells have zero counts
length(which(SB_2fac$Freq <= 10 & SB_2fac$Freq > 0)) 
# 86 out of 720 cells have sparse counts
# Separate data with NA by year
SB_withNA_13 <- SB_withNA[SB_withNA$Y == "2013", -9]
SB_withNA_15 <- SB_withNA[SB_withNA$Y == "2015", -9]
SB_withNA_16 <- SB_withNA[SB_withNA$Y == "2016", -9]
SB_withNA_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + S + A + E + N, 
                              SB_withNA, FUN = sum, na.action = na.pass)


# Separate data by year
SB_2fac_13 <- SB_2fac[SB_2fac$Y == "2013", -9]
length(which(SB_2fac_13$Freq == 0)) 
# 206 out of 236 cells with zero counts
length(which(SB_2fac_13$Freq <= 10 & SB_2fac_13$Freq > 0)) 
# 26 out of 236 cells have counts less than 10

SB_2fac_15 <- SB_2fac[SB_2fac$Y == "2015", -9]
length(which(SB_2fac_15$Freq == 0)) 
# 206 out of 236 cells with zero counts
length(which(SB_2fac_15$Freq <= 10 & SB_2fac_15$Freq > 0)) 
# 28 out of 236 cells have counts less than 10

SB_2fac_16 <- SB_2fac[SB_2fac$Y == "2016", -9]
length(which(SB_2fac_16$Freq == 0)) 
# 206 out of 238 cells with zero counts
length(which(SB_2fac_16$Freq <= 10 & SB_2fac_16$Freq > 0)) 
# 29 out of 238 cells have counts less than 10

# Sum over years - still many zero counts
SB_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + S + A + N + E, SB_2fac, FUN = sum)
length(which(SB_noyear$Freq == 0)) 
# 188 out of 240 cells with zero counts
length(which(SB_noyear$Freq <= 10 & SB_noyear$Freq > 0)) 
# 42 out of 240 cells have counts less than 10

# Sum over all covariates - less zero counts but still quite sparse
SB_nocov <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac, FUN = sum)
length(which(SB_nocov$Freq == 0)) 
# 5 out of 15 cells with zero counts
length(which(SB_nocov$Freq <= 10 & SB_nocov$Freq > 0)) 
# 6 out of 15 cells have counts less than 10

# Prepare data with no covariates for each year
SB_nocov_13 <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac_13, FUN = sum)
SB_nocov_15 <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac_15, FUN = sum)
SB_nocov_16 <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac_16, FUN = sum)




RM_nocov_NA <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_withNA, FUN = sum, na.action = na.pass)
SB_nocov_NA <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_withNA, FUN = sum, na.action = na.pass)
str(RM_nocov_NA)
str(RM_nocov)
mean(RM_nocov$Freq) # 52.74
sd(RM_nocov$Freq) #255.69

### Parametric bootstrap: simulate datasets similar to the base data------
# Use the model with two interactions selected in AIC/BIC as the generating model (up to all two-way interactions) 
# RM.poi.full <- glm(Freq ~ (R1 + R2 + R3 + R4 + R5)^2, data = RM_nocov, family = poisson)
RM.poi.start <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + R2:R4, data = RM_nocov, family = poisson) 
# RM.poi.BIC <- step(RM.poi.start, scope = formula(RM.poi.full), direction = "forward", k = log(nrow(RM_nocov))) 
obs <- sum(RM_nocov$Freq)

## Simulate based on stepwiseBIC model
nsims <- 100
unobs.pred <- predict.glm(RM.poi.start, RM_nocov_NA[is.na(RM_nocov_NA$Freq), ], type = "response") 
npop <- obs + sum(round(unobs.pred))
RM_nocov_fill <- RM_nocov_NA
RM_nocov_fill$Freq[is.na(RM_nocov_fill$Freq)] <- round(unobs.pred) 
RM.poi.fill <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + R2:R4, data = RM_nocov_fill, family = poisson) 
cell.counts <- predict.glm(RM.poi.fill, RM_nocov_fill, type = "response")
# Calculate the probabilities of capture history
cell.prob <- cell.counts / sum(cell.counts)
# Simulation step
set.seed(124)
# Generate bootstrap samples from multinomial distribution with observed cell probabilities
realisations <- as.matrix(rmultinom(nsims, npop, cell.prob), ncol = nsims)
simdata <- realisations[-which(is.na(RM_nocov_NA$Freq)), ] 
# Store simulated datasets (without the NA cells) in a list
lst <- replicate(nsims, RM_nocov, simplify = FALSE) 
for (i in 1:nsims){
  lst[[i]]$Freq <- simdata[, i]
}
# Check the sparsity for simulated datasets 
zero.cell <- sapply(lst, function(x){length(which(x$Freq == 0))}) 
small.cell <- sapply(lst, function(x){length(which(x$Freq <= 5 & x$Freq > 0))}) 
hist(zero.cell/nrow(RM_nocov))
hist(small.cell/nrow(RM_nocov))
mean(zero.cell/nrow(RM_nocov))
mean(small.cell/nrow(RM_nocov))


## Simulate less sparse datasets based on data simulated by the stepwise BIC model----
set.seed(123)
add_random_number <- function(ds, p){ 
  n <- length(which(ds$Freq <= 5))
  select <- sample(which(ds$Freq <= 5), round(n/p), replace = F)  
  ds$Freq[select] <- ds$Freq[select] + sample(2:8, round(n/p), replace = TRUE)
  return(ds)
}
RM_nocov_NA_S <- add_random_number(RM_nocov_NA, 4) # Change the value of p
RM_nocov_S <- RM_nocov_NA_S[-which(is.na(RM_nocov_NA_S$Freq)), ]

# Use the stepwise BIC Poisson model as the generating model (up to all two-way interactions) 
# RM.poi.full_S <- glm(Freq ~ (R1 + R2 + R3 + R4 + R5 + S + A + E + D)^2, data = RM_noyear_S, family = poisson)
RM.poi.start_S <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + R2:R4, data = RM_nocov_S, family = poisson) 
# RM.poi.BIC_S <- step(RM.poi.start_S, scope = formula(RM.poi.full_S), direction = "forward", k = log(nrow(RM_noyear_S))) 
obs_S <- sum(RM_nocov_S$Freq)

## Simulate based on stepwiseBIC model
# nsims <- 1000
nsims <- 100
unobs.pred_S <- predict.glm(RM.poi.start_S, RM_nocov_NA_S[is.na(RM_nocov_NA_S$Freq), ], type = "response") 
npop_S <- obs_S + sum(round(unobs.pred_S))
RM_nocov_fill_S <- RM_nocov_NA_S
RM_nocov_fill_S$Freq[is.na(RM_nocov_fill_S$Freq)] <- round(unobs.pred) 
RM.poi.fill_S <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + R2:R4, data = RM_nocov_fill_S, family = poisson) 
cell.counts <- predict.glm(RM.poi.fill_S, RM_nocov_fill_S, type = "response")
# Calculate the probabilities of capture history
cell.prob <- cell.counts / sum(cell.counts)
# Simulation step
set.seed(124)
# Generate bootstrap samples from multinomial distribution with observed cell probabilities
realisations <- as.matrix(rmultinom(nsims, npop_S, cell.prob), ncol = nsims)
simdata <- realisations[-which(is.na(RM_nocov_NA_S$Freq)), ] 
# Store simulated datasets (without the NA cells) in a list
lst_S <- replicate(nsims, RM_nocov_S, simplify = FALSE) 
for (i in 1:nsims){
  lst_S[[i]]$Freq <- simdata[, i]
}
# Check the sparsity for simulated datasets 
# NA.cell <- sapply(lst, function(x){length(which(is.na(x$Freq)))}) # Always 1 for all simulated datasets
zero.cell_S <- sapply(lst_S, function(x){length(which(x["Freq"] == 0))}) 
small.cell_S <- sapply(lst_S, function(x){length(which(x["Freq"] <= 5 & x["Freq"] > 0))}) 
hist(zero.cell_S/nrow(RM_noyear_S))
hist(small.cell_S/nrow(RM_noyear_S))
mean(zero.cell_S/nrow(RM_noyear_S))
mean(small.cell_S/nrow(RM_noyear_S))





#### Simulation: compare the performance of six methods (no NB if tested with no over-dispersion)-----
### 1. Freq Poisson: stepwise forward selection with BIC------
# (1 hour to run on 1k simulated datasets)
## Define functions to estimate N from a Poisson model
estimate_N <- function(o, model){
  coef <- coef(model)
  unobs <- predict.glm(model, RM_nocov_NA[is.na(RM_nocov_NA$Freq), ], type = "response")
  Nhat <- sum(unobs, o)
  return(Nhat)
}

estimate_m1 <- function(ds){
  ds.poi.start <- glm(Freq ~ R1 + R2 + R3 + R4 + R5, data = ds, family = poisson) 
  ds.poi.full <- glm(Freq ~ (R1 + R2 + R3 + R4 + R5)^2, data = ds, family = poisson)
  NA.len <- length(which(is.na(ds$Freq))) 
  ds.poi.BIC <- step(ds.poi.start, scope = formula(ds.poi.full), direction = "forward") 
  obs <- sum(ds$Freq, na.rm = TRUE) 
  Nhat.m1 <- estimate_N(obs, ds.poi.BIC)
  
  
  # Estimate the BCa CI for each simulated dataset
  # Generate bootstrap samples from multinomial distribution with observed cell probabilities
  nsims.b <- 40
  cell.prob.b <- ds$Freq[!is.na(ds$Freq)]/obs
  # Simulation step
  set.seed(101)
  simdata.b <- as.matrix(rmultinom(nsims.b, obs, cell.prob.b), ncol = nsims.b)
  lst.b <- replicate(nsims.b, ds, simplify = FALSE)
  for (i in 1:nsims.b){
    lst.b[[i]]$Freq[!is.na(lst.b[[i]]$Freq)] <- simdata.b[, i]
  }
  fit.models <- lapply(lst.b, function(x){glm(ds.poi.BIC$formula, family = poisson, data = x)})
  Nhat.b <- sapply(fit.models, function(x) estimate_N(obs, x))
  # Compute the BCa bootstrap CI for each simulated dataset
  CI.m1 <- quantile(Nhat.b, c(0.025, 0.975))
  # library(coxed)
  # CI.m1 <- as.numeric(bca(Nhat.b, conf.level = 0.95))
  
  # Test for Over-dispersion by Cameron & Trivedi (1990), this can decide to use Negtive Binomial model or not
  library(AER)
  dispersion_m1 <- dispersiontest(ds.poi.BIC, trafo = 2) 
  # Return the point estimate, CI & p-value for dispersion test
  results <- c(Nhat.m1, CI.m1, dispersion_m1$p.value)
  return(results)
}

## Assess the performance of this model on simulated datasets 
# For very sparse datasets
re.m1 <- sapply(lst, estimate_m1)
mean(re.m1[1, ])
rowMeans(re.m1[c(2,3), ])
# hist(re.m1[1, ])
# abline(v = npop, col = "red")
# hist(re.m1[4, ])
(sum(which(re.m1[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of mean square error for log of Nhat
log(mean((log(re.m1[1, ]) - log(npop))^2)) 
# Coverage rate
(length(which(re.m1[2, ] <= npop & re.m1[3, ] >= npop)))/nsims


# For less sparse datasets
re.m1 <- sapply(lst_S, estimate_m1)
mean(re.m1[1, ])
rowMeans(re.m1[c(2,3), ])
# hist(re.m1[1, ])
# abline(v = npop, col = "red")
# hist(re.m1[4, ])
(sum(which(re.m1[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of mean square rrror for log of Nhat
log(mean((log(re.m1[1, ]) - log(npop_S))^2)) 
# Coverage rate
(length(which(re.m1[2, ] <= npop_S & re.m1[3, ] >= npop_S)))/nsims



## 2. Freq Poisson: stepwise selection algorithm for sparse data ------
# (50h to run on 1k simulated datasets)
library(SparseMSE)
estimate_m2 <- function(ds){
  estimates <- estimatepopulation(ds, nboot = 20, method = "stepwise", mX = NULL, pthresh = 0.02, iseed = 1234,
                     alpha = c(0.025, 0.975)) 
  # Use the default threshold p-value 0.02
  # nboot defines the number of boostraps for BCabootstrap CI
  obs <- sum(ds$Freq, na.rm = TRUE) 
  Nhat.m2 <- sum(obs, estimates$popest)
  CI.m2 <- estimates$BCaquantiles
  
  # Test for Over-dispersion by Cameron & Trivedi (1990)
  # library(AER)
  # dispersion_m2 <- dispersiontest(estimates$MSEfit$fit, trafo = 1) 
  # results <- c(Nhat_m2, CI_m2, dispersion_m2$p.value)
  # Return the point estimate & CI 
  results <- c(Nhat.m2, CI.m2)
  return(results)
}



re.m2 <- sapply(lst, estimate_m2)
mean(re.m2[1, ])
rowMeans(re.m2[c(2,3), ])
# hist(re.BIC.m2[1, ])
# abline(v = npop.BIC, col = "red")
# hist(re.BIC.m2[4, ])
# (sum(which(re.BIC.m2[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
log(mean((log(re.m2[1, ]) - log(npop))^2)) 
# Coverage rate
(length(which(re.m2[2, ] <= npop & re.m2[3, ] >= npop)))/nsims



re.m2 <- sapply(lst_S, estimate_m2)
mean(re.m2[1, ])
rowMeans(re.m2[c(2,3), ])
# hist(re.BIC.m2[1, ])
# abline(v = npop.BIC, col = "red")
# hist(re.BIC.m2[4, ])
# (sum(which(re.BIC.m2[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
log(mean((log(re.m2[1, ]) - log(npop_S))^2)) 
# Coverage rate
(length(which(re.m2[2, ] <= npop_S & re.m2[3, ] >= npop_S)))/nsims



## 3. Freq Negative Binomial: stepwise forward selection with BIC-----
# (? to run on 1k simulated datasets)
# NB models are used for over-dispersed count data, when the conditional variance exceeds the conditional mean
# It is a generalization of Poisson regression with an extra dispersion parameter to model the over-dispersion
library(MASS)
estimate_m4 <- function(ds){
  ds.nb.start <- glm.nb(Freq ~ LA + NG + PF + GO + GP, data = ds) 
  ds.nb.full <- glm.nb(Freq ~ (LA + NG + PF + GO + GP)^2, data = ds)
  NA.len <- length(which(is.na(ds$Freq))) 
  ds.nb.BIC <- step(ds.nb.start, scope = formula(ds.nb.full), direction = "forward", k = log(nrow(ds) - NA.len)) # Specify k = log(n) to use BIC (deleted 1 missing obs)
  obs <- sum(ds$Freq, na.rm = TRUE) 
  Nhat.m4 <- estimate_N(ds.nb.BIC)
  
  
  # Estimate the BCa CI for each simulated dataset
  # Generate 1000 bootstrap samples from multinomial distribution with observed cell probabilities
  # nsims.b <- 1000
  nsims.b <- 100
  cell.prob.b <- ds$Freq[!is.na(ds$Freq)]/obs
  # Simulation step
  set.seed(102)
  simdata.b <- as.matrix(rmultinom(nsims.b, obs, cell.prob.b), ncol = nsims.b)
  lst.b <- replicate(nsims.b, ds, simplify = FALSE)
  for (i in 1:nsims.b){
    lst.b[[i]]$Freq[!is.na(lst.b[[i]]$Freq)] <- simdata.b[, i]
  }
  fit.models <- lapply(lst.b, function(x){glm.nb(formula(ds.nb.BIC), data = x)})
  Nhat.b <- sapply(fit.models, estimate_N)
  # Compute the BCa bootstrap CI for each simulated dataset
  library(coxed)
  CI.m4 <- as.numeric(bca(Nhat.b, conf.level = 0.95))
  # Return the point estimate & CI 
  results <- c(Nhat.m4, CI.m4)
  return(results)
}

# Assess the performance of this model on simulated datasets
# For datasets based on main effect model
re.main.m4 <- sapply(lst.main, estimate_m4)
hist(re.main.m4[1, ])
abline(v = npop.main, col = "red")
hist(re.main.m4[4, ])
(sum(which(re.main.m4[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
logMSE.main.m4 <- log(mean((log(re.main.m4[1, ]) - log(npop.main))^2)) #
# Coverage rate
covrate.main.m4 <- (length(which(re.main.m4[2, ] <= npop.main & re.main.m4[3, ] >= npop.main)))/nsims

# For datasets based on model with three interactions
re.md.m4 <- sapply(lst.md, estimate_m4)
hist(re.md.m4[1, ])
abline(v = npop.md, col = "red")
hist(re.md.m4[4, ])
(sum(which(re.md.m4[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
logMSE.md.m4 <- log(mean((log(re.md.m4[1, ]) - log(npop.md))^2)) 
# Coverage rate
covrate.md.m4 <- (length(which(re.md.m4[2, ] <= npop.md & re.md.m4[3, ] >= npop.md)))/nsims

# For datasets based on stepwise BIC
re.BIC.m4 <- sapply(lst.BIC, estimate_m4)
hist(re.BIC.m4[1, ])
abline(v = npop.BIC, col = "red")
hist(re.BIC.m4[4, ])
(sum(which(re.BIC.m4[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
logMSE.BIC.m4 <- log(mean((log(re.BIC.m4[1, ]) - log(npop.BIC))^2)) 
# Coverage rate
covrate.BIC.m4 <- (length(which(re.BIC.m4[2, ] <= npop.BIC & re.BIC.m4[3, ] >= npop.BIC)))/nsims


## 4. Bayes GM model with no strata ------
# (70min to run on 1k simulated datasets)
library(dga)
estimate_m3 <- function(ds){
  # The dataframe must be sorted by lexicographical order of cells 
  ds <- ds[order(R1, R2, R3, R4, R5),]
  # Replace the NA counts to 0
  ds <- rbind(c(0, 0, 0, 0, 0, 0), ds) 
  Y <- array(ds$Freq, dim=c(2,2,2,2,2))
  # Nmissing is up to ten times the number of observed values
  Nmissing <- 1:(10 * sum(Y)) 
  num.lists <- 5                 
  delta <- 1/(2 ^ num.lists) # Default value of delta
  data(graphs5)
  weights <- bma.cr(Y, Nmissing, delta, graphs5)
  # Draw a sample of 1000 values from the estimated (averaged) posterior of N
  post_probs <- apply(weights, 2, sum) 
  m0s_sample_from_post <- sample(Nmissing, size = 1000, replace = TRUE, prob = post_probs)
  nhat_sample_from_post <- m0s_sample_from_post + sum(Y) 
  # hist(nhat_sample_from_post, breaks = 50, main = "Posterior Distribution of the Total Population N") # Histogram of posterior samples of N
  Nhat.m3 <- quantile(nhat_sample_from_post, 0.5) # Posterior median of N
  CI.m3 <- quantile(nhat_sample_from_post, c(0.025, 0.975))
  # Return the point estimate & CI 
  results <- c(Nhat.m3, CI.m3)
  return(results)
}

# Assess the performance of this model on simulated datasets
# For datasets based on main effect model
re.m3 <- sapply(lst, estimate_m3)
mean(re.m3[1, ])
rowMeans(re.m3[c(2,3), ])
# hist(re.main.m5[1, ])
# abline(v = npop.main, col = "red")
# hist(re.main.m5[4, ])
# (sum(which(re.main.m5[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
log(mean((log(re.m3[1, ]) - log(npop))^2)) #
# Coverage rate
(length(which(re.m3[2, ] <= npop & re.m3[3, ] >= npop)))/nsims



re.m3 <- sapply(lst_S, estimate_m3)
mean(re.m3[1, ])
rowMeans(re.m3[c(2,3), ])
# hist(re.main.m5[1, ])
# abline(v = npop.main, col = "red")
# hist(re.main.m5[4, ])
# (sum(which(re.main.m5[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
log(mean((log(re.m3[1, ]) - log(npop_S))^2)) #
# Coverage rate
(length(which(re.m3[2, ] <= npop_S & re.m3[3, ] >= npop_S)))/nsims


## 5. Bayes log-linear ------
# (25h to run on 1k simulated datasets)
library(conting)
set.seed(111)
estimate_m6 <- function(ds){
  # Set the maximal model in bict() to include all two-way interactions 
  # Use the default Sabanes Bove and Held prior distribution
  # Request 200k MCMC iterations and save the MCMC output to external ﬁles every 1000 iterations
  ex <- bict(formula = Freq ~ (LA + NG + PF + GO + GP)^2, data = ds, 
                  n.sample = 200000, save = 1000, null.move.prob = 0.5) # Use default values a = b = 10^(-3)
  
  # Set a burn-in phase of 20000 iterations (i.e., discard the ﬁrst 10% of the iterations) 
  # Set thinning at every 5th iteration due to the size of MCMC output 
  # Use the Freeman-Tukey discrepancy statistic to assess model adequacy
  ex_bp <- bayespval(ex, n.burnin = 20000, thin = 5, statistic = "FreemanTukey") # Slow to run
  # Derive an MCMC sample from the posterior distribution of the total population size (taking account of different models) and find its posterior mean and 95% HPDI 
  ex_tot <- total_pop(ex, n.burnin = 20000, thin = 5)
  Nhat.m6 <- median(scot_ex_tot$TOT)
  CI.m6 <- quantile(scot_ex_tot$TOT, c(0.025, 0.975))
  # Return the point estimate & CI 
  results <- c(Nhat.m6, CI.m6)
  return(results)
}

# Assess the performance of this model on simulated datasets
# For datasets based on main effect model
re.main.m6 <- sapply(lst.main, estimate_m6)
hist(re.main.m6[1, ])
abline(v = npop.main, col = "red")
hist(re.main.m6[4, ])
(sum(which(re.main.m6[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
logMSE.main.m6 <- log(mean((log(re.main.m6[1, ]) - log(npop.main))^2)) #
# Coverage rate
covrate.main.m6 <- (length(which(re.main.m6[2, ] <= npop.main & re.main.m6[3, ] >= npop.main)))/nsims

# For datasets based on model with three interactions
re.md.m6 <- sapply(lst.md, estimate_m6)
hist(re.md.m6[1, ])
abline(v = npop.md, col = "red")
hist(re.md.m6[4, ])
(sum(which(re.md.m6[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
logMSE.md.m6 <- log(mean((log(re.md.m6[1, ]) - log(npop.md))^2)) 
# Coverage rate
covrate.md.m6 <- (length(which(re.md.m6[2, ] <= npop.md & re.md.m6[3, ] >= npop.md)))/nsims

# For datasets based on stepwise BIC
re.BIC.m6 <- sapply(lst.BIC, estimate_m6)
hist(re.BIC.m6[1, ])
abline(v = npop.BIC, col = "red")
hist(re.BIC.m6[4, ])
(sum(which(re.BIC.m6[4, ] < 0.05)))/nsims # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of Mean Square Error for log of Nhat
logMSE.BIC.m6 <- log(mean((log(re.BIC.m6[1, ]) - log(npop.BIC))^2)) 
# Coverage rate
covrate.BIC.m6 <- (length(which(re.BIC.m6[2, ] <= npop.BIC & re.BIC.m6[3, ] >= npop.BIC)))/nsims



## 6. Bayes non-parametric latent class model ------
# (2h to run on 1k simulated datasets)
library(LCMCR)
set.seed(111)
estimate_m5 <- function(ds){
  # LCMCR requires that the list-membership columns be factors
  for (col in c('R1', 'R2', 'R3', 'R4', 'R5')) {
    ds[[col]] <- as.factor(ds[[col]])
  }
  sampler <- lcmCR(captures = ds, tabular = TRUE, K = 32, a_alpha = 0.25, b_alpha = 0.25,
                   seed = "auto", buffer_size = 10000, thinning = 100) # Default values of a and b
  # Start with K* = 32, equal to the number of unique capture histories
  # Thinning is set for the tracing buffer
  N <- lcmCR_PostSampl(sampler, burnin = 10000, samples = 5000, thinning = 100)
  Nhat.m5 <- quantile(N, 0.5)
  CI.m5 <- quantile(N, c(0.025, 0.975))
  # Return the point estimate & CI 
  results <- c(Nhat.m5, CI.m5)
  return(results)
}

# Assess the performance of this model on simulated datasets
# For datasets based on main effect model
re.m5 <- sapply(lst, estimate_m5)
mean(re.m5[1,])
rowMeans(re.m5[c(2,3),])
# Compute the log of Mean Square Error for log of Nhat
log(mean((log(re.m5[1, ]) - log(npop))^2)) #
# Coverage rate
(length(which(re.m5[2, ] <= npop & re.m5[3, ] >= npop)))/nsims

# Assess the performance of this model on simulated datasets
# For datasets based on main effect model
re.m5 <- sapply(lst_S, estimate_m5)
# Compute the log of Mean Square Error for log of Nhat
log(mean((log(re.m5[1, ]) - log(npop_S))^2)) #
# Coverage rate
(length(which(re.m5[2, ] <= npop_S & re.m5[3, ] >= npop_S)))/nsims