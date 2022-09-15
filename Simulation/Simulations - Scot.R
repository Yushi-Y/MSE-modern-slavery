## Data setup
library(conting)
library(mefa)
library(data.table)
data("ScotPWID")
scot <- ScotPWID # This is our base data

## Data Handling ------ 
# Replace the table entries (characters) to 0 and 1 
scot <- data.frame(lapply(scot, as.character), stringsAsFactors = FALSE)
for (i in c("un", "GGC", "Male", "Young")){
  scot[scot == i] <- 0
}
for (i in c("obs", "Rest", "Female", "Old")){
  scot[scot == i] <- 1
}
str(scot)
Freq <- as.numeric(scot$y)
scot <- data.frame(lapply(scot[, -1], as.factor), Freq)

# Marginalise the table over the covariates
Scot <- data.frame(lapply(scot, as.numeric))
Scot[, -8] <- Scot[, -8] - 1
# Remove individuals not recorded on any list
Scot <- Scot[apply(Scot[,1:4], 1, sum) > 0, ]
# Repeat each capture history by the time of its frequencies (to make the table the right form)
Scot <- Scot[rep(seq_len(nrow(Scot)), times = Scot$Freq), -8]
overlaps <- Scot[, 1:4]
# Count the number of same rows in the dataframe
scot.list <- as.data.frame(setDT(overlaps)[, list(Freq = .N), names(overlaps)])
scot.list


## Simulate datasets similar to the base data ------
# Fit Poisson Log-linear Models on the base data (up to all two-way interactions) 
# Use stepwise forward selection with highest BIC as the generating model
scot.poi.start <- glm(Freq ~ S1 + S2 + S3 + S4, data = scot, family = poisson) # Start with the independence model (with no interactions)
scot.poi <- glm(Freq ~ (S1 + S2 + S3 + S4 + Region + Gender + Age)^2, data = scot, family = poisson)
BICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi), direction = "forward", k = log(nrow(scot) - 8)) # Specify k = log(n) to use BIC (deleted 8 missing obs)
coef <- coef(BICmodel.poi)

# Parametric bootstrap: simulate new data sets similar to the structure of the base data (here is scot)
# In practice, generate bootstrap samples directly from multinomial distribution with observed cell probabilities
obs <- sum(scot$Freq, na.rm = TRUE)
nsims <- 1000
# scot[which(is.na(scot$Freq)), ] # Check all the NA cells
unobs.pred <- c(exp(coef[1]), exp(sum(coef[c(1, 8)])), exp(sum(coef[c(1, 7)])), exp(sum(coef[c(1, 7, 8, 13)])),
                exp(sum(coef[c(1, 6)])), exp(sum(coef[c(1, 6, 8, 12)])),
                exp(sum(coef[c(1, 6, 7)])), exp(sum(coef[c(1, 6:8, 12:13)]))) # Match the order with the order of NA cells
cell.counts <- numeric(2^(ncol(scot) - 1))
cell.counts[which(is.na(scot$Freq))] <- unobs.pred
cell.counts[which(!is.na(scot$Freq))] <- BICmodel.poi$fitted.values
npop <- obs + round(sum(unobs.pred))
# Calculate the probs of capture history
cell.prob <- cell.counts / sum(cell.counts)
# Simulation step
set.seed(001)
realisations <- as.matrix(rmultinom(nsims, npop, cell.prob), ncol = nsims)
simdata <- realisations[which(!is.na(scot$Freq)), ] # Drop the rows that are supposed to have NA cells
lst <- replicate(nsims, scot, simplify = FALSE)
for (i in 1:nsims){
  lst[[i]]$Freq[!is.na(lst[[i]]$Freq)] <- simdata[, i]
}
# Check the sparsity for simulated datasets 
# NA.cell <- sapply(lst, function(x){length(which(is.na(x$Freq)))}) # Always 8 for all simulated datasets
zero.cell <- sapply(lst, function(x){length(which(x$Freq == 0))}) 
small.cell <- sapply(lst, function(x){length(which(x$Freq <= 10 & x$Freq > 0))}) 
hist(zero.cell)
hist(small.cell) # Both similar to the scot dataset


### Simulation 1: compare the four methods
## 1. Freq Poisson log-linear model ------
# (1 hour to run on 1k simulated datasets)
# Define functions to estimate N from a Poisson model
estimate_N <- function(model){
  coef <- coef(model)
  # Could change the covariate names below
  coef.sub <- coef[c("Region1", "Gender1", "Age1", "Region1:Gender1", "Region1:Age1", "Gender1:Age1")]
  coef.sub[is.na(coef.sub)] <- 0
  # Could change the index number below
  un1 <- exp(coef[1])
  un234 <- sum(exp(coef[1] + coef.sub[1:3]))
  un5 <- exp(coef[1] + sum(coef.sub[c(1, 2, 4)]))
  un6 <- exp(coef[1] + sum(coef.sub[c(1, 3, 5)]))   
  un7 <- exp(coef[1] + sum(coef.sub[c(2, 3, 6)]))
  un8<- exp(coef[1] + sum(coef.sub))
  unobs <- sum(un1, un234, un5, un6, un7, un8)
  Nhat <- sum(unobs, obs)
  return(Nhat)
}

estimate_poisson <- function(ds){
  scot.poi.start <- glm(Freq ~ S1 + S2 + S3 + S4, data = ds, family = poisson) 
  scot.poi <- glm(Freq ~ (S1 + S2 + S3 + S4 + Region + Gender + Age)^2, data = ds, family = poisson)
  NA.len <- length(which(is.na(ds$Freq))) 
  BICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi), direction = "forward", k = log(nrow(ds) - NA.len)) # Specify k = log(n) to use BIC (deleted 8 missing obs)
  obs <- sum(ds$Freq, na.rm = TRUE) # 5670
  Nhat_1 <- estimate_N(BICmodel.poi)

  # Non-parametric bootstrap: estimate the CI 
  # In practice, directly generate bootstrap samples from multinomial distribution with observed cell probabilities
  ##### Any other ways to do bootstrap? e.g. double bootstrap method
  nsims.b <- 1000
  cell.prob.b <- ds$Freq[!is.na(ds$Freq)]/obs
  # Simulation step
  set.seed(002)
  simdata.b <- as.matrix(rmultinom(nsims.b, obs, cell.prob.b), ncol = nsims.b)
  lst.b <- replicate(nsims.b, ds, simplify = FALSE)
  for (i in 1:nsims.b){
    lst.b[[i]]$Freq[!is.na(lst.b[[i]]$Freq)] <- simdata.b[, i]
  }
  fit.models <- lapply(lst.b, function(x){glm(BICmodel.poi$formula, family = poisson, data = x)})
  Nhat.b <- sapply(fit.models, estimate_N)
  # Percentiles of bootstrap estimates of N
  CI <- as.numeric(quantile(Nhat.b, c(0.025, 0.975)))
  
  # Test for Over-dispersion by Cameron & Trivedi (1990)
  library(AER)
  dispersion <- dispersiontest(BICmodel.poi, trafo = 1) 
  # Return the point estimate, CI & p-value for dispersion test
  results <- c(Nhat_1, CI, dispersion$p.value)
  return(results)
}

# Assess the performance of Poisson models
results.poi <- sapply(lst, estimate_poisson)
hist(results.poi[1, ])
abline(v = npop,col = "red")
hist(results.poi[4, ])
which(results.poi[4, ] < 0.05) # 0, simulated datasets based on this base data has no evidence of over-dispersion
# Compute the log of MSE for log of Nhat
logMSE.poi <- log(mean((log(poi.results[1, ]) - log(npop))^2)) #
# Coverage rate
cov_rate.poi <- length(which(results.poi[2, ] >= npop & results.poi[3, ] <= npop))/nsims
##### Try alternative def of coverage rate (wrt a dist on N - Paper 4)

##### Try to simulate (sparse & ) over-dispersed data and compare poisson & NB models, poisson & sparse poisson. Is there sparse NB?
##### How to compare SparseMSE with the other methods? Use data e.g. UKdat_5 in SparseMSE



## 2. Bayes GM model with no strata ------
# (? to run on 1k simulated datasets)


## 3. Bayes 
