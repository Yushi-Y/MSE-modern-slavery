## Data setup
library(conting)
library(mefa)
library(rsq)
data("ScotPWID")
scot <- ScotPWID

## Data Handling -----
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
# Examine the data sparsity: count the number of NA, 0, non-zero small counts in 'Freq'
length(which(is.na(scot$Freq))) # 8 out of 128 unknown cells
length(which(scot$Freq == 0)) # 32 out of 128 cells with zero counts
length(which(scot$Freq <= 10 & scot$Freq > 0)) # 47 out of 128 cells have counts less than 10
##### Produce histograms???
       

## Frequentist: Poisson Log-linear Models -----
# Include all two-way interactions
# Use glm function
scot.poi.full <- glm(Freq ~ (S1 + S2 + S3 + S4 + Region + Gender + Age)^2, data = scot, family = poisson)
summary(scot.poi)

# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
scot.poi.start <- glm(Freq ~ S1 + S2 + S3 + S4 + Region + Gender + Age, data = scot, family = poisson) # Start with the independence model (with no interactions)
NA.len <- length(which(is.na(scot$Freq))) 
AICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi.full), direction = "forward") 
BICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi.full), direction = "forward", k = log(nrow(scot) - NA.len)) # Specify k = log(n) to use BIC (deleted 8 missing obs)
# BICmodel <- step(scot.glm, direction = "backward", k = log(nrow(Scot) - 8))
summary(AICmodel.poi)
summary(BICmodel.poi) # Parameter estimates
# confint(BICmodel) # Inaccurate CI of MLEs with small counts b\c normal approximation does not hold (GLM notes) - use parametric bootstrap 
rsq.kl(BICmodel.poi) # 0.994, very high. Is it good to use?
# Count the number of zeros predicted by the BIC model - very close to the real zeros - no evidence of structural zeros
# Hence no need for Zero-inflated model
sum(round(BICmodel.poi$fitted.values) == 0) # 31, true value is 32
# Could use bootstrap to do chi-squared goodness-of-fit tests based on deviance for the top 5 BIC models (Paper 3) because inaccurate p-values for small Poisson counts - see later


# Total Population Point Estimate -----
Obs <- sum(scot$Freq, na.rm = TRUE) # 5670
Unobs <- predict.glm(BICmodel.poi, scot[is.na(scot$Freq), ], type = "response") # A vector of unobserved counts
Nhat <- sum(Unobs, Obs) # 32605.86


# Non-parametric bootstrap: estimate the CI -----
# In practice, directly generate bootstrap samples from multinomial dist with observed cell probs
library(coxed)
nsims <- 1000
cell.prob <- scot$Freq[!is.na(scot$Freq)]/Obs
set.seed(123)
realisations <- rmultinom(nsims, Obs, cell.prob)
lst <- replicate(nsims, scot, simplify = FALSE)
for (i in 1:nsims){
  lst[[i]]$Freq[!is.na(lst[[i]]$Freq)] <- realisations[, i]
}
fit.models <- lapply(lst, function(x){glm(BICmodel.poi$formula, family = poisson, data = x)})
Nhat.bootstrap <- sapply(fit.models, function(x){sum(obs, predict.glm(x, scot[is.na(scot$Freq), ], type = "response"))})
# Estimate the CI using percentiles of bootstrap estimates of N
quantile(Nhat.bootstrap, c(0.025, 0.975)) # (28027.73, 38761.59)
# Alterntivly, calculate the BCa bootstrap CI
bca(Nhat.bootstrap, conf.level = 0.95) # (28746.56.57, 39241.60), narrower than above


# Parametric bootstrap for goodness-of-fit: estimate the chi-square dist -----
# Simulation random samples based on the stepwise forward BIC model
# In practice, directly generate bootstrap samples from multinomial dist with observed cell probs
obs <- sum(scot$Freq, na.rm = TRUE)
nsims <- 1000
unobs.pred <- predict.glm(BICmodel.poi, scot[is.na(scot$Freq), ], type = "response") 
cell.counts.p <- predict.glm(BICmodel.poi, scot, type = "response") 
npop <- obs + round(sum(unobs.pred))
# Calculate the probabilities of capture history
cell.prob.p <- cell.counts.p / sum(cell.counts.p)
# Simulation step
set.seed(101)
realisations.p <- as.matrix(rmultinom(nsims, npop, cell.prob.p), ncol = nsims)
simdata <- realisations.p[which(!is.na(scot$Freq)), ] # Drop the rows that are supposed to have NA cells
lst.p <- replicate(nsims, scot, simplify = FALSE)
for (i in 1:nsims){
  lst.p[[i]]$Freq[!is.na(lst.p[[i]]$Freq)] <- simdata[, i]
}
# Fit the BIC model to the simulated datasets
fit.models.p <- lapply(lst.p, function(x){glm(BICmodel.poi$formula, family = poisson, data = x)})
deviance.bootstrap <- sapply(fit.models.p, function(x){x$deviance})
1 - ecdf(deviance.bootstrap)(BICmodel.poi$deviance) # 0.111, suggest good fit
# Note the 5% p-value criteria could be relaxed in light of large n 



# Alternatively, use Rcapture package to fit log-linear models ----
# Take too long to run b\c many possible models (with two-way interactions) for 7 lists
# library(Rcapture)
# Scot <- data.frame(lapply(scot, as.character))
# Scot <- data.frame(lapply(Scot, as.numeric))
# result <- closedpMS.t(Scot, dfreq = TRUE, h = 'Poisson', maxorder = 2, stopiflong = FALSE)
# print(result)
# models <- as.data.frame(result$results)
# models <- models[order(models$BIC), ]



# Lasso Regression for Variable Selection ----
library(glmnet)
# First step: using .^2 to include all pairwise interactions
f <- as.formula(y ~ .^2)
y <- scot[is.na(scot$Freq) == FALSE, 8]
# Second step: using 'model.matrix' function to take advantage of f
x <- model.matrix(f, scot[is.na(scot$Freq) == FALSE, -8])[, -1]
# Fit the model using cross-validation to select the optimal lambda (the default loss is Poisson deviance)
cvfit <- cv.glmnet(x, y, family = "poisson") 
plot(cvfit)
cvfit$lambda.min # The lambda that gives the minimum corss-validated error
coef.lasso <- coef(cvfit, s = "lambda.min") # s is the penalty parameter lambda
coef.lasso # Under the cross-validated model, only one variable is dropped - the focus for lasso is prediction, not which variables included
# Unobs <- exp(coef[1]) + sum(exp(coef[1] + coef[6:8])) + exp(sum(coef[c(1, 6, 8, 12)])) +
  # exp(sum(coef[c(1, 7, 8, 13)])) + exp(sum(coef[c(1, 6, 7)])) + 
  # exp(sum(coef[c(1, 6:8, 12:13)])) # 26935.86 
# Nhat <- sum(Unobs, Obs) # 32605.86
##### No standard way to compute CI for lasso - bootstrap for CI is not good as each bootstrap sample could keep a different set of variables 
  
  
  
## Test for Over-dispersion by Cameron & Trivedi (1990) -----
library(AER)
dispersiontest(BICmodel.poi, trafo = 2)  
# trafo = 1 and trafo = 2 yield the linear and quadratic formulations of f
# trafo = 2 corresponds to a negative binomial (NB) model with quadratic variance function (called NB2 by Cameron and Trivedi, 2005)
# trafo = 1 corresponds to a NB model with linear variance function (called NB1 by Cameron and Trivedi, 2005) or quasi-Poisson model with dispersion parameter
# p-value is 0.78, no evidence of over-dispersion
# Note that in summary(BICmodel.poi), the ratio of Residual deviance and residual degree of freedom is roughly 1 (1.05)
# This also suggests no evidence of over-dispersion (why? check formula in GLM notes)
# When there is over-dispersion, this means there is heterogeneity not fully captured by the covariates

## Frequentist: Negative Binomial Models -----
# NB models are used for over-dispersed count data, that is when the conditional variance exceeds the conditional mean
# It is a generalization of Poisson regression with an extra parameter to model the over-dispersion
# Include all two-way interactions
library(MASS)
scot.nb.full <- glm.nb(Freq ~ (S1 + S2 + S3 + S4 + Region + Gender + Age)^2, data = scot)
summary(scot.nb)

# Model Selection using 'step' function, used forward regression
scot.nb.start <- glm.nb(Freq ~ S1 + S2 + S3 + S4 + Region + Gender + Age, data = scot) # Start with the independence model 
AICmodel.nb <- step(scot.nb.start, scope = formula(scot.nb.full), direction = "forward") 
BICmodel.nb <- step(scot.nb.start, scope = formula(scot.nb.full), direction = "forward", k = log(nrow(scot) - NA.len)) # Specify k = log(n) to use BIC (deleted 8 missing obs)
# BICmodel <- step(scot.glm, direction = "backward", k = log(nrow(Scot) - 8))
summary(AICmodel.nb)
summary(BICmodel.nb) # Here NB model has the same variables as Poisson model, with very similar parameter estimates

# Total Population Point Estimate -----
Obs <- sum(scot$Freq, na.rm = TRUE) # 5670
Unobs <- predict.glm(BICmodel.nb, scot[is.na(scot$Freq), ], type = "response") # A vector of unobserved counts
Nhat <- sum(Unobs, Obs) # 32605.86, same as Poisson model as no over-dispersion


# Non-parametric bootstrap: estimate the CI -----
# In practice, directly generate bootstrap samples from multinomial dist with observed cell probs
library(coxed)
nsims <- 1000
cell.prob <- scot$Freq[!is.na(scot$Freq)]/Obs
set.seed(123)
realisations <- rmultinom(nsims, Obs, cell.prob)
lst <- replicate(nsims, scot, simplify = FALSE)
for (i in 1:nsims){
  lst[[i]]$Freq[!is.na(lst[[i]]$Freq)] <- realisations[, i]
}
fit.models <- lapply(lst, function(x){glm.nb(formula(ds.nb.BIC), data = x)})
Nhat.bootstrap <- sapply(fit.models, function(x){sum(obs, predict.glm(x, scot[is.na(scot$Freq), ], type = "response"))})
# Estimate the CI using percentiles of bootstrap estimates of N
quantile(Nhat.bootstrap, c(0.025, 0.975)) # (28027.73, 38761.59)
# Alterntivly, calculate the BCa bootstrap CI
bca(Nhat.bootstrap, conf.level = 0.95) # (28746.56.57, 39241.60), narrower than above


# Compare RMSE of Poisson & NB models
library(ModelMetrics)
rmse(AICmodel.poi)
rmse(BICmodel.poi)
# Manual calculation of RMSE
sqrt(mean((scot$Freq-AICmodel.nb$fitted.values)^2))
sqrt(mean((scot$Freq-BICmodel.nb$fitted.values)^2))

# Visual inspection of model fit for both Poisson and NB models
# install.packages("countreg", repos = "http://R-Forge.R-project.org")
library(countreg)
library(ggplot2)
library(cowplot)
# Generate the hanging rootograms (bars hanging from the red line) for the two models
# root.poi.AIC <- rootogram(AICmodel.poi, style = "hanging", plot = FALSE)
# root.nb.AIC <- rootogram(AICmodel.nb, style = "hanging", plot = FALSE)
root.poi.BIC <- rootogram(BICmodel.poi, style = "hanging", plot = FALSE)
root.nb.BIC <- rootogram(BICmodel.nb, style = "hanging", plot = FALSE)
ylims <- ylim(-2, -7)  # Common scale for comparison, adjust for your data
# plot_grid(autoplot(root.poi.AIC) + ylims, autoplot(root.nb.AIC) + ylims, ncol = 2, labels = "auto")
plot_grid(autoplot(root.poi.BIC) + ylims, autoplot(root.nb.BIC) + ylims, ncol = 2, labels = "auto")
# Interpretation: see https://www.r-bloggers.com/2016/06/rootograms/
# Expected counts - red line; observed counts - bars. If bars fall above 0 - over-predicts; below 0 - under-predicts
# Generally poor fit, althogh zero count is nicely estimated (good agreement), a large number of consecutive under-predictions, especially for larger counts
# Use AIC and BIC models give (almost) the same rootograms

----------------------------
## Data processing: add 'Year' in data -----
# Why not include Year (with at least three levels) as a variable? 
# In Davina's method, used Year as polynomial variables not indicator variables - not compatible with other methods
# Also in Recapture package, only 0 and 1 are allowed for table entries (i.e. only two-level indicator variables)

# Modify the dataset to include 'Year'
myScot <- rep(scot, each = 3)
# Create the new columns of 'Year' & 'Freq'
myScot$Year <- rep(c(2006, 2007, 2008), nrow(scot))
myScot$Year <- as.factor(myScot$Year)
myScot$Freq <- NA
for (i in 1:nrow(scot)){
  if (is.na(scot$y[i]) == FALSE){
    myScot$Freq[c(3*i - 2, 3*i - 1, 3*i)] <- rmultinom(1, scot$y[i], c(0.23, 0.33, 0.44))
  }
}
myScot <- myScot[, -1]
str(myScot)

# The final data with all covariates as factors (despite 'Freq')
scot_y <- data.frame(lapply(myScot[, -9], as.factor))
scot_y$Freq <- as.numeric(myScot$Freq)
str(scot_y)


## Frequentist: Poisson Log-linear Models
# Include all two-way interactions
# Use the glm function
scot_y.glm <- glm(Freq ~ (S1 + S2 + S3 + S4 + Region + Gender + Age + Year)^2, data = scot_y, family = poisson)
summary(scot.glm)
# Specify k = log(n) in the 'step' function to get BIC
# I use both forward & backward regression (result the same with only backward)
BICmodel_y <- stepAIC(scot_y.glm, direction = "both", k = log(n))
summary(BICmodel_y)
confint(BICmodel_y) # Inaccurate CI of MLEs with small counts (GLM notes) - need bootstrap 
rsq.kl(BICmodel_y)



