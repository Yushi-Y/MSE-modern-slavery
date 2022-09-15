## This is to test for over-dispersion in real datasets

## Load the Data
setwd("~/desktop/Dissertation/Data")
UK <- read.csv(file = "UK6.csv")
NL <- read.csv(file = "Netherlands.csv")
NO <- read.csv(file = "New Orleans.csv")
WE <- read.csv(file = "Western.csv")
KO <- read.csv(file = "Kosovo.csv")


### UK data----
head(UK)
## Frequentist: Poisson Log-linear Models 
# Include all two-way interactions
UK.poi <- glm(Freq ~ (LA + NG + PF + GO + GP + NCA)^2, data = UK, family = poisson)
# summary(UK.poi)
# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
UK.poi.start <- glm(Freq ~ LA + NG + PF + GO + GP + NCA, data = UK, family = poisson) # Start with the independence model (with no interactions)
AIC.UK.poi <- step(UK.poi.start, scope = formula(UK.poi), direction = "forward") 
BIC.UK.poi <- step(UK.poi.start, scope = formula(UK.poi), direction = "forward", k = log(nrow(UK))) # Specify k = log(n) to use BIC (deleted 8 missing obs)
summary(AIC.UK.poi)
summary(BIC.UK.poi) # Return the same model


## Test for Over-dispersion by Cameron & Trivedi (1990)
library(AER)
dispersiontest(AIC.UK.poi, trafo = 2) # p-value is 1, no evidence of over-dispersion
dispersiontest(BIC.UK.poi, trafo = 2) # p-value is 1, no evidence of over-dispersion
# Note that in summary(BICmodel.poi), the ratio of Residual deviance and residual degree of freedom is roughly 1 (1.07)
# This also suggests no evidence of over-dispersion (why? check formula in GLM notes)
# When there is over-dispersion, this means there is heterogeneity not fully captured by the covariates



### Netherlands data----
head(NL)
## Frequentist: Poisson Log-linear Models 
# Include all two-way interactions
NL.poi <- glm(Freq ~ (I + K + O + P + R + Z)^2, data = NL, family = poisson)
# summary(NL.poi)
# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
NL.poi.start <- glm(Freq ~ I + K + O + P + R + Z, data = NL, family = poisson) # Start with the independence model (with no interactions)
AIC.NL.poi <- step(NL.poi.start, scope = formula(NL.poi), direction = "forward") 
BIC.NL.poi <- step(NL.poi.start, scope = formula(NL.poi), direction = "forward", k = log(nrow(NL))) # Specify k = log(n) to use BIC (deleted 8 missing obs)
summary(AIC.NL.poi)
summary(BIC.NL.poi) # AIC returns two more insignificant pairwise interactions


## Test for Over-dispersion by Cameron & Trivedi (1990)
dispersiontest(AIC.NL.poi, trafo = 2) # p-value is 0.259, no evidence of over-dispersion
dispersiontest(BIC.NL.poi, trafo = 2) # p-value is 0.217, no evidence of over-dispersion



### New Orleans data----
head(NO)
## Frequentist: Poisson Log-linear Models 
# Include all two-way interactions
NO.poi <- glm(Freq ~ (A + B + C + D + E + F + G + H)^2, data = NO, family = poisson)
# summary(NO.poi)
# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
NO.poi.start <- glm(Freq ~ A + B + C + D + E + F + G + H, data = NO, family = poisson) # Start with the independence model (with no interactions)
AIC.NO.poi <- step(NO.poi.start, scope = formula(NO.poi), direction = "forward") 
BIC.NO.poi <- step(NO.poi.start, scope = formula(NO.poi), direction = "forward", k = log(nrow(NO))) # Specify k = log(n) to use BIC (deleted 8 missing obs)
summary(AIC.NO.poi)
summary(BIC.NO.poi) # AIC returns two more insignificant pairwise interactions


## Test for Over-dispersion by Cameron & Trivedi (1990)
dispersiontest(AIC.NO.poi, trafo = 2) # p-value is 1, no evidence of over-dispersion
dispersiontest(BIC.NO.poi, trafo = 2) # p-value is 1, no evidence of over-dispersion



### Western data----
head(WE)
## Frequentist: Poisson Log-linear Models 
# Include all two-way interactions
WE.poi <- glm(Freq ~ (A + B + C + D + E)^2, data = WE, family = poisson)
# summary(WE.poi)
# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
WE.poi.start <- glm(Freq ~ A + B + C + D + E, data = WE, family = poisson) # Start with the independence model (with no interactions)
AIC.WE.poi <- step(WE.poi.start, scope = formula(WE.poi), direction = "forward") 
BIC.WE.poi <- step(WE.poi.start, scope = formula(WE.poi), direction = "forward", k = log(nrow(WE))) # Specify k = log(n) to use BIC (deleted 8 missing obs)
summary(AIC.WE.poi)
summary(BIC.WE.poi) # Return the same model


## Test for Over-dispersion by Cameron & Trivedi (1990)
dispersiontest(AIC.WE.poi, trafo = 2) # p-value is 1, no evidence of over-dispersion
dispersiontest(BIC.WE.poi, trafo = 2) # p-value is 1, no evidence of over-dispersion


### Kosovo data----
head(KO)
## Frequentist: Poisson Log-linear Models 
# Include all two-way interactions
KO.poi <- glm(Freq ~ (EXH + ABA + OSCE + HRW)^2, data = KO, family = poisson)
# summary(WE.poi)
# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
KO.poi.start <- glm(Freq ~ EXH + ABA + OSCE + HRW, data = KO, family = poisson) # Start with the independence model (with no interactions)
AIC.KO.poi <- step(KO.poi.start, scope = formula(KO.poi), direction = "forward") 
BIC.KO.poi <- step(KO.poi.start, scope = formula(KO.poi), direction = "forward", k = log(nrow(KO))) # Specify k = log(n) to use BIC (deleted 8 missing obs)
summary(AIC.KO.poi)
summary(BIC.KO.poi) # Return the same model


## Test for Over-dispersion by Cameron & Trivedi (1990)
dispersiontest(AIC.KO.poi, trafo = 2) # p-value is 0.3669, there is no evidence of over-dispersion
dispersiontest(BIC.KO.poi, trafo = 2) # p-value is 0.3669, there is no evidence of over-dispersion



## Frequentist: Negative Binomial Models -----
# NB models are used for over-dispersed count data, that is when the conditional variance exceeds the conditional mean
# It is a generalization of Poisson regression with an extra parameter to model the over-dispersion
# Include all two-way interactions
library(MASS)
KO.nb <- glm.nb(Freq ~ (EXH + ABA + OSCE + HRW)^2, data = KO)
summary(KO.nb)

# Model Selection using 'step' function
# Used forward regression (result slight different to 'both' or backward' regression)
KO.nb.start <- glm.nb(Freq ~ EXH + ABA + OSCE + HRW, data = KO) # Start with the independence model 
AIC.KO.nb <- step(KO.nb.start, scope = formula(KO.nb), direction = "forward") 
BIC.KO.nb <- step(KO.nb.start, scope = formula(KO.nb), direction = "forward", k = log(nrow(KO))) # Specify k = log(n) to use BIC
summary(AIC.KO.nb)
summary(BIC.KO.nb) # Return the same model
# Here NB model has the same variables as Poisson model, with slightly different parameter estimates


# Compare RMSE of Poisson & NB models
library(ModelMetrics)
rmse(AIC.KO.poi) 
rmse(BIC.KO.poi) # 20.510
# sqrt(mean((KO$Freq-BIC.KO.poi$fitted.values)^2))
sqrt(mean((KO$Freq-AIC.KO.nb$fitted.values)^2))
sqrt(mean((KO$Freq-BIC.KO.nb$fitted.values)^2)) # 33.740, larger than Poisson?
# Visual inspection of model fit for both Poisson and NB models
library(countreg)
library(ggplot2)
library(cowplot)
# Generate the hanging rootograms (bars hanging from the red line) for the two models
# root.poi.AIC <- rootogram(AIC.KO.poi, style = "hanging", plot = FALSE)
# root.nb.AIC <- rootogram(AIC.KO.nb, style = "hanging", plot = FALSE)
root.poi.BIC <- rootogram(BIC.KO.poi, style = "hanging", plot = FALSE)
root.nb.BIC <- rootogram(BIC.KO.nb, style = "hanging", plot = FALSE)
ylims <- ylim(-2, 2)  # Common scale for comparison,adjust for your data
# plot_grid(autoplot(root.poi.AIC) + ylims, autoplot(root.nb.AIC) + ylims, ncol = 2, labels = "auto")
plot_grid(autoplot(root.poi.BIC) + ylims, autoplot(root.nb.BIC) + ylims, ncol = 2, labels = "auto")
# Interpretation: see https://www.r-bloggers.com/2016/06/rootograms/
# Expected counts - red line; observed counts - bars. If bars fall above 0 - over-predicts; below 0 - under-predicts
# Generally poor fit for both Poisson & NB, consecutive under-predictions
# Poisson & NB give very similar rootograms

# Given that counts are not sparse, we can also compare Poisson & NB using AIC and likelihood test
# AIC
AIC.KO.poi$aic # 161.3
AIC.KO.nb$aic # 154.7, so NB does slightly better in AIC
# LRT
library(pscl)
odTest(AIC.KO.nb, 0.05) # p-value = 0.002, so evidence of using NB over Poisson

