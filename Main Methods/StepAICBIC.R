## Take AIC model as an example
# Fit on both years' data and treat 'Y' as a dummy variable
library(coxed)
RM_2fac$Y <- as.factor(RM_2fac$Y)
RM_withNA$Y <- as.factor(RM_withNA$Y)
RM.poi.full <- glm(Freq ~ (R1 + R2 + R3 + R4 + R5 + S + A + E + D + Y)^2, data = RM_2fac, family = poisson)
RM.poi.start <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D + Y, data = RM_2fac, family = poisson) 
AIC.poi <- step(RM.poi.start, scope = formula(RM.poi.full), direction = "forward") 
summary(AIC.poi) 
# Check for structural zeros
sum(RM_2fac$Freq == 0) 
sum(round(AIC.poi$fitted.values) == 0) 

# Total population point estimate 
Obs <- sum(RM_2fac$Freq) 
Obs_15 <- sum(RM_2fac_15$Freq)
Obs_16 <- sum(RM_2fac_16$Freq)
Unobs_AIC <- predict(AIC.poi, RM_withNA[is.na(RM_withNA$Freq), ], type = "response") 
Unobs_AIC_15 <- predict(AIC.poi, RM_withNA[(is.na(RM_withNA$Freq) & RM_withNA$Y == "2015"), ], type = "response") 
Unobs_AIC_16 <- predict(AIC.poi, RM_withNA[(is.na(RM_withNA$Freq) & RM_withNA$Y == "2016"), ], type = "response") 
Nhat_AIC <- sum(Unobs_AIC, Obs) 
Nhat_AIC_15 <- sum(Unobs_AIC_15, Obs_15) 
Nhat_AIC_16 <- sum(Unobs_AIC_16, Obs_16) 

# Estimate the CI: non-parametric bootstrap
nsims <- 50
cell.prob <- RM_2fac$Freq/Obs
set.seed(111)
realisations <- rmultinom(nsims, Obs, cell.prob)
lst <- replicate(nsims, RM_2fac, simplify = FALSE)
for (i in 1:nsims){
  lst[[i]]$Freq <- realisations[, i]
}
fit.models_AIC <- lapply(lst, function(x){glm(AIC.poi$formula, family = poisson, data = x)})
Nhat.bootstrap_AIC <- sapply(fit.models_AIC, function(x){sum(Obs, predict.glm(x, RM_withNA[is.na(RM_withNA$Freq), ], type = "response"))})
# Estimate the CI using percentiles of bootstrap estimates of N
quantile(Nhat.bootstrap_AIC, c(0.025, 0.975)) 
# Alternatively, calculate the BCa bootstrap CI
CI_AIC <- bca(Nhat.bootstrap_AIC, conf.level = 0.95) 

Nhat.bootstrap_AIC_15 <- sapply(fit.models_AIC, function(x){sum(Obs_15, predict.glm(x, RM_withNA[(is.na(RM_withNA$Freq) & RM_withNA$Y == "2015"), ], type = "response"))})
quantile(Nhat.bootstrap_AIC_15, c(0.025, 0.975)) 
CI_AIC_15 <- bca(Nhat.bootstrap_AIC_15, conf.level = 0.95) 

Nhat.bootstrap_AIC_16 <- sapply(fit.models_AIC, function(x){sum(Obs_16, predict.glm(x, RM_withNA[(is.na(RM_withNA$Freq) & RM_withNA$Y == "2016"), ], type = "response"))})
quantile(Nhat.bootstrap_AIC_16, c(0.025, 0.975)) 
CI_AIC_16 <- bca(Nhat.bootstrap_AIC_16, conf.level = 0.95) 