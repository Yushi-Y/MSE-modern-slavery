## Data Setup
library(LCMCR)
library(conting)
data("ScotPWID")
scot <- ScotPWID


## Data Handling -----
# Replace the table entries to 0 and 1 
scot <- data.frame(lapply(scot, as.character), stringsAsFactors = FALSE)
for (i in c("un", "GGC", "Male", "Young")){
  scot[scot == i] <- 0
}
for (i in c("obs", "Rest", "Female", "Old")){
  scot[scot == i] <- 1
}
# Marginalise the table over the covariates
colnames(scot)[1] <- "Freq"
scot <- data.frame(lapply(scot, as.numeric))
# Remove individuals not recorded on any list
Scot <- scot[apply(scot[,2:5], 1, sum) > 0, ]
# Repeat each capture history by the time of its frequencies (to make the table the right form)
Scot <- Scot[rep(seq_len(nrow(Scot)), times = Scot$Freq), -1]
lcmcr.Scot <- Scot[, 1:4]
# LCMCR requires that the list-membership columns be factors
for (col in c('S1', 'S2', 'S3', 'S4')) {
  lcmcr.Scot[[col]] <- as.factor(lcmcr.Scot[[col]])
}


## Bayesian estimation using non-parametric latent class models (no covariates) -----
# Set the seed for reproducibility (of MCMC samples)
set.seed(1)
sampler <- lcmCR(captures = lcmcr.Scot, tabular = FALSE, K = 16, a_alpha = 0.25, b_alpha = 0.25,
                 seed = "auto", buffer_size = 10000, thinning = 100) 
# Start with K* = 16, equal to the number of unique capture histories
# Thinning is set for the tracing buffer
##### Change values of a and b and check for sensitivity of estimates (if they change a lot)

# Inference
N <- lcmCR_PostSampl(sampler, burnin = 10000, samples = 5000, thinning = 100)
# par(mfrow = c(1, 1))
hist(N, breaks = 50,
     main = "Posterior Distribution of Total Population") # Bimodal. Wht???

qs <- quantile(N, c(0.025, 0.5, 0.975))
print(qs) # Posterior median is 25464
CI <- paste('[', round(qs['2.5%'], 0), ', ', round(qs['97.5%'], 0), ']', sep='')
print(CI) # 95% CI is (18446, 42435)
##### LCMCR uses MCMC to estimate the parameters, so need to apply convergence diagnostics to check for failure to converge