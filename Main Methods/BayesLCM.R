library(LCMCR)
set.seed(111)
for (col in c('R1', 'R2', 'R3', 'R4', 'R5')) {
  RM_2fac_15[[col]] <- as.factor(RM_2fac_15[[col]])
}
# Split the dataset with 4 covariates to 16 separate datasets
ds_split <- split(RM_2fac_15[, c(1:5, 10)], rep(1:16, each = 31)) 
# Start with K* = 2^5 = 32, equal to the number of unique capture histories
# Use the default values of a = b = 0.25
# Set thinning for the tracing buffer
sampler <- lapply(ds_split, function(x) lcmCR(captures = x, tabular = TRUE, K = 32, a_alpha = 0.25, b_alpha = 0.25, seed = "auto", buffer_size = 10000, thinning = 100)) 
# Draw 500k samples with 50k burnins
N_samples <- lapply(sampler, function(y) lcmCR_PostSampl(y, burnin = 500000, samples = 5000, thinning = 100))
Nhat <- sum(sapply(N_samples, function(z) quantile(z, 0.5)))
CI <- rowSums(sapply(N_samples, function(z) quantile(z, c(0.025, 0.975)))) 