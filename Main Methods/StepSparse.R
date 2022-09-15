library(SparseMSE)
# Split the dataset to 16 separate tables with 4 covariates
RM15_split <- split(RM_2fac_15[, c(1:5, 10)], rep(1:16, each = 2^5 - 1)) 
# Convert dataframes to matrices for SparseMSE package
RM15_split <- lapply(RM15_split, function(x) data.matrix(x)) 
# Remove two stratified tables with zero counts in all cells
# Fit the model to each remaining stratified table
estimates <- lapply(RM15_split[-c(13,15)], function(x){
  estimatepopulation(x, nboot = 50, method = "stepwise", mX = NULL, pthresh = 0.02, iseed = 123, alpha = c(0.025, 0.975))}) 
# Use the default threshold p-value 0.02 and 100 bootstrap samples 

Nhat <- sum(Obs, sapply(estimates, function(x) x$popest)) 
CI <- rowSums(sapply(estimates, function(x) x$BCaquantiles), na.rm = TRUE) 
p.est <- sum(sapply(estimates, function(x) x$popest)[c(1:5, 7)])
CI <- CI + p.est 