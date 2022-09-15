## Data setup
library(conting)
library(mefa)
library(data.table)
library(coxed)
data("ScotPWID")
scot <- ScotPWID

# Data Handling -----
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
overlaps <- Scot[, 1:4]
# Count the number of same rows in the dataframe
scot.list <- as.data.frame(setDT(overlaps)[, list(Freq = .N), names(overlaps)])
str(scot.list)

# Examine the data sparsity: count the number of NA, 0, non-zero small counts in 'Freq'
length(which(scot.list$Freq == 0)) # 0 out of 16 cells with zero counts
length(which(scot.list$Freq <= 10 & scot.list$Freq > 0)) # 4 out of 16 cells have counts less than 10
# The data is no longer sparse after marginalisation - not really suitable to use in SparseMSE


## Frequentist: Poisson Log-linear Models - now handles sparsity -----
# Consider all two-way interactions
library(SparseMSE)
# Note that here scot.list data is not really sparse
estimates <- estimatepopulation(scot.list, nboot = 1000, method = "stepwise", mX = NULL, pthresh = 0.02, 
                   iseed = 1234, alpha = c(0.025, 0.975))
estimates$popest
estimates$BCaquantiles[1]
# nboot defines the number of boostraps for BCabootstrap CI
# mX = NULL means starting with main effect model
# alpha defines the quantiles of BCa bootstrap CI that we want to return
# Point estimate: 34828.35; 95% CI:(28151.44, 43821.36) - slightly higher estimates than Poisson & narrower CI
# This gives the model selected by the stepwise approach, parameter estimates, AIC, point estimate & BCa bootstrap CI

##### Use different pthresh values (threshold p-value) based on e.g. simulation (Sparsity Paper S4.3)