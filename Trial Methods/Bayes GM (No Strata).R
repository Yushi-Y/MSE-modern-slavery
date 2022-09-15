## Data setup
library(conting)
library(mefa)
library(plyr)
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
colnames(scot)[1] <- "Freq"
scot <- data.frame(lapply(scot, as.numeric))
# Remove individuals not recorded on any list
Scot <- scot[apply(scot[,2:5], 1, sum) > 0, ]
# Repeat each capture history by the time of its frequencies (to make the table the right form)
Scot <- Scot[rep(seq_len(nrow(Scot)), times = Scot$Freq), -1]


## Bayesian Graphical Models (stratify by one covariate) -----
library(dga)
overlaps <- Scot[, 1:4]
obs <- sum(scot$Freq, na.rm = TRUE)

# Stratify by ONE covariate only, say Region (can also by Gender or Age)
strata <- make.strata(overlaps, locations = Scot$Region)
rownames(strata$overlap.counts) <- c("GGC", "Rest")
# check to make sure that all strata are OK
check <- check.strata(strata)
# Look at the stratas, just to make sure
par(mfrow = c(1, 2), mar = rep(1, 4))
for (i in 1:nrow(strata$overlap.counts)) {
  venn4(strata$overlap.counts[i, ],
        main = rownames(strata$overlap.counts)[i],
        cex.main = 1
  )
}

# Load the graphs to make the estimates
data(graphs4)
# Select expansion factor defining the largest number of unrecorded elements
# So that max(Nmissing) = sum(Y)*fac
fac <- 8

# Set the prior
##### Change the value of delta?
##### Change the logprior (on Nmissing) in bma.cr()?

num.lists <- 4
delta <- 1 / 2 ^ num.lists

# Loop over stratas to calculate posterior distributions of the total population size for each stratum
par(mfrow = c(1, 2), mar = rep(1.75, 4))
n <- nrow(strata$overlap.counts)
for (i in 1:n) {
  Nmissing <- 0:(sum(strata$overlap.counts[i, ]) * fac)
  Y <- array(strata$overlap.counts[i, ], dim = rep(2, num.lists))
  weights <- bma.cr(Y, Nmissing, delta, graphs4)
  plotPosteriorN(weights, Nmissing + sum(strata$overlap.counts[i, ]),
                 main = rownames(strata$overlap.counts)[i])
}
##### Sum of posterior medians across 2 stratas is around 25000


# Try no stratification -----
counts <- count(overlaps)
counts[nrow(counts) + 1,] <- c(1, 1, 1, 1, 0)
counts <- rbind(c(0, 0, 0, 0, 0), counts) # Replace the NA counts to 0
Y <- array(counts$freq, dim=c(2,2,2,2))  # The dataframe must be sorted by the S1, S2, S3, S4 columns so that this array is shaped properly

Nmissing <- 1:(10 * sum(Y))        # Nmissing is a sequence from 1 to (say) ten times the number of observed values.
num.lists <- 4                  
delta <- 1/(2 ^ num.lists) # default value of delta
data(graphs4)
weights <- bma.cr(Y, Nmissing, delta, graphs4)
par(mfrow = c(1, 1), mar = rep(1.75, 4))
plotPosteriorN(weights, (Nmissing + sum(Y)))
##### The posterior median is around 32000


-----# haven't look into this ------
# Draw a sample of 1000 values from the estimated (averaged) posterior of the total population
post_probs <- apply(weights, 2, sum) # Average over all possible models for each value of Nmissing
m0s_sample_from_post <- sample(Nmissing, size = 1000, replace = TRUE, prob = post_probs)
nhat_sample_from_post <- m0s_sample_from_post + sum(Y) # Add the observed value

hist(nhat_sample_from_post, breaks = 50,
     main = "Posterior Distribution of Estimated Killings",
     xlim=c(1500, 4500))
# This is very similar to the graph before

qs <- quantile(nhat_sample_from_post, c(0.025, 0.5, 0.975))
CI <- paste('[', round(qs['2.5%'], 0), ', ', round(qs['97.5%'], 0), ']', sep='') # (2066, 2758)
##### Compute the Bayes estimate that minimises the relative squared error loss function (Paper 4)
##### How to get posterior model probablities?
