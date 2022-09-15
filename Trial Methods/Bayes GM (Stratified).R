## Data setup
library(conting)
library(mefa)
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
colnames(scot)[1] <- "Freq"
scot <- data.frame(lapply(scot, as.numeric))
# Remove individuals not recorded on any list
Scot <- scot[apply(scot[,2:5], 1, sum) > 0, ]
# Repeat each capture history by the time of its frequencies (to make the table the right form)
Scot <- Scot[rep(seq_len(nrow(Scot)), times = Scot$Freq), -1]


## Bayesian Graphical Models (stratify by three covariates) -----
library(dga)
library(chron)
# Stratify by Region * Gender * Age - each strata is analysed separately (no info shared across stratas)
overlaps <- Scot[, 1:4]
obs <- sum(scot$Freq, na.rm = TRUE)
# Use the dates variable to stratify for Region
dates <- paste(
  sample(1:28, obs, replace = TRUE),"-",
  sample(1:4, obs, replace = TRUE),"-", Scot$Region
)
dates <- chron(dates, format = c(dates = "d-m-y"))

strata <- make.strata(overlaps, dates = dates, date.defs = "yearly", locations = Scot$Gender, demographics = Scot$Age)
rownames(strata$overlap.counts) <- c("GGC Male Young", "GGC Male Old", "GGC Female Young", "GGC Female Old", "Rest Male Young", "Rest Male Old", "Rest Female Young", "Rest Female Old")

# check to make sure that all strata are OK
# Lots of zero counts - 32 out of 120
check <- check.strata(strata)

# Look at the stratas, just to make sure
par(mfrow = c(4, 2), mar = rep(1, 4))
for (i in 1:nrow(strata$overlap.counts)) {
  venn4(strata$overlap.counts[i, ],
        main = rownames(strata$overlap.counts)[i],
        cex.main = 1
  )
}


# Load the graphs to make the estimates
data(graphs4)
# Select expansion factor defining the largest number of unrecorded elements
# This makes Nmissing <- 0:(sum(Y)*fac)
fac <- 20
# Set the prior
##### Change the prior: change delta (Delta is the hyper-parameter for the hyper-Dirichlet prior); 
##### Change the prior: change the arguments 'logprior', 'log.prior.model.weights' in bma.cr()
num.lists <- 4
delta <- 1 / 2 ^ num.lists

# Loop over stratas to calculate posterior distribution of the total population 
# Produce posterior plots of population in each strata
n <- nrow(strata$overlap.counts)
par(mfrow = c(2, 2), mar = rep(1.75, 4))
for (i in 1:(n/2)) {
  Nmissing <- 0:(sum(strata$overlap.counts[i, ]) * fac)
  Y <- array(strata$overlap.counts[i, ], dim = rep(2, num.lists))
  weights <- bma.cr(Y, Nmissing, delta, graphs4)
  plotPosteriorN(weights, Nmissing + sum(strata$overlap.counts[i, ]),
                 main = rownames(strata$overlap.counts)[i])
}

for (i in (n/2 + 1):n) {
  Nmissing <- 0:(sum(strata$overlap.counts[i, ]) * fac)
  Y <- array(strata$overlap.counts[i, ], dim = rep(2, num.lists))
  weights <- bma.cr(Y, Nmissing, delta, graphs4)
  plotPosteriorN(weights, Nmissing + sum(strata$overlap.counts[i, ]),
                 main = rownames(strata$overlap.counts)[i])
}
# Sum of posterior medians across 8 stratas is around 23000


# Estimate the population in each strata and sum them up
Nhat <- numeric(length = n)
CI <- matrix(nrow = 2, ncol = n)
for (i in 1:n) {
  Nmissing <- 0:(sum(strata$overlap.counts[i, ]) * fac)
  Y <- array(strata$overlap.counts[i, ], dim = rep(2, num.lists))
  weights <- bma.cr(Y, Nmissing, delta, graphs4)
  # Draw a sample of 1000 values from the estimated (averaged) posterior of each strata
  post_probs <- apply(weights, 2, sum) # Average over all possible models for each value of Nmissing
  m0s_sample_from_post <- sample(Nmissing, size = 1000, replace = TRUE, prob = post_probs)
  nhat_sample_from_post <- m0s_sample_from_post + sum(strata$overlap.counts[i, ]) # Add the observed counts
  # hist(nhat_sample_from_post, breaks = 50, main = "Posterior Distribution of the Population in Strata") # Histogram of posterior for population in the strata
  Nhat[i] <- quantile(nhat_sample_from_post, 0.5) # Here returns the posterior median of N
  CI[, i] <- quantile(nhat_sample_from_post, c(0.025, 0.975))
}
# Return the point estimate & CI of total population N
Nhat_sum <- sum(Nhat)
CI_sum <- rowSums(CI)
##### Compute the point estimate that minimises the relative squared error loss function (Paper 4)
##### How to get posterior model probabilities? - rows of 'weights'?
