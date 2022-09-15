# Use both year's data and stratify by Y, S, E
library(dga)
library(chron)
RM_dga <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + E + Y, RM_2fac, sum)
RM_dga[RM_dga == "2015"] <- 0
RM_dga[RM_dga == "2016"] <- 1
# Repeat each capture history by the time of its frequencies 
RM_DGA <- RM_dga[rep(seq_len(nrow(RM_dga)), times = RM_dga$Freq), -9] # Remove the Freq column
# Stratify by Y*S*E - each strata is analysed separately 
overlaps <- RM_DGA[, 1:5]
obs <- sum(RM_dga$Freq)
# Use the dates variable to stratify for Y
dates <- paste(
  sample(1:28, obs, replace = TRUE),"-",
  sample(1:4, obs, replace = TRUE),"-", RM_DGA$Y)
dates <- chron(dates, format = c(dates = "d-m-y"))
strata <- make.strata(overlaps, dates = dates, date.defs = "yearly", locations = RM_DGA$S, demographics = RM_DGA$E)

data(graphs5)
# Nmissing <- 0:(sum(Y)*fac)
fac <- 10
# Set the prior
num.lists <- 5
delta <- 1 / 2 ^ num.lists # Default value of delta
# Estimate the population in each strata and sum them up
n <- nrow(strata$overlap.counts)
Nhat <- numeric(length = n)
CI <- matrix(nrow = 2, ncol = n)
for (i in 1:n) {
  Nmissing <- 0:(sum(strata$overlap.counts[i, ]) * fac)
  Y <- array(strata$overlap.counts[i, ], dim = rep(2, num.lists))
  weights <- bma.cr(Y, Nmissing, delta, graphs5)
  # Draw a sample of 10000 values from the averaged posterior of each strata
  post_probs <- apply(weights, 2, sum) 
  m0s_sample_from_post <- sample(Nmissing, size = 10000, replace = TRUE, prob = post_probs)
  nhat_sample_from_post <- m0s_sample_from_post + sum(Y) 
  Nhat[i] <- quantile(nhat_sample_from_post, 0.5) # Posterior median of N
  CI[, i] <- quantile(nhat_sample_from_post, c(0.025, 0.975))
}