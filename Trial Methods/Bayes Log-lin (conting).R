setwd("~/desktop/Dissertation/R Code Files")

## Data Setup
library(conting)
data("ScotPWID")
scot <- ScotPWID

## Bayesian log-linear model -----
# Set the seed for reproducibility of MCMC samples
set.seed(1)
# Set the maximal model in bict() to include all two-way interactions 
# Use the default Sabanes Bove and Held prior distribution
# Request 2 million MCMC iterations and save the MCMC output to external ﬁles every 1000 iterations
scot_ex <- bict(formula = y ~ (S1 + S2 + S3 + S4 + Age + Gender + Region) ^ 2, data = scot, 
                n.sample = 20000, null.move.prob = 0.5) # Use default a = , b = 
##### change the values of a and b in bict()
##### change the value of null.move.prob in bict()
scot_ex
list.files()
# Check the acceptance Rate
accept_rate(scot_ex)
# Informal Convergence Analysis of MCMC - Traceplot
ts.plot(scot_ex$BETA[, 1], ylab = "Intercept parameter", xlab = "Iteration number")
##### Perform thinning on the graph

# Use a burn-in phase of 200000 iterations (i.e., discard the ﬁrst 10% of the iterations) 
# Set thining at every 5th iteration due to the size of MCMC output 
# Use the Freeman-Tukey discrepancy statistic to assess model adequacy
scot_ex_bp <- bayespval(scot_ex, n.burnin = 200000, thin = 5, statistic = "FreemanTukey") # take two hours to run
scot_ex_bp
# Bayesian p value is 0.18 which indicates that there is no evidence of an inadequate model
##### Change the discrepancy statistic

# Calculate the posterior probabilities of the log-linear parameters in the maximal model and present those with probability greater than 0.70
inter_probs(scot_ex, n.burnin = 200000, thin = 5, cutoff = 0.70)

# Derive an MCMC sample from the posterior distribution of the total population size (taking account of different models) and find its posterior mean and 95% HPDI (Highest Posterior Density Interval)
scot_ex_tot <- total_pop(scot_ex, n.burnin = 200000, thin = 5)
scot_ex_tot # gives the posterior mean & 95% HPDI
# Alternatively, produce a different point estimate, e.g. posterior median 
round(median(scot_ex_tot$TOT))
# Or the one minimises relative squared error loss function 
round(mean(1/scot_ex_tot$TOT) / mean(1/(scot_ex_tot$TOT^2)), 0)
# Histogram of the posterior of N
hist(scot_ex_tot) # Bimodal due to the presence/absence of an interaction term (Conting Paper p23)
##### Produce estimates for each strata - see Y0 (sampled values of the missing and censored cell counts) returned from bict.fit
##### Produce the posterior distribution graph with multiple models as in Bayes GM
##### Model-averaged estimates of N?

##### Test for over-dispersion for Bayes Poisson log-linear model 




