library(conting)
set.seed(111)
# Use the default Sabanes Bove and Held prior distribution with a = b = 10^(-3)
# Request 200k MCMC iterations
RM_withNA$Y <- as.factor(RM_withNA$Y)
ex <- bict(formula = Freq ~ (R1 + R2 + R3 + R4 + R5 + S + A + E + D + Y)^2, data = RM_withNA, 
           n.sample = 200000, null.move.prob = 0.5) 
# Check the acceptance Rate
accept_rate(ex)
# Convergence analysis of MCMC - traceplot
ts.plot(ex$BETA[, 1], ylab = "Intercept parameter", xlab = "Iteration number")
# Calculate posterior probabilities of the log-linear parameters  
# Present those with probability greater than 0.70
inter_probs(ex, n.burnin = 200000, thin = 50, cutoff = 0.70)

# Derive MCMC samples from the model-averaged posterior distribution of N
ex_tot <- total_pop(ex, n.burnin = 20000, thin = 50)
median(scot_ex_tot$TOT)
quantile(scot_ex_tot$TOT, c(0.025, 0.975))
hist(ex_tot)

# Use the Freeman-Tukey discrepancy statistic to assess model adequacy
# Set burn-in of 20000 iterations, set thinning at every 5th iteration 
ex_bp <- bayespval(ex, n.burnin = 20000, thin = 50, statistic = "FreemanTukey") 