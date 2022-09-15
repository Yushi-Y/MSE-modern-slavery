library(MCMCpack)
# install.packages("remotes")
# remotes::install_github("bernardsilverman/modslavmse")
library(modslavmse)

# Reproduce Silverman 2018 work
data(UKdat, UKdat_5, UKdat_4, Ned, Ned_5, NewOrl, NewOrl_5, Kosovo) # the datasets used in the paper
# make_AIC_stepwise_table1()  # Table 5
# make_AIC_stepwise_table2()  # Table 6
# make_allmodels_plots_script() # Figures 1, 2, 3 and 4
make_MCMCfit_tables_script() # Tables 7, 8, 9, 10, 11 and 12
make_MCMCeffects_table_script() #   Table 13
# make_madyor_table_script()     #  Table 14 and Figures 5 and 6, as well as some numbers in the text
# make_LCMCR_table_script()      # Table 15

bayes.poisson.out <- MCMCfit(UKdat, priorprec = 1, sigthresh = 2, totalpopest = T, 
                             mcmc = 100000, burnin = 5000, thin = 5, seed = 123) # Use lambda = 1, tau = 2 suggested by Sliverman(2020)
head(bayes.poisson.out)
plot(bayes.poisson.out) # Cannot show graph
summary(bayes.poisson.out)
# Nhat is 12363, CI is (9756, 12924) by the posterior median & CI of the (intercept)
##### Convergence analysis, e.g. trace and density plots for log-linear parameters