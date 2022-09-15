## Parametric bootstrap: simulate datasets similar to the base data
# Use the base Poisson model
RM.poi.start <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D + S:E + A:D, data = RM_noyear, family = poisson) 
obs <- sum(RM_noyear$Freq)
# Simulate data based on the base Poisson model
nsims <- 500
unobs.pred <- predict.glm(RM.poi.start, RM_noyear_NA[is.na(RM_noyear_NA$Freq), ], type = "response") 
npop <- obs + sum(round(unobs.pred))
RM_noyear_fill <- RM_noyear_NA
RM_noyear_fill$Freq[is.na(RM_noyear_fill$Freq)] <- round(unobs.pred) 
RM.poi.fill <- glm(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D + S:E + A:D, data = RM_noyear_fill, family = poisson) 
cell.counts <- predict.glm(RM.poi.fill, RM_noyear_fill, type = "response")
# Calculate the probabilities of all capture histories
cell.prob <- cell.counts / sum(cell.counts)
# Simulation step
set.seed(124)
# Generate bootstrap samples from multinomial distribution 
realisations <- as.matrix(rmultinom(nsims, npop, cell.prob), ncol = nsims)
simdata <- realisations[-which(is.na(RM_noyear_NA$Freq)), ] 
lst <- replicate(nsims, RM_noyear, simplify = FALSE) 
for (i in 1:nsims){
  lst[[i]]$Freq <- simdata[, i]
}
# Check the sparsity for simulated datasets 
zero.cell <- sapply(lst, function(x){length(which(x$Freq == 0))}) 
small.cell <- sapply(lst, function(x){length(which(x$Freq <= 5 & x$Freq > 0))}) 
mean(zero.cell/nrow(RM_noyear))
mean(small.cell/nrow(RM_noyear))

## Simulate less sparse datasets based on the base data
set.seed(123)
add_random_number <- function(ds, p){ 
  n <- length(which(ds$Freq <= 5))
  select <- sample(which(ds$Freq <= 5), round(n/p), replace = F)  
  ds$Freq[select] <- ds$Freq[select] + sample(2:8, round(n/p), replace = TRUE)
  return(ds)
}
RM_noyear_NA_S <- add_random_number(RM_noyear_NA, 1.2) # Change the value of p
RM_noyear_S <- RM_noyear_NA_S[-which(is.na(RM_noyear_NA_S$Freq)), ]