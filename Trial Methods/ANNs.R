# ANNs
# devtools::install_github("bips-hb/neuralnet")
library(neuralnet)

## Simulate datasets similar to ScotPWID ----
## Data setup
library(conting)
library(mefa)
library(data.table)
data("ScotPWID")
scot <- ScotPWID # This is our base data

## Data Handling 
# Replace the table entries (characters) to 0 and 1 
scot <- data.frame(lapply(scot, as.character), stringsAsFactors = FALSE)
for (i in c("un", "GGC", "Male", "Young")){
  scot[scot == i] <- 0
}
for (i in c("obs", "Rest", "Female", "Old")){
  scot[scot == i] <- 1
}
str(scot)
Freq <- as.numeric(scot$y)
scot <- data.frame(lapply(scot[, -1], as.factor), Freq)


# Fit Poisson Log-linear Models on the base data (up to all two-way interactions) 
# Use stepwise forward regression to find the model with highest BIC as the generating model
scot.poi.start <- glm(Freq ~ S1 + S2 + S3 + S4, data = scot, family = poisson) # Start with the independence model (with no interactions)
scot.poi <- glm(Freq ~ (S1 + S2 + S3 + S4 + Region + Gender + Age)^2, data = scot, family = poisson)
BICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi), direction = "forward", k = log(nrow(scot) - 8)) # Specify k = log(n) to use BIC (deleted 8 missing obs)
coef <- coef(BICmodel.poi)

# Parametric bootstrap: simulate new datasets similar to the structure of the base data (here is scot)
# In practice, directly generate bootstrap samples from multinomial dist with observed cell probs
obs <- sum(scot$Freq, na.rm = TRUE)
nsims <- 1000
# scot[which(is.na(scot$Freq)), ] # Check all the NA cells
unobs.pred <- c(exp(coef[1]), exp(sum(coef[c(1, 8)])), exp(sum(coef[c(1, 7)])), exp(sum(coef[c(1, 7, 8, 13)])),
                exp(sum(coef[c(1, 6)])), exp(sum(coef[c(1, 6, 8, 12)])),
                exp(sum(coef[c(1, 6, 7)])), exp(sum(coef[c(1, 6:8, 12:13)]))) # Match the order with the order of NA cells
cell.counts <- numeric(2^(ncol(scot) - 1))
cell.counts[which(is.na(scot$Freq))] <- unobs.pred
cell.counts[which(!is.na(scot$Freq))] <- BICmodel.poi$fitted.values
npop <- obs + round(sum(unobs.pred))
# Calculate the probs of capture history
cell.prob <- cell.counts / sum(cell.counts)
# Simulation step
set.seed(111)
realisations <- as.matrix(rmultinom(nsims, npop, cell.prob), ncol = nsims)
npop_sim <- colSums(realisations) # The true total population sizes for simulated data samples
simdata <- realisations[which(!is.na(scot$Freq)), ] # Drop the rows that are supposed to have NA cells
lst <- replicate(nsims, scot, simplify = FALSE)
for (i in 1:nsims){
  lst[[i]]$Freq[!is.na(lst[[i]]$Freq)] <- simdata[, i]
}
# Check the sparsity for simulated datasets 
# NA.cell <- sapply(lst, function(x){length(which(is.na(x$Freq)))}) # Always 8 for all simulated datasets
zero.cell <- sapply(lst, function(x){length(which(x$Freq == 0))}) 
small.cell <- sapply(lst, function(x){length(which(x$Freq <= 10 & x$Freq > 0))}) 
hist(zero.cell)
hist(small.cell) # Both similar to the scot dataset

set.seed(123)
index <- sample(1:1000, 1) # 128
npop_mydata <- colSums(realisations)[index] # 32606
mydata <- lst[[index]]
mydata <- mydata[which(!is.na(mydata$Freq)), ]

# Split the data into training and test set
set.seed(123)
sample <- sample.split(mydata$Freq, SplitRatio = .80)
train.data <- subset(mydata, sample == TRUE)
test.data  <- subset(mydata, sample == FALSE)

## Fit an ANN
# Convert all the qualitative variables (factors) to binary ("dummy") variables
m <- model.matrix( 
  ~ Freq + S1 + S2 + S3 + S4 + Region + Gender + Age, 
  data = train.data
)
head(m)
set.seed(1111)
nn <- neuralnet( 
  Freq ~ S11 + S21 + S31 + S41 + Region1 + Gender1 + Age1, # Note that the variable names change in m
  data = m, hidden = c(2,2), threshold = 0.3, rep = 4, learningrate = 0.2, act.fct = "relu")
plot(nn, rep = "best")
# nn$result.matrix
m_test <- model.matrix( 
  ~ Freq + S1 + S2 + S3 + S4 + Region + Gender + Age, 
  data = test.data
)
pred <- predict(nn, m_test, rep = 4)
results <- data.frame(actual = test.data$Freq, prediction = pred)
results
attach(results)
deviation <- (actual-predicted)/actual
accuracy=1-abs(mean(deviation))


