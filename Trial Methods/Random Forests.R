## Bagged Trees & Random Forests
library(caret)
library(randomForest)
library(e1071)
library(caTools)
library(ModelMetrics)

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


## Fit RF on one of the simulated dataset ----
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

set.seed(123)
# Define the control - evaluate the model with a grid search of 10 folder
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")
# Fit the model on the training set
rf_default <- train(Freq ~ .,
                    data = train.data,
                    method = "rf",
                    metric = "RMSE",
                    trControl = trControl) 
print(rf_default) 
# The algorithm uses 500 trees and tested three different values of mtry: 2, 4, 7
# The final value used by the model was mtry = 4 with an RMSE of 75.98

# More detailed grid search for mtry
set.seed(1234)
tuneGrid <- expand.grid(.mtry = c(1:7))
rf_mtry <- train(Freq~.,
                 data = train.data,
                 method = "rf",
                 metric = "RMSE",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 5, # Minimum size of terminal nodes 
                 ntree = 300)
print(rf_mtry)
best_mtry <- rf_mtry$bestTune$mtry # 3
min(rf_mtry$results$RMSE) # 76.10

# Detailed grid search for maxnode give the best_mtry
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(5: 30)) {
  set.seed(1234)
  rf_maxnode <- train(Freq ~.,
                      data = train.data,
                      method = "rf",
                      metric = "RMSE",
                      tuneGrid = tuneGrid, # Uses best_mtry
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 5,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)
# maxnode = 25 gives the lowest mean RMSE

# Detailed grid search for maxnode give the best_mtry
store_maxtrees <- list()
for (ntree in c(200, 250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 1500, 2000, 2500, 3000)) {
  set.seed(5678)
  rf_maxtrees <- train(Freq ~.,
                       data = train.data,
                       method = "rf",
                       metric = "RMSE",
                       tuneGrid = tuneGrid, # Uses best_mtry
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 5,
                       maxnodes = 25,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
# ntree = 1500 gives the lowest mean RMSE and mean MAE

# Final model
fit_rf <- train(Freq ~.,
                data = train.data,
                method = "rf",
                metric = "RMSE",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                nodesize = 5,
                ntree = 1500,
                maxnodes = 25)
print(fit_rf)
plot(fit_rf$finalModel)

# Evaluate the fit - make predictions on the test data
prediction <-predict(fit_rf, test.data)
predict(fit_rf, c(0,0,0,0,1,1,1,0))
# Compute the RMSE for prediction 
rmse(predictions, test.data$Freq) # 27.23
rmse(predict(fit_rf, mydata), mydata$Freq) # 76.90

# Variable importance
varImp(fit_rf) 


# Compare with RMSE using Poisson log-linear models ----
# Model Selection using 'step' function
# Used forward selection (result slight different to 'both' or backward' regression)
scot.poi.start <- glm(Freq ~ S1 + S2 + S3 + S4, data = mydata, family = poisson) # Start with the independence model (with no interactions)
AICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi), direction = "forward") 
BICmodel.poi <- step(scot.poi.start, scope = formula(scot.poi), direction = "forward", k = log(nrow(mydata))) # Specify k = log(n) to use BIC
# summary(AICmodel.poi)
# summary(BICmodel.poi)
rmse(AICmodel.poi$fitted.values, mydata$Freq) # 3.81
rmse(BICmodel.poi$fitted.values, mydata$Freq) # 4.32
# Conclusion: RF gives a very bad fit for model