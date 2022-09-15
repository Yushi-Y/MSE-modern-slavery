library(mefa)

## Load the data
setwd("~/desktop/Dissertation/Data")
RM <- read.csv(file = "Romania 2015-2016.csv")[, -1]
SB <- read.csv(file = "Serbia 2013-2016.csv")[, -1]


## Data Handling for Romania data -----
# Convert all categorical variables to factors to fit Poisson/NB models
RM_3fac <- RM
RM_3fac[, 6:10] <- data.frame(lapply(RM_3fac[, 6:10], as.factor))
str(RM_3fac)

# Separate data by each covariate
# rm(sum)
aggregate(Freq ~ S, RM_3fac, FUN = sum) 
aggregate(Freq ~ A, RM_3fac, FUN = sum) 
aggregate(Freq ~ A + Y, RM_3fac, FUN = sum) 
aggregate(Freq ~ E, RM_3fac, FUN = sum) 
aggregate(Freq ~ D, RM_3fac, FUN = sum) 
aggregate(Freq ~ D + Y, RM_3fac, FUN = sum) 

# Separate data by year 
aggregate(Freq ~ Y, RM_3fac, FUN = sum) 

RM_3fac_15 <- RM_3fac[RM_3fac$Y == "2015", ]
length(which(RM_3fac_15$Freq == 0)) 
length(which(RM_3fac_15$Freq <= 10 & RM_3fac_15$Freq > 0)) 

RM_3fac_16 <- RM_3fac[RM_3fac$Y == "2016", ]
length(which(RM_3fac_16$Freq == 0)) 
length(which(RM_3fac_16$Freq <= 10 & RM_3fac_16$Freq > 0)) 


# For two-factor treatment, replace all table entries to 0 and 1 
for (i in c("male", "adult", "other", "internal")){
  RM[RM == i] <- 0
}
for (i in c("beggary")){
  RM[RM == i] <- 0 
  # Treat 'beggary' & 'other' together as 'non-sexual' exploitation
}
for (i in c("female", "minor", "sexual", "external")){
  RM[RM == i] <- 1
}
str(RM)
# Convert all factors to numeric to fit other models
RM <- data.frame(lapply(RM, as.numeric))


# For comparison across models, only allow 0 and 1 for exploitation type
# As packages for Bayes threshold approach/Bayes log-linear model only work with tables of 0 and 1 
RM_2fac <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D + Y, RM, FUN = sum)
# Replace all unobserved cells with NA 
NArows <- which(RM_2fac$R1==0 & RM_2fac$R2==0 & RM_2fac$R3==0 & RM_2fac$R4==0 & RM_2fac$R5==0)
RM_withNA <- RM_2fac
RM_withNA$Freq[NArows] <- NA
# Remove all NA cells
RM_2fac <- RM_withNA[-which(is.na(RM_withNA$Freq)), ]
length(which(RM_2fac$Freq == 0)) 
length(which(RM_2fac$Freq <= 10 & RM_2fac$Freq > 0)) 
# Separate data with NA by year
RM_withNA_15 <- RM_withNA[RM_withNA$Y == "2015", -10]
RM_withNA_16 <- RM_withNA[RM_withNA$Y == "2016", -10]
RM_withNA_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D, 
                              RM_withNA, FUN = sum, na.action = na.pass)

# Separate data by year
RM_2fac_15 <- RM_2fac[RM_2fac$Y == "2015", -10]
length(which(RM_2fac_15$Freq == 0)) 
length(which(RM_2fac_15$Freq <= 10 & RM_2fac_15$Freq > 0)) 

RM_2fac_16 <- RM_2fac[RM_2fac$Y == "2016", -10]
length(which(RM_2fac_16$Freq == 0)) 
length(which(RM_2fac_16$Freq <= 10 & RM_2fac_16$Freq > 0)) 

# Sum over years - still many zero counts
RM_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D, RM_2fac, FUN = sum)
length(which(RM_noyear$Freq == 0)) 
length(which(RM_noyear$Freq <= 10 & RM_noyear$Freq > 0)) 

# Sum over all covariates - less zero counts but still quite sparse
RM_nocov <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac, FUN = sum)
length(which(RM_nocov$Freq == 0)) 
length(which(RM_nocov$Freq <= 10 & RM_nocov$Freq > 0)) 

# Prepare data with no covariates for each year
RM_nocov_15 <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac_15, FUN = sum)
RM_nocov_16 <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac_16, FUN = sum)




## Data Handling for Serbia data -----
# Replace all unobserved cells with NA 
# Convert all categorical variables to factors to fit Poisson/NB models
SB_3fac <- SB
SB_3fac[, 5:8] <- data.frame(lapply(SB_3fac[, 5:8], as.factor))
# Change covariate names 
levels(SB_3fac$S) <- c("female", "male")
levels(SB_3fac$A) <- c("adult", "minor")
levels(SB_3fac$N) <- c("non-Serbian", "Serbian")
levels(SB_3fac$E) <- c("beggary", "other", "sexual")
str(SB_3fac)

# Separate data by each covariate
aggregate(Freq ~ S, SB_3fac, FUN = sum) 
# 365 female, 509 male
aggregate(Freq ~ A, SB_3fac, FUN = sum) 
# 626 adult, 248 minor
aggregate(Freq ~ A + Y, SB_3fac, FUN = sum) 
aggregate(Freq ~ E, SB_3fac, FUN = sum) 
# 230 sexual, 589 other, 55 beggary
aggregate(Freq ~ N, SB_3fac, FUN = sum) 
# 776 Serbian, 98 non-Serbian
aggregate(Freq ~ N + Y, SB_3fac, FUN = sum) 
# Separate data by year 
aggregate(Freq ~ Y, SB_3fac, FUN = sum) 
# 156 obs for 2013, 480 obs for 2014
# 96 obs for 2015, 142 obs for 2016

SB_3fac_13 <- SB_3fac[SB_3fac$Y == "2013", ]
length(which(SB_3fac_13$Freq == 0)) 
# 346 out of 384 cells with zero counts
length(which(SB_3fac_13$Freq <= 10 & SB_3fac_13$Freq > 0)) 
# 34 out of 384 cells have counts less than 10

SB_3fac_14 <- SB_3fac[SB_3fac$Y == "2014", ]
length(which(SB_3fac_14$Freq == 0)) 
# 351 out of 384 cells with zero counts
length(which(SB_3fac_14$Freq <= 10 & SB_3fac_14$Freq > 0)) 
# 28 out of 384 cells have counts less than 10

SB_3fac_15 <- SB_3fac[SB_3fac$Y == "2015", ]
length(which(SB_3fac_15$Freq == 0)) 
# 352 out of 384 cells with zero counts
length(which(SB_3fac_15$Freq <= 10 & SB_3fac_15$Freq > 0)) 
# 31 out of 384 cells have counts less than 10

SB_3fac_16 <- SB_3fac[SB_3fac$Y == "2016", ]
length(which(SB_3fac_16$Freq == 0)) 
# 346 out of 384 cells with zero counts
length(which(SB_3fac_16$Freq <= 10 & SB_3fac_16$Freq > 0)) 
# 36 out of 384 cells have counts less than 10

# For two-factor treatment, replace all table entries to 0 and 1 
SB <- SB_3fac
SB[, 5:8] <- data.frame(lapply(SB[, 5:8], as.character))
for (i in c("male", "adult", "other", "Serbian")){
  SB[SB == i] <- 0
}
for (i in c("beggary")){
  SB[SB == i] <- 0 
  # Treat 'beggary' & 'other' together as 'non-sexual' exploitation
}
for (i in c("female", "minor", "sexual", "non-Serbian")){
  SB[SB == i] <- 1
}
str(SB)
# Convert all factors to numeric to fit other models
SB <- data.frame(lapply(SB, as.numeric))


# For comparison across models, only allow 0 and 1 for exploitation type
# As packages for Bayes threshold approach/Bayes log-linear model only work with tables of 0 and 1 
SB_2fac <- aggregate(Freq ~ R1 + R2 + R3 + R4 + S + A + N + E + Y, SB, FUN = sum)
# Remove year 2014 data
SB_2fac <- SB_2fac[-which(SB_2fac$Y == '2014'), ]
# Replace all unobserved cells with NA 
NArows <- which(SB_2fac$R1==0 & SB_2fac$R2==0 & SB_2fac$R3==0 & SB_2fac$R4==0)
SB_withNA <- SB_2fac
SB_withNA$Freq[NArows] <- NA
# Remove all NA cells
SB_2fac <- SB_withNA[-which(is.na(SB_withNA$Freq)), ]
length(which(SB_2fac$Freq == 0)) 
# 625 out of 720 cells have zero counts
length(which(SB_2fac$Freq <= 10 & SB_2fac$Freq > 0)) 
# 86 out of 720 cells have sparse counts
# Separate data with NA by year
SB_withNA_13 <- SB_withNA[SB_withNA$Y == "2013", -9]
SB_withNA_15 <- SB_withNA[SB_withNA$Y == "2015", -9]
SB_withNA_16 <- SB_withNA[SB_withNA$Y == "2016", -9]
SB_withNA_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + S + A + E + N, 
                              SB_withNA, FUN = sum, na.action = na.pass)


# Separate data by year
SB_2fac_13 <- SB_2fac[SB_2fac$Y == "2013", -9]
length(which(SB_2fac_13$Freq == 0)) 
# 206 out of 236 cells with zero counts
length(which(SB_2fac_13$Freq <= 10 & SB_2fac_13$Freq > 0)) 
# 26 out of 236 cells have counts less than 10

SB_2fac_15 <- SB_2fac[SB_2fac$Y == "2015", -9]
length(which(SB_2fac_15$Freq == 0)) 
# 206 out of 236 cells with zero counts
length(which(SB_2fac_15$Freq <= 10 & SB_2fac_15$Freq > 0)) 
# 28 out of 236 cells have counts less than 10

SB_2fac_16 <- SB_2fac[SB_2fac$Y == "2016", -9]
length(which(SB_2fac_16$Freq == 0)) 
# 206 out of 238 cells with zero counts
length(which(SB_2fac_16$Freq <= 10 & SB_2fac_16$Freq > 0)) 
# 29 out of 238 cells have counts less than 10

# Sum over years - still many zero counts
SB_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + S + A + N + E, SB_2fac, FUN = sum)
length(which(SB_noyear$Freq == 0)) 
# 188 out of 240 cells with zero counts
length(which(SB_noyear$Freq <= 10 & SB_noyear$Freq > 0)) 
# 42 out of 240 cells have counts less than 10

# Sum over all covariates - less zero counts but still quite sparse
SB_nocov <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac, FUN = sum)
length(which(SB_nocov$Freq == 0)) 
# 5 out of 15 cells with zero counts
length(which(SB_nocov$Freq <= 10 & SB_nocov$Freq > 0)) 
# 6 out of 15 cells have counts less than 10

# Prepare data with no covariates for each year
SB_nocov_13 <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac_13, FUN = sum)
SB_nocov_15 <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac_15, FUN = sum)
SB_nocov_16 <- aggregate(Freq ~ R1 + R2 + R3 + R4, SB_2fac_16, FUN = sum)


# Repeat each capture history by the time of its frequencies
# RM_spread <- RM[rep(seq_len(nrow(RM)), times = RM$Freq), -1]

# Sort dataframe in lexigraphical order
# order(..., na.last = TRUE, decreasing = FALSE)

# Try Poisson & NB fitting with three levels for reference
