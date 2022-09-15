## Load the data
RM <- read.csv(file = "Romania 2015-2016.csv")[, -1]

## Data Handling for Romania data -----
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
# Convert all factors to numeric to fit other models
RM <- data.frame(lapply(RM, as.numeric))

# To compare across models, only allow 0 and 1 for exploitation type
RM_2fac <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D + Y, RM, FUN = sum)
# Replace all unobserved cells with NA 
NArows <- which(RM_2fac$R1==0 & RM_2fac$R2==0 & RM_2fac$R3==0 & RM_2fac$R4==0 & RM_2fac$R5==0)
RM_withNA <- RM_2fac
RM_withNA$Freq[NArows] <- NA
# Remove all NA cells
RM_2fac <- RM_withNA[-which(is.na(RM_withNA$Freq)), ]
# Separate data with NA by year
RM_withNA_15 <- RM_withNA[RM_withNA$Y == "2015", -10]
RM_withNA_16 <- RM_withNA[RM_withNA$Y == "2016", -10]
RM_withNA_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D, RM_withNA, FUN = sum, na.action = na.pass)
# Separate data by year
RM_2fac_15 <- RM_2fac[RM_2fac$Y == "2015", -10]
RM_2fac_16 <- RM_2fac[RM_2fac$Y == "2016", -10]

# Sum over years - still many zero counts
RM_noyear <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5 + S + A + E + D, RM_2fac, FUN = sum)
# Sum over all covariates - less zero counts but still quite sparse
RM_nocov <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac, FUN = sum)
# Prepare data with no covariates for each year
RM_nocov_15 <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac_15, FUN = sum)
RM_nocov_16 <- aggregate(Freq ~ R1 + R2 + R3 + R4 + R5, RM_2fac_16, FUN = sum)
