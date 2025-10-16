### Script for generating mean, SD, and association stats with PD for each biomarker
## ==== Set up tables and environment ====
library(data.table)

# Load data
bio <- fread("marker.txt") # marker.txt contains the sample FID/IID and biomarker measure for each person for the given biomarker
colnames(bio) <- c("FID","IID","bio")
cov <- fread("covar.txt") # basic covar file

# Format data
merged <- merge(bio,cov,by="FID",all.x=TRUE)
merged <- merged[!is.na(merged$bio),]
merged <- merged[!is.na(merged$Status1),]
merged$bio <- as.numeric(merged$bio)
merged$Status2 <- as.numeric(merged$Status2) # Status2 marks SWEDDs with the value "3"
merged$Status1 <- as.numeric(merged$Status1) # Status1 marks prodromals with the value "3"
merged$Age <- as.numeric(merged$Age)
merged$Sex <- as.numeric(merged$Sex)
merged$PC1 <- as.numeric(merged$PC1)
merged$PC2 <- as.numeric(merged$PC2)
merged$PC3 <- as.numeric(merged$PC3)
merged$PC4 <- as.numeric(merged$PC4)
merged$PC5 <- as.numeric(merged$PC5)

# Recode disease status and remove SWEDD
merged$Status2[merged$Status2 == 1] <- 0
merged$Status2[merged$Status2 == 2] <- 1
merged <- merged[merged$Status2 != 3,]

## ==== Get mean, median and SD ====
mean <- mean(merged$bio)
median <- median(merged$bio)
sd <- sd(merged$bio)
mean_median_sd <- cbind(mean, median, sd)

# Save
write.table(mean_median_sd, file = "means.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## ==== Run regression of biomarker x status ==== 
# basic covs, prodromal included, prodromal removed
# prodromal included with controls
fit_basic <-  glm(Status2 ~ bio, data = merged, family = binomial) # no covariates
sum <- data.frame(summary(fit_basic)$coefficients)
basic_p_hc <- sum["bio", ]

fit_p_hc <- glm(Status2 ~ bio + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5, data = merged, family = binomial) # basic covariates
sum <- data.frame(summary(fit_p_hc)$coefficients)
all_p_hc <- sum["bio", ]

# prodromals excluded
merged2 <- merged[merged$Status1 != 3,] 

fit_basic <-  glm(Status2 ~ bio, data = merged2, family = binomial)
sum <- data.frame(summary(fit_basic)$coefficients)
basic_no_p <- sum["bio", ]

fit_p_hc <- glm(Status2 ~ bio + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5, data = merged2, family = binomial)
sum <- data.frame(summary(fit_p_hc)$coefficients)
all_no_p <- sum["bio", ]

# Save
prodromal <- rbind(basic_p_hc, all_p_hc)
no_p <- rbind(basic_no_p, all_no_p)

write.table(prodromal, file = "prodromal_regression.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(no_p, file = "prodromal-removed_regression.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
