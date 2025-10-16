### Creat fullSTATS and plots from metal output
## ==== Load data ====
require(data.table)
require(qqman)
require(ggplot2)
require(dplyr)
assoc <- fread("meta.tab") # this is your output from metal, with CHR and BP split into separate columns
colnames(assoc) <- c("CHR","BP","SNP","A1","A2","Freq1","FreqSE","MinFreq","MaxFreq","Effect","StdErr","P","Direction","HetISq","HetChiSq","HetDf","HetPval")

## ==== Prepare summary statistics ====
# Format
assoc$minorAllele <- ifelse(assoc$Freq1 <= 0.5, as.character(assoc$A1), as.character(assoc$A2))
assoc$majorAllele <- ifelse(assoc$Freq1 <= 0.5, as.character(assoc$A2), as.character(assoc$A1))
assoc$beta <- ifelse(assoc$Freq1 <= 0.5, assoc$Effect, -(assoc$Effect))
assoc$maf <- ifelse(assoc$Freq1 <= 0.5, assoc$Freq1, 1 - assoc$Freq1)
assoc$se <- assoc$StdErr
assoc$OR <- exp(assoc$beta)
assoc$zscore = assoc$beta/assoc$se

# Save
fullSTATS = assoc[,c("CHR","BP","SNP","minorAllele","majorAllele","maf","beta","se","OR","zscore","P")]
write.table(fullSTATS, "fullSTATS.meta.tab", quote = F, sep = "\t", row.names = F)

## ==== Prepare Manhattan and QQ plots ====
gwasResults = fullSTATS[,c("SNP", "CHR", "BP", "P", "zscore")]
test1 <- subset(gwasResults, P < 1e-4)
test2 <- subset(gwasResults, P > 1e-4)
length <- (nrow(test2))/2
test3 <- sample_n(test2, length, replace = FALSE)
subGWAS <- rbind(test1, test3)

# Plot and save QQ plot 
tiff(file="QQ.meta.tiff", width = 8, height = 6, unit = "in", res = 800, pointsize = 6)
qq(subGWAS$P)
dev.off()

# Plot and save Manhattan plot
tiff(file="ManH.meta.tiff", width = 8, height = 6, unit = "in", res = 800, pointsize = 6)
manhattan(subGWAS)
dev.off()
