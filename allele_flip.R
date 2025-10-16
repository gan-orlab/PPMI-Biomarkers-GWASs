### Allele matching/flipping of summary stats
# Make sure all columns / column headers match
## ==== Load data ====
library(data.table)
library(dplyr)
meta1 <- fread("meta1.tab")
meta2 <- fread("meta2.tab")
    
# ==== Align to get all matching SNPs ====
merged <- merge(meta1, meta2, by.x=c('CHR','BP','minorAllele','majorAllele'), by.y=c('CHR','BP','minorAllele','majorAllele'), )
snps <- merged$SNP.x

# Pull only perfect matches
matches <- meta1 %>% filter(SNP %in% snps)
write.table(matches, "matches.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# Pull only mismatches, flip beta and freq
mismatches <- meta1 %>% filter(!SNP %in% snps)
mismatches$BETA <- -(mismatches$beta)
mismatches$FREQ <- 1-(mismatches$maf)
write.table(mismatches, "mismatches.txt", quote = F, sep = "\t", row.names = F, col.names = T)
