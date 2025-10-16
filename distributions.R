### Creat distribution graphs and run shapiro test
## ==== Load data ====
library(data.table)

line <- read.table("marker.txt")
dat <- fread("biomarker_data.txt")
dat <- na.omit(dat)
colnames(dat) <- c("FID","IID","bio")

## ==== Run Shapiro ====
res <- shapiro.test(dat$bio)
p <- res$p.value
p_table <- data.frame(line, p)
write.table(p_table, file = "shapiro.txt", sep = "\t", row.names = F, col.names = F, quote = F)

## ==== Make histogram ====
tiff(file="bio_hist.tiff", width = 8, height = 6, unit = "in", res = 800, pointsize = 6)
hist(dat$bio, col = "mediumpurple3")
dev.off()
