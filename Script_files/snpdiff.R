
library(tidyr)
library(devtools)

#setwd("")
AD <- read.table("out.AD.FORMAT", sep="\t", header = T)
#ADR <- read.table("out.ADR.FORMAT", sep="\t", header = T)
#ADF <- read.table("out.ADF.FORMAT", sep="\t", header = T)

cols <- c("chrom", "pos", "normal", "tumor")
colnames(AD) <- cols

#split the columns at a delimiter (:)  using the separate function from the tidyr package
AD_normal <- AD %>%
  
  separate(col = normal,
           into= c("nAD1", "nAD2", "nAD3", "nAD4"),
           sep = ",")
AD_all <- AD_normal %>%
  
  separate(col = tumor,
           into= c("tAD1", "tAD2", "tAD3", "tAD4"),
           sep = ",")
AD_cal <- AD_all[names(AD_all) %in% c("chrom", "pos","nAD1", "nAD2","tAD1", "tAD2")]

AD_cal$nAD1 <- as.numeric(AD_cal$nAD1)
AD_cal$nAD2 <- as.numeric(AD_cal$nAD2)
AD_cal$tAD1 <- as.numeric(AD_cal$tAD1)
AD_cal$tAD2 <- as.numeric(AD_cal$tAD2)

AD_cal$ncount <- AD_cal$nAD1 + AD_cal$nAD2
AD_cal$nvaf <- AD_cal$nAD2 / AD_cal$ncount
AD_cal$tcount <- AD_cal$tAD1 + AD_cal$tAD2
AD_cal$tvaf <- AD_cal$tAD2 / AD_cal$tcount

snpdiff <- AD_cal[names(AD_cal) %in% c("chrom", "pos", "ncount", "nvaf", "tcount", "tvaf")]
write.table(snpdiff, file="snpdiff.txt", sep= "\t", quote=F,col.names = T, row.names = F)

