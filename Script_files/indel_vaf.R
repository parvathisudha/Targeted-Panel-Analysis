library(tidyr)
library(devtools)

#setwd("")
TAR <- read.table("out.TAR.FORMAT", sep="\t", header = T)
TIR <- read.table("out.TIR.FORMAT", sep="\t", header = T)
#split the columns at a delimiter (:)  using the separate function from the tidyr package
TAR_tumor <- TAR %>%
  
  separate(col = TUMOR,
           into= c("tTAR1", "tTAR2"),
           sep = ",")

TIR_tumor <- TIR %>%
  
  separate(col = TUMOR,
           into= c("tTIR1", "tTIR2"),
           sep = ",")

indel <- merge(TAR_tumor,TIR_tumor, by= c("CHROM", "POS"))

indel$tTAR1 <- as.numeric(indel$tTAR1) 
indel$tTAR2 <- as.numeric(indel$tTAR2) 

indel$tTIR1 <- as.numeric(indel$tTIR1) 
indel$tTIR2 <- as.numeric(indel$tTIR2) 

vaf <- indel$tTIR1/(indel$tTIR1+indel$tTAR1)
vaf
indel$vaf <- vaf
write.table(indel, file="indel_vaf.txt", sep="\t", quote=F, col.names = T, row.names = F)
