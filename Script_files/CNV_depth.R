Tumor <- read.table("Tumor_Sample.targetcoverage.cnn", sep="\t", header = T)
Normal <- read.table("Normal_Sample.targetcoverage.cnn", sep="\t", header = T)

t_cols <- c("chrom", "start","end", "id", "tumor_mean", "t_log2")
n_cols <- c("chrom", "start","end", "id", "normal_mean", "n_log2")

colnames(Tumor) <- t_cols
colnames(Normal) <- n_cols

Depth <- merge(Normal,Tumor, by= c("chrom", "start","end","id"), sort = F)

depth <- Depth[names(Depth) %in% c("chrom", "start","end", "id", "normal_mean", "tumor_mean")]
depth$tumor_depth <- NA
write.table(depth, file="depth.txt", sep="\t", quote=F, col.names = T, row.names = F)
