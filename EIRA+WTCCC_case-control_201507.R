setwd("H:/plink_win64/")

p1 <- read.table("p1.txt")
p2 <- read.table("p2.txt")
p3 <- read.table("p3.txt")
p4 <- read.table("p4.txt")
p5 <- read.table("p5.txt")
p6 <- read.table("p6.txt")
p7 <- read.table("p7.txt")
p8 <- read.table("p8.txt")
p9 <- read.table("p9.txt")
p10 <- read.table("p10.txt")
p11 <- read.table("p11.txt")
p12 <- read.table("p12.txt")
p13 <- read.table("p13.txt")
p14 <- read.table("p14.txt")
p15 <- read.table("p15.txt")
p16 <- read.table("p16.txt")

marker <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16)
write.table(marker, "Eira+wtccc_all-marker_case-ctrl.txt", quote=F, col.names=F, row.names=F)

#univariate logistic regression for summerized 56 markers
p1 <- read.table("plink.P1.assoc.logistic", header=T)
p2 <- read.table("plink.P2.assoc.logistic", header=T)
p3 <- read.table("plink.P3.assoc.logistic", header=T)
p4 <- read.table("plink.P4.assoc.logistic", header=T)
p5 <- read.table("plink.P5.assoc.logistic", header=T)
p6 <- read.table("plink.P6.assoc.logistic", header=T)
p7 <- read.table("plink.P7.assoc.logistic", header=T)
p8 <- read.table("plink.P8.assoc.logistic", header=T)
p9 <- read.table("plink.P9.assoc.logistic", header=T)
p10 <- read.table("plink.P10.assoc.logistic", header=T)
p11 <- read.table("plink.P11.assoc.logistic", header=T)
p12 <- read.table("plink.P12.assoc.logistic", header=T)
p13 <- read.table("plink.P13.assoc.logistic", header=T)
p14 <- read.table("plink.P14.assoc.logistic", header=T)
p15 <- read.table("plink.P15.assoc.logistic", header=T)
p16 <- read.table("plink.P16.assoc.logistic", header=T)
###################################################
p1 <- p1[which(p1$TEST=="ADD"),]
p2 <- p2[which(p2$TEST=="ADD"),]
p3 <- p3[which(p3$TEST=="ADD"),]
p4 <- p4[which(p4$TEST=="ADD"),]
p5 <- p5[which(p5$TEST=="ADD"),]
p6 <- p6[which(p6$TEST=="ADD"),]
p7 <- p7[which(p7$TEST=="ADD"),]
p8 <- p8[which(p8$TEST=="ADD"),]
p9 <- p9[which(p9$TEST=="ADD"),]
p10 <- p10[which(p10$TEST=="ADD"),]
p11 <- p11[which(p11$TEST=="ADD"),]
p12 <- p12[which(p12$TEST=="ADD"),]
p13 <- p13[which(p13$TEST=="ADD"),]
p14 <- p14[which(p14$TEST=="ADD"),]
p15 <- p15[which(p15$TEST=="ADD"),]
p16 <- p16[which(p16$TEST=="ADD"),]

###################################################
eira_wtccc <- data.frame(p1[,2],p1[,9],p2[,9],p3[,9],p4[,9],p5[,9],p6[,9],p7[,9],p8[,9],p9[,9],p10[,9],p11[,9],p12[,9],p13[,9],p14[,9],p15[,9],p16[,9])
dis_or <- data.frame(p1[,2],p1[,7],p2[,7],p3[,7],p4[,7],p5[,7],p6[,7],p7[,7],p8[,7],p9[,7],p10[,7],p11[,7],p12[,7],p13[,7],p14[,7],p15[,7],p16[,7])
rep_pvalue <- data.frame(p1[,2],p1[,9],p2[,9],p3[,9],p4[,9],p5[,9],p6[,9],p7[,9],p8[,9],p9[,9],p10[,9],p11[,9],p12[,9],p13[,9],p14[,9],p15[,9],p16[,9])
rep_or <- data.frame(p1[,2],p1[,7],p2[,7],p3[,7],p4[,7],p5[,7],p6[,7],p7[,7],p8[,7],p9[,7],p10[,7],p11[,7],p12[,7],p13[,7],p14[,7],p15[,7],p16[,7])
rep_pvalue <- rep_pvalue[,-11]
rep_or <- rep_or[,-11]

write.table(eira_wtccc, "Eira+wtccc_all-marker_case-ctrl.txt", quote=F, row.names=F, col.names=F, sep="\t")

eira_wtccc_bi <- ifelse(eira_wtccc[,2:17]<0.05/56, 1, 0)
eira_wtccc_bi <- data.frame(eira_wtccc[,1], eira_wtccc_bi[,1:16])
write.table(eira_wtccc_bi, "Eira+wtccc_all-marker_bi_case-ctrl.txt", quote=F, row.names=F, col.names=F, sep="\t")

#Heatmap to cluster abs based on significance
library(pheatmap)
test <- eira_wtccc_bi[,2:17]
ab <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm15", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit")
row.names(test) <- eira_wtccc_bi[,1]
names(test) <- ab
tmp <- eira_wtccc[,2:17]
tmp <- -log10(tmp)
row.names(tmp) <- eira_wtccc[,1]
names(tmp) <- ab
pheatmap(test, scale="column", clustering_distance_col="correlation", cluster_rows=T, cluster_cols=T, fontsize=12, fontsize_rows=3)
pheatmap(tmp, scale="none", clustering_distance_col="euclidean", cluster_rows=T, cluster_cols=T, fontsize=12, fontsize_rows=3)
#clustering_distance_col="euclidean"? or "correlation" if use p-values
#for binary data, scale="none" -> does not change if scale="column"
#for jaccard distance, clustering_distance_col=sig_dist
pheatmap(sig, scale="none", clustering_distance_col=sig_dist, cluster_rows=F, cluster_cols=T, fontsize=12, fontsize_rows=3)

#Heatmap based on jaccard distance of significance
setwd("H:/Project Results/")

library(openxlsx)
library(permute)
library(lattice)
library(vegan)

sig <- read.xlsx("eira+wtccc_1509_zm.xlsx", sheet="eira+wtccc_case-ctrl_summary", rowNames=T)
#For narac, removing ptm15: 
sig <- sig[,-10]
tmpr <- t(sig)
sig_dist <- vegdist(tmpr, method="jaccard")

sig <- sig[-10,]
sig_table <- sapply(sig, as.numeric)
row.names(sig_table) <- row.names(sig)
names(sig_table) <- names(sig)
pheatmap(sig_table, scale="none", clustering_distance_col="euclidean", cluster_rows=T, cluster_cols=T, fontsize=12, fontsize_rows=3)
