setwd("H:/plink-1.07-dos/")

p1 <- read.table("p1.txt")
p2 <- read.table("p2.txt")
p3 <- read.table("p3.txt")
p4 <- read.table("p4.txt")
p5 <- read.table("p5.txt")
p6 <- read.table("p6.txt")
p7 <- read.table("p7.txt")
p8 <- read.table("p8.txt")
p9 <- read.table("p9.txt")
p11 <- read.table("p11.txt")
p12 <- read.table("p12.txt")
p13 <- read.table("p13.txt")
p14 <- read.table("p14.txt")
p16 <- read.table("p16.txt")

marker <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p11, p12, p13, p14, p16)
write.table(marker, "Narac_all-marker_case-only.txt", quote=F, col.names=F, row.names=F)

#univariate logistic regression for summerized 16 markers
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

narac <- data.frame(p1[,2],p1[,9],p2[,9],p3[,9],p4[,9],p5[,9],p6[,9],p7[,9],p8[,9],p9[,9],p10[,9],p11[,9],p12[,9],p13[,9],p14[,9],p15[,9],p16[,9])
write.table(narac, "Narac_all-marker_case-only.txt", quote=F, row.names=F, col.names=F, sep="\t")

narac_bi <- ifelse(narac[,2:17]<0.05/16, 1, 0)
narac_bi <- data.frame(narac[,1], narac_bi[,1:16])
write.table(narac_bi, "Narac_all-marker_bi_case-only.txt", quote=F, row.names=F, col.names=F, sep="\t")

#Heatmap to cluster abs based on significance
library(pheatmap)
test <- narac_bi[,2:17]
ab <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm15", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit")
row.names(test) <- narac_bi[,1]
names(test) <- ab
tmp <- narac[,2:17]
tmp <- -log10(tmp)
row.names(tmp) <- narac[,1]
names(tmp) <- ab
#Heatmap based on p-value
pheatmap(tmp, scale="column", clustering_distance_col="correlation", cluster_rows=T, cluster_cols=T, fontsize=15, fontsize_rows=5)
#Heatmap based on significance
pos_ctrl <- c(rep(1,16))
test_new <- rbind(test, pos_ctrl)
row.names(test_new)[17] <- "positive_control"
pheatmap(test_new, scale="column", clustering_distance_col="correlation", cluster_rows=T, cluster_cols=T, fontsize=15, fontsize_rows=5)







