setwd("H:/plink_win64/")

library(permute)
library(lattice)
library(vegan)
library(pheatmap)

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
p17 <- read.table("p17.txt")
p18 <- read.table("p18.txt")
p19 <- read.table("p19.txt")
p20 <- read.table("p20.txt")
p21 <- read.table("p21.txt")
p22 <- read.table("p22.txt")
p23 <- read.table("p23.txt")
p24 <- read.table("p24.txt")
p25 <- read.table("p25.txt")
p26 <- read.table("p26.txt")
p27 <- read.table("p27.txt")

marker <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27)
marker <- marker[!duplicated(marker[,1]),1]
write.table(marker, "EIRA_26abs_case-ctrl.txt", quote=F, col.names=F, row.names=F)

#univariate logistic regression for all significant markers#
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
p17 <- read.table("plink.P17.assoc.logistic", header=T)
p18 <- read.table("plink.P18.assoc.logistic", header=T)
p19 <- read.table("plink.P19.assoc.logistic", header=T)
p20 <- read.table("plink.P20.assoc.logistic", header=T)
p21 <- read.table("plink.P21.assoc.logistic", header=T)
p22 <- read.table("plink.P22.assoc.logistic", header=T)
p23 <- read.table("plink.P23.assoc.logistic", header=T)
p24 <- read.table("plink.P24.assoc.logistic", header=T)
p25 <- read.table("plink.P25.assoc.logistic", header=T)
p26 <- read.table("plink.P26.assoc.logistic", header=T)
p27 <- read.table("plink.P27.assoc.logistic", header=T)
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
p17 <- p17[which(p17$TEST=="ADD"),]
p18 <- p18[which(p18$TEST=="ADD"),]
p19 <- p19[which(p19$TEST=="ADD"),]
p20 <- p20[which(p20$TEST=="ADD"),]
p21 <- p21[which(p21$TEST=="ADD"),]
p22 <- p22[which(p22$TEST=="ADD"),]
p23 <- p23[which(p23$TEST=="ADD"),]
p24 <- p24[which(p24$TEST=="ADD"),]
p25 <- p25[which(p25$TEST=="ADD"),]
p26 <- p26[which(p26$TEST=="ADD"),]
p27 <- p27[which(p27$TEST=="ADD"),]
###################################################
eira_pvalue <- data.frame(p1[,2],p1[,9],p2[,9],p3[,9],p4[,9],p5[,9],p6[,9],p7[,9],p8[,9],p9[,9],p10[,9],p11[,9],p12[,9],p13[,9],p14[,9],p15[,9],p16[,9],p17[,9],p18[,9],p19[,9],p20[,9],p21[,9],p22[,9],p23[,9],p24[,9],p25[,9],p26[,9],p27[,9])
eira_or <- data.frame(p1[,2],p1[,7],p2[,7],p3[,7],p4[,7],p5[,7],p6[,7],p7[,7],p8[,7],p9[,7],p10[,7],p11[,7],p12[,7],p13[,7],p14[,7],p15[,7],p16[,7],p17[,7],p18[,7],p19[,7],p20[,7],p21[,7],p22[,7],p23[,7],p24[,7],p25[,7],p26[,7],p27[,7])
eira_bi <- ifelse(eira_pvalue[,2:28]<0.05/23, 1, 0)
eira_sig <- data.frame(eira_pvalue[,1], eira_bi[,1:27])

write.table(eira_pvalue, "EIRA_26abs_case-ctrl_pvalue.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(eira_or, "EIRA_26abs_case-ctrl_OR.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(eira_sig, "EIRA_26abs_case-ctrl_sig.txt", quote=F, row.names=F, col.names=F, sep="\t")

#Heatmap based on Jaccard distance#
ab27 <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm15", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit", "Fib 72 cit", "Fib 74 cit", "Carb-CEP-1", "CPP-3", "CPP-5", "Pept-Z1", "Pept-Z2", "Pept-1", "Pept-5", "Bla 26", "Ab-neg")
ab26 <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit", "Fib 72 cit", "Fib 74 cit", "Carb-CEP-1", "CPP-3", "CPP-5", "Pept-Z1", "Pept-Z2", "Pept-1", "Pept-5", "Bla 26", "Ab-neg")
sig <- eira_sig[,2:28]
row.names(sig) <- eira_sig[,1]
names(sig) <- ab27
#sig <- sig[,-10]
tsig <- t(sig)
sig_dist <- vegdist(tsig, method="jaccard")
pheatmap(sig, scale="none", clustering_distance_col=sig_dist, cluster_rows=F, cluster_cols=T, fontsize=15, fontsize_rows=3)

#Recode sig based on direction of ORs#
eira_sig_new <- eira_sig
eira_sig_new[10,9:10] <- c(-1,-1)
eira_sig_new[10,12:15] <- c(-1,-1,-1,-1)
eira_sig_new[c(7,10,23),14] <- c(-1,-1,-1)
eira_sig_new[20,18] <- c(-1)
eira_sig_new[10,22] <- c(-1)
eira_sig_new[c(1,2,4,6,7,9,11:14,16:23,25),28] <- rep(-1,19)

narac_sig_new <- narac_sig
narac_sig_new[1:3,28] <- c(-1,-1,-1)
narac_sig_new[5,c(9,10,14,15)] <- c(-1,-1,-1,-1)
narac_sig_new[6,28] <- -1
narac_sig_new[9,18] <- -1
narac_sig_new[9:14,28] <- c(-1,-1,-1,-1,-1,-1)

write.table(eira_sig_new, "EIRA_26abs_case-ctrl_sig+OR.txt", quote=F, row.names=F, col.names=F, sep="\t")
##############################################
sig <- eira_sig_new[,2:28]
row.names(sig) <- eira_sig_new[,1]
names(sig) <- ab27
#sig <- sig[,-27]
#sig <- sig[-8,]
tsig <- t(sig)
sig_dist <- vegdist(tsig, method="jaccard")
pheatmap(sig, scale="none", clustering_distance_col=sig_dist, cluster_rows=F, cluster_cols=T, fontsize=15, fontsize_rows=3)
