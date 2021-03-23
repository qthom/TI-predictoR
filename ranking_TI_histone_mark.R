#load the script TIF_pipeline to "Process_TIF.R"


Wang1 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Wang 2015 (PMID 26100864)/Wang2015_H3K4me1_merged.bedgraph.gz')
Wang3 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Wang 2015 (PMID 26100864)/Wang2015_H3K4me3_merged.bedgraph.gz')

H3k27me3 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Zhou 2017 (PMID 28403905)/Zhou2017_H3K27me3_merged.bedgraph.gz')
H3k9me2 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Bewick 2016 (PMID 27457936) - H3K56ac ChIP-Seq/H3K9me2_macs.bedgraph.gz')

wang1 <- dropSeqlevels(Wang1,c("Pt","Mt"),pruning.mode ="coarse")
wang3 <- dropSeqlevels(Wang3,c("Pt","Mt"),pruning.mode ="coarse")
h3k27me3 <- dropSeqlevels(H3k27me3,c("Pt","Mt"),pruning.mode ="coarse")
h3k9me2 <- dropSeqlevels(H3k9me2,c("Pt","Mt"),pruning.mode ="coarse")

wang.a1 <- resize(wang1, width = 1, fix = "center")
wang.a3 <- resize(wang3, width = 1, fix = "center")
h3k27me3.1 <- resize(h3k27me3, width = 1, fix = "center")
h3k9me2.1 <- resize(h3k9me2, width = 1, fix = "center")

coding_genes_n0big <- coding_genes_nO[width(coding_genes_nO) > 1000, ]
shrinked_genes <- extend(coding_genes_n0big,upstream=-350,downstream=-150)
shrinked_genes <- dropSeqlevels(shrinked_genes,c("Mt","Pt"),pruning.mode ="coarse")
TSS.2 <- resize(coding_genes_n0big, width = 1 , fix = "start")
extended_TSS.2 <- extend(TSS.2,upstream=100,downstream=100)
extended_TSS.2 <- dropSeqlevels(extended_TSS.2,c("Mt","Pt"),pruning.mode ="coarse")
data <- list("h3k4me1" = wang.a1,"h3k4me3" = wang.a3,"h3k27me3" = h3k27me3.1,"h3k9me2" = h3k9me2.1)
seqinfo(wang.a1) <- seqinfo(extended_TSS.2)
seqinfo(wang.a3) <- seqinfo(extended_TSS.2)
seqinfo(h3k27me3.1) <- seqinfo(extended_TSS.2)
seqinfo(h3k9me2.1) <- seqinfo(extended_TSS.2)
K4ME.tss <- getOverlappingScores_v2(extended_TSS.2,data)

dat <- as.data.frame(K4ME.tss)

dat$r_h3k4me1 <- NA
dat$r_h3k4me1[order(dat$h3k4me1)] <- 1:nrow(dat)/nrow(dat)
dat$r_h3k4me3 <- NA
dat$r_h3k4me3[order(dat$h3k4me3)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k9me2 <- NA
dat$r_h3k9me2[order(dat$h3k9me2)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k27me3 <- NA
dat$r_h3k27me3[order(dat$h3k27me3)] <- (1:nrow(dat))/nrow(dat)
dat$scores <- NA
dat$scores <- (dat$r_h3k4me1/rowSums(dat[14:15:16]))
write.table(dat, file="histone_forTI.csv", quote=F, sep="\t", row.names=F, col.names=T)

head(dat)
dat$ratio <- dat$h3k4me1/dat$h3k4me3
dat.2 <- dat[order(-dat$h3k4me1,-dat$ratio,-dat$r_h3k9me2,-dat$r_h3k27me3),] 

write.table(dat.2, file="histone_forTI_rank2.csv", quote=F, sep="\t", row.names=F, col.names=T)
