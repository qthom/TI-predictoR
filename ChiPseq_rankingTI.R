#Wang1 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Wang 2015 (PMID 26100864)/Wang2015_H3K4me1_merged.bedgraph.gz')
#Wang3 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Wang 2015 (PMID 26100864)/Wang2015_H3K4me3_merged.bedgraph.gz')

#load the relevant histone marks
Ina1 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Inagaki Kakutani 2017 (PMID 28100676)/H3K4me1_treat_pileup.bedgraph.gz',format="bedgraph")
Ina3 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Inagaki Kakutani 2017 (PMID 28100676)/H3K4me3_treat_pileup.bedgraph.gz',format="bedgraph")
Ina2 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Inagaki Kakutani 2017 (PMID 28100676)/H3K4me2_treat_pileup.bedgraph.gz',format="bedgraph")
k36m2 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Luo Lam 2013 (PMID 22962860)/Luo2012_H3K36me2_treat_pileup.bdg.gz',format="bedgraph")    
k36m3 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Luo Lam 2013 (PMID 22962860)/Luo2012_H3K36me3_treat_pileup.bdg.gz',format="bedgraph")
H3k27me3 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Liu Weigel 2016 (PMID 27225844)/Liu2016_h3k27me3_treat_pileup.bdg.gz', format='bedgraph')
H3k9me2 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/ChIP-Seq/Inagaki Kakutani 2017 (PMID 28100676)/H3K9me2_treat_pileup.bedgraph.gz')


#process them
Ina1 <- dropSeqlevels(Ina1,c("Pt","Mt"),pruning.mode ="coarse")
Ina3 <- dropSeqlevels(Ina3,c("Pt","Mt"),pruning.mode ="coarse")
Ina2 <- dropSeqlevels(Ina2,c("Pt","Mt"),pruning.mode ="coarse")
k36m2 <- dropSeqlevels(k36m2,c("Pt","Mt"),pruning.mode ="coarse")
k36m3 <- dropSeqlevels(k36m3,c("Pt","Mt"),pruning.mode ="coarse")
h3k27me3 <- dropSeqlevels(H3k27me3,c("Pt","Mt"),pruning.mode ="coarse")
h3k9me2 <- dropSeqlevels(H3k9me2,c("Pt","Mt"),pruning.mode ="coarse")

Ina.a1 <- resize(Ina1, width = 1, fix = "center")
Ina.a3 <- resize(Ina3, width = 1, fix = "center")
Ina.a2 <- resize(Ina2, width = 1, fix = "center")
k36.a2 <- resize(k36m2, width = 1, fix = "center")
k36.a3 <- resize(k36m3, width = 1, fix = "center")
h3k27me3.1 <- resize(h3k27me3, width = 1, fix = "center")
h3k9me2.1 <- resize(h3k9me2, width = 1, fix = "center")

#down select genes up to 1000bp
coding_genes_n0big <- coding_genes[width(coding_genes) >= 800, ]
#shrinks genes to gene bodies
shrinked_genes <- extend(coding_genes_n0big,upstream=-150,downstream=-150)
shrinked_genes <- dropSeqlevels(shrinked_genes,c("Mt","Pt"),pruning.mode ="coarse")
#define TSSs for the calculation
TSS.2 <- resize(coding_genes_n0big, width = 1 , fix = "start")
extended_TSS.2 <- extend(TSS.2,upstream=150,downstream=150)
extended_TSS.2 <- dropSeqlevels(extended_TSS.2,c("Mt","Pt"),pruning.mode ="coarse")
#create list of histone marks for calculation
data <- list("h3k4me1" = Ina.a1,"h3k4me3" = Ina.a3, "h3k4me2" = Ina.a2, "h3k27me3" = h3k27me3.1,"h3k9me2" = h3k9me2.1,"h3k36me2" = k36.a2,"h3k36me3" = k36.a3)

#calculate the TPMs for each regions of genes
K4ME.tss <- getOverlappingScores_v2(extended_TSS.2,data)
K4ME.genebody <- getOverlappingScores_v2(shrinked_genes,data)

#combine the genebody values to the TSS genomic range

K4ME.tss$h3k4me1.gb <- K4ME.genebody$h3k4me1
K4ME.tss$h3k4me3.gb <- K4ME.genebody$h3k4me3
K4ME.tss$h3k4me2.gb <- K4ME.genebody$h3k4me2
K4ME.tss$h3k27me3.gb <- K4ME.genebody$h3k27me3
K4ME.tss$h3k9me2.gb <- K4ME.genebody$h3k9me2

#calculates FPKM from TPM

K4ME.tss$h3k4me1.gb.norm <- K4ME.tss$h3k4me1.gb/width(K4ME.tss)
K4ME.tss$h3k4me3.gb.norm <- K4ME.tss$h3k4me3.gb/width(K4ME.tss)
K4ME.tss$h3k4me2.gb.norm <- K4ME.tss$h3k4me2.gb/width(K4ME.tss)
K4ME.tss$h3k27me3.gb.norm <- K4ME.tss$h3k27me3.gb/width(K4ME.tss)
K4ME.tss$h3k9me2.gb.norm <- K4ME.tss$h3k9me2.gb/width(K4ME.tss)

#gr <- gr_wt
#source(/home/bcm215/Desktop/Quentin/scripts/TIF-Seq_2019/TIF-forh3k4me1.R)
#K4ME.tss$overlap_candidatesTIF <- countOverlaps(K4ME.tss,candidates,type=c("any"),ignore.strand=FALSE)

dat <- as.data.frame(K4ME.tss)

#rank the genes by H3k4me marks ups and down (high h3k4me1 and low H3k4me3 and other repressive marks)
dat$r_h3k4me1 <- NA
dat$r_h3k4me1[order(dat$h3k4me1)] <- 1:nrow(dat)/nrow(dat)
dat$r_h3k4me3 <- NA
dat$r_h3k4me3[order(-dat$h3k4me3)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k4me2 <- NA
dat$r_h3k4me2[order(-dat$h3k4me2)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k9me2 <- NA
dat$r_h3k9me2[order(-dat$h3k9me2)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k27me3 <- NA
dat$r_h3k27me3[order(-dat$h3k27me3)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k36me2 <- NA
dat$r_h3k36me2[order(dat$h3k36me2)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k36me3 <- NA
dat$r_h3k36me3[order(dat$h3k36me3)] <- (1:nrow(dat))/nrow(dat)

dat$r_h3k4me1.gb.norm <- NA
dat$r_h3k4me1.gb.norm[order(dat$h3k4me1.gb.norm)] <- 1:nrow(dat)/nrow(dat)
dat$r_h3k4me3.gb.norm <- NA
dat$r_h3k4me3.gb.norm[order(-dat$h3k4me3.gb.norm)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k4me2.gb.norm <- NA
dat$r_h3k4me2.gb.norm[order(-dat$h3k4me2.gb.norm)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k9me2.gb.norm <- NA
dat$r_h3k9me2.gb.norm[order(-dat$h3k9me2.gb.norm)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k27me3.gb.norm <- NA
dat$r_h3k27me3.gb.norm[order(-dat$h3k27me3.gb.norm)] <- (1:nrow(dat))/nrow(dat)

#calculate h3k4me1/h3k4me3 ratio and rank it

dat$ra13 <- dat$h3k4me1/dat$h3k4me3
dat$ra12 <- dat$h3k4me1/dat$h3k4me2
dat$ra36_23 <- dat$h3k36me2/dat$h3k36me3

dat$r_ra13 <- NA
dat$r_ra13[order(dat$ra13)] <- (1:nrow(dat))/nrow(dat)
dat$r_ra12 <- NA
dat$r_ra12[order(dat$ra12)] <- (1:nrow(dat))/nrow(dat)
dat$r_ra36_23 <- NA
dat$r_ra36_23[order(dat$ra36_23)] <- (1:nrow(dat))/nrow(dat)

#Import pNET-Seq data for gene expression ranking
pNET <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Maxim/OUR DATA/Our plaNET-Seq/Merged/Merged runs and replicates/Bedgraph/WT_merged.bedgraph.gz', format ='bedgraph')
pnet <- dropSeqlevels(pNET,c("Pt","Mt"),pruning.mode ="coarse")

seqlevels(pnet) <- c("1", "2", "3", "4", "5")

data <- list("pNET" = pnet)
Norm <- getOverlappingScores_v2(shrinked_genes,data)
txdb2 <- dropSeqlevels(txdb,c("6","7"),pruning.mode ="coarse")
Norm$normpNET <- Norm$pNET/width(Norm)

dat$pNET <- Norm$normpNET
dat$r_pNET <- NA
dat$r_pNET[order(dat$pNET)] <- (1:nrow(dat))/nrow(dat)

#subset genes with high ration and high h3k4me1

dat.2 <- dat[dat$ra13 > 1,]
dat.2 <- dat.2[dat.2$ra12 > 1,]
dat.2 <- dat.2[dat.2$h3k4me1 > 30,]


#calculate a score to find best candidates

s1 = 20 #ratio me1/me3
s1.2 = 10 #ratio me1/me2
s2 = 12 #h3k4me1
s3 = 10 #h3k4me3
s4 = 6 #h3k27me3
s5 = 6 #h3k9me2
s6 = 6 #h3k4me2
s7 = 10 #h3k36me3
s8 = 5 #h3k36me2 
s9 = 12 #ratio h3k36me2/me3

dat.2$score2 <- (s1*(dat.2$r_ra13) + s1.2*(dat.2$r_ra12) + s2*(dat.2$r_h3k4me1) + s3*(dat.2$r_h3k4me3) + s4*(dat.2$r_h3k27me3) + s5*(dat.2$r_h3k9me2) + s6*(dat.2$r_h3k4me2) + s7*(dat.2$r_h3k36me3) + s8*(dat.2$r_h3k36me2) + s9*(dat.2$r_ra36_23))


tair_ann <- import.gff3("/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/sequencing_files/TAIR10/Arabidopsis_thaliana.TAIR10.26.gff3")
mapping <- data.frame("geneID" = sub("gene:", "", tair_ann$ID), "name"=tair_ann$external_name)
mapping$gene_id <- mapping$geneID
dat.4 <- left_join(x=dat.2, y=mapping, by="gene_id", na_matches = "never")
dat.4 <- dat.4[complete.cases(dat.4$gene_id), ]

#add phenotypes

#pheno <- read.csv("/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/phenotypes.csv",header=TRUE, sep=",")
#dat.5 <- left_join(x=dat.4, y=pheno, by="gene_id", na_matches = "never")

#export table

write.table(dat.4, file="histone_forTI.csv", quote=F, sep="\t", row.names=F, col.names=T)














source("/home/bcm215/Desktop/Quentin/scripts/CAGEfightR.R")

df.3$fact <- match(df.3$gene_id, fa$geneID,nomatch=0) 
a <- df.3$fact != 0
df.3$fact[a] <- 1

write.table(df.3, file="H3K4ME_longupstream.csv", quote=F, sep="\t", row.names=F, col.names=T)

ndf.3 <- df.3[order(-df.3$ratio),]
nb <- nrow(ndf.3)
nb.up <- round(0.25*nb)
h3.14 <- df.3[1:nb.up,]
h3.24 <- df.3[nb.up:(nb.up*2),]
h3.34 <- df.3[(nb.up*2):(nb.up*3),]
h3.44 <- df.3[(nb.up*3):nb,]


coding_genes_nO$highratio14 <- match(coding_genes_nO$gene_id,h3.14$gene_id,nomatch=0)
genes_14 <- coding_genes_nO[coding_genes_nO$highratio14 > 0,]
coding_genes_nO$highratio24 <- match(coding_genes_nO$gene_id,h3.24$gene_id,nomatch=0)
genes_24 <- coding_genes_nO[coding_genes_nO$highratio24 > 0,]
coding_genes_nO$highratio34 <- match(coding_genes_nO$gene_id,h3.34$gene_id,nomatch=0)
genes_34 <- coding_genes_nO[coding_genes_nO$highratio34 > 0,]
coding_genes_nO$highratio44 <- match(coding_genes_nO$gene_id,h3.44$gene_id,nomatch=0)
genes_44 <- coding_genes_nO[coding_genes_nO$highratio44 > 0,]

TSS_14 <- resize(genes_14, width = 1 , fix = "start")
TSS_24 <- resize(genes_24, width = 1 , fix = "start")
TSS_34 <- resize(genes_34, width = 1 , fix = "start")
TSS_44 <- resize(genes_44, width = 1 , fix = "start")

makebed(TSS_14,"TSS_HR_14.bed")
makebed(TSS_24,"TSS_HR_24.bed")
makebed(TSS_34,"TSS_HR_34.bed")
makebed(TSS_44,"TSS_HR_44.bed")