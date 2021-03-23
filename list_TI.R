#this script makes the table for identifying the best TI targets


source("/home/bcm215/Desktop/Quentin/scripts/TIF_pipeline.R")
source("/home/bcm215/Desktop/Quentin/scripts/CAGEfightR.R")


gr <- gr_hen2
source("/home/bcm215/Desktop/Quentin/scripts/seqclusters_forTI.R")
gr <- TI
source("/home/bcm215/Desktop/Quentin/scripts/annotate.R")
TI_hen2 <- gr

gr <- gr_wt
source("/home/bcm215/Desktop/Quentin/scripts/seqclusters_forTI.R")
gr <- TI
source("/home/bcm215/Desktop/Quentin/scripts/annotate.R")
TI_wt <- gr

df <- data.frame("geneID" = genes$gene_id)

fact_genes <- fa$geneID
hen2prox_genes <- hen$geneID
ssrp1_genes <- ss$geneID
spt16_genes <- sp$geneID
hen2TI_genes <- TI_hen2$geneID
wtTI_genes <- TI_wt$geneID

df$CAP_fact <- match(df$geneID,fact_genes,nomatch=0)
df$CAP_hen2 <- match(df$geneID,hen2prox_genes,nomatch=0)
df$CAP_ssrp1 <- match(df$geneID,ssrp1_genes,nomatch=0)
df$CAP_spt16 <- match(df$geneID,spt16_genes,nomatch=0)
#df$TIF_wt <- match(df$geneID,wtTI_genes,nomatch=0)
#df$TIF_hen2 <- match(df$geneID,hen2TI_genes,nomatch=0)

cfa <- df$CAP_fact != 0
df$CAP_fact[cfa] <- 1
che <- df$CAP_hen2 != 0
df$CAP_hen2[che] <- 1
css <- df$CAP_ssrp1 != 0
df$CAP_ssrp1[css] <- 1
csp <- df$CAP_spt16 != 0
df$CAP_spt16[csp] <- 1
#twt <- df$TIF_wt != 0
#df$TIF_wt[twt] <- 1
#the <- df$TIF_hen2 != 0
#df$TIF_hen2[the] <- 1

df$sum <- rowSums(df[2:5])
df2<- df[order(df$sum),]
df3 <- df[sum > 4]


#write.table(df3, file="best_TI_candidates", quote=F, sep="\t", row.names=T, col.names=T
