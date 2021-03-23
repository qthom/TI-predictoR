# Load the required libraries

library(CAGEfightR)						# version 0.99.0
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BiocParallel)
register(MulticoreParam(4), default=T)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)

source("/home/bcm215/Desktop/Quentin/scripts/Nielsen_et_al_2018-master/assignTxType_custom.R")

# Load the TSS-Seq BigWig files (see the 01-Alignment_of_5Cap-Seq_data.sh pipeline):

setwd("/home/bcm215/Desktop/sequencing_files/CAP-Seq/bedgraph_capseq/expanded_ssrp1/")
bw_plus_filenames <- list.files(".", pattern="fw_cov3_expanded.bw$")
bw_minus_filenames <- list.files(".", pattern="rev_cov3_expanded.bw$")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- sub('fw_cov3_expanded.bw', '', bw_plus_filenames)
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# Make the design matrix:
design_matrix <- data.frame("Name"=sample_names, "BigWigPlus"=bw_plus_filenames, "BigWigMinus"=bw_minus_filenames,
row.names=sample_names, genotype=(c("col0","fact","fact","col0","fact","fact")),fact = (c("col0","ssrp1","spt16","col0","ssrp1","spt16")))


# Quantify all tag clusters (TCs):
ctss <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, design=design_matrix, genome=seqinfo(Athaliana))

# Call candidate TSS:
tss <- quickTSSs(ctss)


# Annotate TSS and enhancers by genomic features (observe that a custom assignTxType() function is used):
rowRanges(tss)$txType <- suppressWarnings(assignTxType_custom(rowRanges(tss), txdb=txdb, asFactor=TRUE))

# Combine candidate TSS and enhancers into a single RangedSummarizedExperiment object:
rowRanges(tss)$clusterType <- "TSS"

rse <- combineClusters(tss, tss, removeIfOverlapping="object1")

# Remove low expressed TCs:
rse <- subsetBySupport(rse, inputAssay = "counts", outputColumn = "support", unexpressed = 0, minSamples = 1) # n = 96232

# Annotate TCs by gene IDs:
rse <- suppressWarnings(assignGeneID(rse, geneModels=txdb))

# Annotate TCs by gene names:
tair_ann <- import.gff3("/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/sequencing_files/TAIR10/Arabidopsis_thaliana.TAIR10.26.gff3")
mapping <- data.frame("geneID" = sub("gene:", "", tair_ann$ID), "name"=tair_ann$external_name)
tmp <- left_join(x=tibble(geneID=rowRanges(rse)$geneID), y=mapping, by="geneID", na_matches = "never")
mcols(rse) <- DataFrame(mcols(rse), tmp[,-1])
rm(tmp)

# For DE calling consider only strong peaks (TPM >= 1 in at least 2 samples):
rse2 <- rse[rowSums(cpm(assay(rse)) >= 1) >= 2]	# n = 25964

# Use DESeq2:
dds <- DESeqDataSet(rse2, design = ~genotype)
dds <- DESeq(dds)

dds.2 <- DESeqDataSet(rse2, design = ~fact)
dds.2 <- DESeq(dds.2)

# Extract DE results:

fact <- results(dds, contrast = c("genotype", "fact", "col0"))
summary(fact)
ssrp1 <- results(dds.2, contrast = c("fact", "ssrp1", "col0"))
summary(ssrp1)
spt16 <- results(dds.2, contrast = c("fact", "spt16", "col0"))
summary(spt16)

# Combine DE results:

df2 <- as.data.frame(fact)[,c(2,6)]
colnames(df2) <- c("log2FC_fact", "padj_fact")
df3 <- as.data.frame(ssrp1)[,c(2,6)]
colnames(df3) <- c("log2FC_ssrp1", "padj_ssrp1")
df4 <- as.data.frame(spt16)[,c(2,6)]
colnames(df4) <- c("log2FC_spt16", "padj_spt16")
m <- cbind(df2, df3, df4)
m[is.na(m)] <- 1			# replace padj=NA with padj=1

# Expand DE results to the original TC number:
m$key <- rownames(m)
orig <- data.frame("key"=names(rowRanges(rse)))
n <- left_join(orig, m, by=c("key"), all.x=TRUE)

# Add DE results to the RSE object:
mcols(rse) <- cbind(mcols(rse), n[-1])
mc <- mcols(rse)


mcols(rse)$fact_de <- "nonDE"
fact_up <- !is.na(mc$padj_fact) & mc$padj_fact<=0.1 & mc$log2FC_fact > 0
mcols(rse)$fact_de[fact_up] <- "up"
fact_down <- !is.na(mc$padj_fact) & mc$padj_fact<=0.01 & mc$log2FC_fact < 0
mcols(rse)$fact_de[fact_down] <- "down"

mcols(rse)$ssrp1_de <- "nonDE"
ssrp1_up <- !is.na(mc$padj_ssrp1) & mc$padj_ssrp1<=0.1 & mc$log2FC_ssrp1 > 0
mcols(rse)$ssrp1_de[ssrp1_up] <- "up"
ssrp1_down <- !is.na(mc$padj_ssrp1) & mc$padj_ssrp1<=0.1 & mc$log2FC_ssrp1 < 0
mcols(rse)$ssrp1_de[ssrp1_down] <- "down"

mcols(rse)$spt16_de <- "nonDE"
spt16_up <- !is.na(mc$padj_spt16) & mc$padj_spt16<=0.1 & mc$log2FC_spt16 > 0
mcols(rse)$spt16_de[spt16_up] <- "up"
spt16_down <- !is.na(mc$padj_spt16) & mc$padj_spt16<=0.1 & mc$log2FC_spt16 < 0
mcols(rse)$spt16_de[spt16_down] <- "down"

mc <- mcols(rse)

fact_only <- rse[mc$spt16_de == "up" | mc$ssrp1_de == "up"]
ssrp1_only <- rse[mc$ssrp1_de == "up"]
spt16_only <- rse[mc$spt16_de == "up"]

fact_upregulation <- fact_only[mcols(fact_only)$txType == "promoter" | mcols(fact_only)$txType == "fiveUTR"]
fa <- mcols(fact_upregulation)
ssrp1_upregulation <- ssrp1_only[mcols(ssrp1_only)$txType == "promoter" | mcols(ssrp1_only)$txType != "fiveUTR"]
ss <- mcols(ssrp1_upregulation)
spt16_upregulation <- spt16_only[mcols(spt16_only)$txType == "promoter"]
sp <- mcols(spt16_upregulation)
