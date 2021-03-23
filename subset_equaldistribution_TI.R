df <- data.frame("geneID" = genes$gene_id)
fact_genes <- fact_upregulated$geneID
wt_genes <- notfact_upregulated$geneID
df$wt_genes <- match(df$geneID, wt_genes,nomatch=0)
df$fact_genes <- match(df$geneID, fact_genes,nomatch=0)


MyData <- read.csv(file="/home/bcm215/Desktop/Quentin/TIF-Seq/analysis/gene_expression_levels.csv", header = F, sep=",")
MyData$geneID <- MyData$V1
MyData$V1 <- NULL

df2 <- merge(df, MyData, by='geneID', all=T)

df2 <- df2[complete.cases(df2$V2), ]

fact<- data.frame("geneID" = df2$geneID[df2$fact_genes > 0 ], "exp" = df2$V2[df2$fact_genes > 0])
wt <- data.frame("geneID" = df2$geneID[df2$wt_genes > 0 ], "exp" = df2$V2[df2$wt_genes > 0])
fact <- fact[complete.cases(fact$exp), ]
wt <- wt[complete.cases(wt$exp), ]

hen2_none <- merge(fact, wt, by='exp', all=T)
hen2_none <- hen2_none[complete.cases(hen2_none$geneID.x), ]
hen2_none$dup.x <- duplicated(hen2_none$geneID.x)
hen2_none$dup.y <- duplicated(hen2_none$geneID.y)
genes_fact  <- hen2_none$geneID.x[hen2_none$dup.x == FALSE & hen2_none$dup.y == FALSE]
genes_wt  <- hen2_none$geneID.y[hen2_none$dup.x == FALSE & hen2_none$dup.y == FALSE]
control <- data.frame("genes_fact" = genes_fact, "genes_wt" = genes_wt)
control <- control[complete.cases(control$genes_wt), ]

TSS.fact <- genes[control$genes_fact]
TSS.fact <- resize(TSS.fact, width = 1 , fix = "start")
TSS.wt <- genes[control$genes_wt]
TSS.wt <- resize(TSS.wt, width = 1 , fix = "start")

TTS.hen2 <- TTS[control$genes_hen2]
TTS.none <- genes[control$genes_none]
TTS.none <- resize(TSS.none, width = 1 , fix = "end")


center.hen2 <- genes[control$genes_hen2]
center.hen2 <- resize(center.hen2, width = 1 , fix = "center")
center.none <- genes[control$genes_none]
center.none <- resize(center.none, width = 1 , fix = "center")

grl <- list("Hen2 genes" = TSS.hen2 , "Control genes" = TSS.none)
#grl <- list("Hen2 genes" = TTS.hen2 , "Control genes" = TTS.none)
#grl <- list("Hen2 genes" = center.hen2 , "Control genes" = center.none)


grl <- list("hen2 end"= hen2.tts)
hen2.tts <- resize(small_hen2, width = 1 , fix = "end")
hen2.tss <- resize(small_hen2, width = 1 , fix = "start")
