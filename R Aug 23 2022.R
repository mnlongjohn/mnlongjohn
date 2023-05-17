#   Differential expression analysis using z-scores with limma

#Load packages
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE31376", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8227", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("4X4X3544445X31XX14X13100000111111222215344444XX5X3",
               "4X534X133XX22222333300000000")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

design <- model.matrix(~0+factor(c(4,4,3,5,4,4,4,4,5,3,1,1,4,1,3,1,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,1,5,3,4,4,4,4,4,5,3,4,5,3,4,1,3,3,2,2,2,2,2,3,3,3,3,0,0,0,0,0,0,0,0)))
head(design)

colnames(design) <- c("HC","HDP","MRD","BCR","ETV6","TCF3")
design

write.csv(exprs(gset), file = "forzscores4.csv")
Zscores313 <- read.csv("forzscores4.csv", header = T)
Zscores313

colnames(Zscores313)
row.names(Zscores313)

#Calculate Z.scores for each gene across the samples (per rows) according to info from https://www.biostars.org/p/136952/
#Calculate Z scores on log 2 transformed data, which is in line with Wang et al, 2014 PMID: 15231529
#Scale() function in T calculates z-score per column, hence to get z-scores from rows (genes), transpose first (and back afterwards) 
#How to transpose according to info from https://www.biostars.org/p/201179/
#Intro to Z-scores https://www.r-bloggers.com/how-to-compute-the-z-score-with-r/
t.Zscores313 <- t(Zscores313)
write.csv(t.Zscores313, file = "t.Zscores313.csv")

#Open on excel, delete x and additional rows, save as csv and read in
t.Zscores313 <- read.csv("t.Zscores313.2.csv", header = F)
t.Zscores313

#Calculate Z scores and write into csv
Z.scores.313 <- scale(t.Zscores313)
write.csv(Z.scores.313, file = "forzscores3.3.csv")

#Add columns and rows and read in csv
z.scores.313 <- read.csv("forzscores3.4.csv")
z.scores.313
Z.scores.313 <- t(z.scores.313)
write.csv(Z.scores.313, file = "zscoresforheatmaps31376.csv")

#Using z-scores to identify differentially expressed miRNAs according to http://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html

#Z scores for identifying DGE https://transmart-app.readthedocs.io/en/latest/advanced_workflow.html#z-score-calculation
#I don't understand this script though

#Presenting DGE data using heatmaps, plotting expression across multiple genes from different groups, etc according to https://hbctraining.github.io/GCC-BOSC-2018/lessons/data_visualization.html
# volcano plots https://www.rdocumentation.org/packages/a4Base/versions/1.20.0/topics/volcanoPlot
#volcano plots https://www.r-bloggers.com/using-volcano-plots-in-r-to-visualize-microarray-and-rna-seq-results/
#volcano plot in limma http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/volcanoplot.html
#volcano plot https://www.biostars.org/p/309962/
#volcano plot https://www.biostars.org/p/214100/
#Advanced flow for heatmaps, volcano plots, etc https://transmart-app.readthedocs.io/en/latest/advanced_workflow.html#z-score-calculation