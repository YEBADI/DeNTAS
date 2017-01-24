##### RNA-seq analysis of FPKM data using limma & edgeR

rm (list =ls())  # clear the workkspace

#### installing the required packages if needed
#install.packages('limma')
#install.packages('edgeR')
#install.packages('dendextend')
#install.packages('gplots')
#install.packages('calibrate')

### loading the required packages
library(limma)
library(edgeR)
library(dendextend)
library(gplots)
library(calibrate)


# Reading in the variables set by python
myArgs <- commandArgs(trailingOnly = TRUE)
genes  <- read.table(file=myArgs[1], fill = TRUE)
myArgs <- myArgs[2:length(myArgs)]

x <- length(myArgs)
files  <- c(myArgs[1:(x/2)])
group  <- as.factor(c(myArgs[((x/2)+1):x]))

################ loading in the data ####################

### readDGE function creates a DGEList-object containing 1) sample info  & 2) expression matrix 
x <- readDGE(files, columns=c(1,3))
# class(x) # inspecting the data
# dim(x) # inspecting the data
# x # inspecting the data

### adding additional experimetal group information
samplenames <- colnames(x)
x$samples$group <- group
#x$samples$group

### adding in titles to the gene reference table
colnames(genes) <- c('ref_seqID', 'gene', 'symbol')
rownames(genes) <- genes$ref_seqID
#head(genes) # inspecting the data

########### writing out a data matrix for the user to download #############

df <- data.frame(x$counts)
GeneID <- genes[rownames(x),3]
df <- cbind(df, GeneID)
write.table(df, file = "results/tables/1.data_matrix.txt")


#################### data filtering #######################

### changing our FPKM values to Log2(fpkm) and setting log2(0) to 0
x$counts <- log2(x$counts)
x$counts[is.infinite(x$counts)] <- 0
# head(x$counts) # inspecting the data

table(rowSums(x$counts==0)==9)  # this shows us the genes with no expression in any samples
# Genes must be expressed in at least one group or in >3 samples to be kept for downstream analysis.
# dim(x) # before filtering (18907)
lfpkm <- x$counts
keep.exprs <- rowSums(lfpkm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
# dim(x) # after filtering (14021)


################## Unsupervised exploratory data analysis ##################
##########PCA & HCA on total data set (ie BEFORE DGE-analysis)##############

### PCA plot
pdf('results/graphs/1.PCA_all_data.pdf')
lfpkm <- t(x$counts)
Xpca <-prcomp(lfpkm)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)
dev.off()


### HCA dendrogram
pdf('results/graphs/2.HC_all_data.pdf')
hc <-hclust(dist(lfpkm))
dend <-as.dendrogram(hc)
colorCodes <-c(Uninfected = "black", Hendra_8hrs  = "red", Hendra_24hrs = "green")
labels(dend) <- group[order.dendrogram(dend)]
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend)
dev.off()


############ Differntial Gene Expression Analysis ###############

# building the design matrix (ie setting sample groups)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

# building the contrast matrix 
contr.matrix <- makeContrasts(
  Uninf_vs_H_8hrs = Uninfected - Hendra_8hrs, 
  Uninf_vs_H_24hrs = Uninfected - Hendra_24hrs, 
  H_8hrs_vs_H_24hrs = Hendra_8hrs - Hendra_24hrs, 
  levels = colnames(design))
contr.matrix

##### we will start here swapping v for our fpkm table
fit <- lmFit(x$counts, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit, trend=TRUE)    # we will use trend=TRUE here
summary(decideTests(efit))    # this shows if genes are up/downegulated or unchanged
dt <- decideTests(efit)   

#### this selects those genes that are differentially expressed from each comparrision
DE_UvsT8 <- which(dt[,1]!=0)
DE_UvsT24 <- which(dt[,2]!=0)
DE_T8vsT24 <- which(dt[,3]!=0)
DE_all <- c(which(dt[,1]!=0), which(dt[,2]!=0) , which(dt[,3]!=0))
DE_all <- unique(DE_all)

#### this looks up their gene names
de_UvsT8 <- genes[DE_UvsT8,]
de_UvsT24 <- genes[DE_UvsT24,]
de_T8vsT24 <- genes[DE_T8vsT24,]
de_all <- genes[DE_all,]

### this writes out the list of DE genes for each comparision
write.table(de_UvsT8, file="results/tables/2.de_UvsT8.txt", row.names=FALSE)
write.table(de_UvsT24, file="results/tables/3.de_UvsT24.txt", row.names=FALSE)
write.table(de_T8vsT24, file="results/tables/4.de_T8vsT24.txt", row.names=FALSE)
write.table(de_all, file="results/tables/5.de_all.txt", row.names=FALSE)

### this writes out a table describing the limma derived model
write.fit(efit, dt, file="results/tables/6.de_limma_model.txt")


######################### Volcano Plots ######################
# 1) make the plot of logFC vs -log10(P.Value)
# 2) highlight those genes with adj.P.Val < 0.05 & |logFC|>2) in red
# 3) label significant genes using the genes lookup table


### Examining individual DE genes - topTreat displays they in descending order of significance
Uninf_vs_H_8hrs <- topTreat(efit, coef=1, n=Inf)
Uninf_vs_H_24hrs <- topTreat(efit, coef=2, n=Inf)
H_8hrs_vs_H_24hrs <- topTreat(efit, coef=3, n=Inf)
# head(Uninf_vs_H_8hrs)    # - inspecting the data 
# head(Uninf_vs_H_24hrs)   # - inspecting the data 
# head(H_8hrs_vs_H_24hrs)  # - inspecting the data 


#### Uninfected vs Hendra_8hrs
pdf('results/graphs/3.Uninf_vs_H_8hrs_Volcano.pdf')
with(Uninf_vs_H_8hrs, plot(logFC, -log10(P.Value), pch=20, main="Uninfected vs Hendra_8hrs",
                                xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
with(subset(Uninf_vs_H_8hrs, adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                                pch=20, , cex =0.5,col="red"))
with(subset(Uninf_vs_H_8hrs, adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                                labs=genes[rownames(subset(Uninf_vs_H_8hrs,
                                adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.8))
dev.off()

#### Uninfected vs Hendra_24hrs
pdf('results/graphs/4.Uninf_vs_H_24hrs_Volcano.pdf')
with(Uninf_vs_H_24hrs, plot(logFC, -log10(P.Value), pch=20, main="Uninfected vs Hendra_24hrs",
                                xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
with(subset(Uninf_vs_H_24hrs, adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                                pch=20, , cex =0.5,col="red"))
with(subset(Uninf_vs_H_24hrs, adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                                labs=genes[rownames(subset(Uninf_vs_H_24hrs,
                                adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.4))
dev.off()

#### Hendra_8hrs vs Hendra_24hrs
pdf('results/graphs/5.H_8hrs_vs_H_24hrs_Volcano.pdf')
with(H_8hrs_vs_H_24hrs, plot(logFC, -log10(P.Value), pch=20, main="Hendra_8hrs vs Hendra_24hrs",
                                     xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
with(subset(H_8hrs_vs_H_24hrs, adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                                     pch=20, , cex =0.5,col="red"))
with(subset(H_8hrs_vs_H_24hrs, adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                                     labs=genes[rownames(subset(H_8hrs_vs_H_24hrs,
                                     adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.4))
dev.off()

############ Plotting a Heatmap of DE_Genes ###########

# the x$count data matrix and genes lookup table are subset by the DE_all genes list
pdf('results/graphs/6.Topgenes_Heatmap.pdf')
heatmap.2(x$counts[DE_all[1:60],], col=redgreen(30), scale="row", 
          labRow=genes[DE_all[1:60],3], labCol=group,
          key=T, keysize=1, density.info="none",
          trace= "none", cexCol=0.9, cexRow = 0.5, 
          margin=c(6,9))
dev.off()

################## Unsupervised exploratory data analysis ##################
############### PCA & HCA on genes identified by DGE-analysis ##############

### PCA plot
pdf('results/graphs/7.PCA_topgenes.pdf')
lfpkm <- t(x$counts[DE_all,])
Xpca <-prcomp(lfpkm, scale= TRUE)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)
dev.off()

### HCA dendrogram
pdf('results/graphs/8.HC_topgenes.pdf')
hc <-hclust(dist(lfpkm))
dend <-as.dendrogram(hc)
colorCodes <-c(Basal = "black", LP = "red", ML = "green")
labels(dend) <- group[order.dendrogram(dend)]
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend)
dev.off()

sessionInfo()

