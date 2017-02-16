# Could have raw_input from user for number of columns and use this object throughout so it overcomes the 3 groups barrier?
#alternatively if this proves to be too difficult; to meet the deadline we can have a set of hardcoded R scripts for various group numbers
#and call the appropriate script according to the number of groups the end-user picks in FLASK

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
myArgs
# useful for working
# myArgs <- c("blast/1_gene_list.txt", "blast/blast_10s0r1.fasta.txt", "blast/blast_11s0r2.fasta.txt",
#              "blast/blast_12s0r3.fasta.txt", "blast/blast_13s8r1.fasta.txt", "blast/blast_14s8r2.fasta.txt" ,
#              "blast/blast_15s8r3.fasta.txt", "blast/blast_16s24r1.fasta.txt", "blast/blast_17s24r2.fasta.txt" ,
#              "blast/blast_18s24r3.fasta.txt", "Uninfected", "Uninfected", "Uninfected", "Hendra_8hrs",
#              "Hendra_8hrs", "Hendra_8hrs", "Hendra_24hrs", "Hendra_24hrs", "Hendra_24hrs")
# # myArgs
genes  <- read.table(file=myArgs[1], fill = TRUE)
myArgs <- myArgs[2:length(myArgs)]

x <- length(myArgs)
files  <- c(myArgs[1:(x/2)])
group  <- as.factor(c(myArgs[((x/2)+1):x]))

# extracting variable names for downstream use

group_IDs <- unique(group)
group_1 <- toString(group_IDs[1])
group_2 <- toString(group_IDs[2])
group_3 <- toString(group_IDs[3])

################ loading in the data ####################

### readDGE function creates a DGEList-object containing 1) sample info  & 2) expression matrix 
x <- readDGE(files, columns=c(1,3))
# class(x) # inspecting the data
# dim(x) # inspecting the data
# inspecting the data

### adding additional experimetal group information
samplenames <- colnames(x)
x$samples$group <- group

### adding in titles to the gene reference table
colnames(genes) <- c('ref_seqID', 'gene', 'symbol')
rownames(genes) <- genes$ref_seqID
#head(genes) # inspecting the data

########### writing out a data matrix for the user to download #############

df <- data.frame(x$counts)
GeneID <- genes[rownames(x),3]
df <- cbind(df, GeneID)
write.table(df, file = "static/results/tables/1.data_matrix.tsv",  quote=FALSE, sep='\t', row.names=FALSE)

#################### data filtering #######################

### changing our FPKM values to Log2(fpkm) and setting log2(0) to 0
x$counts <- log2(x$counts)
x$counts[is.infinite(x$counts)] <- 0
# head(x$counts) # inspecting the data

#table(rowSums(x$counts==0)==9)  # this shows us the genes with no expression in any samples
# Genes must be expressed in at least one group or in >3 samples to be kept for downstream analysis.
# dim(x) # before filtering (18907)
lfpkm <- x$counts
keep.exprs <- rowSums(lfpkm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
# dim(x) # after filtering (14021)


################## Unsupervised exploratory data analysis ##################
##########PCA & HCA on total data set (ie BEFORE DGE-analysis)##############

### PCA plot
png('static/results/graphs/1.PCA_all_data.png')
lfpkm <- t(x$counts)
Xpca <-prcomp(lfpkm)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)
dev.off()


### HCA dendrogram
png('static/results/graphs/2.HC_all_data.png')
hc <-hclust(dist(lfpkm))
dend <-as.dendrogram(hc)
colorCodes <-c(group_1 = "black", group_2  = "red", group_3 = "green")
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
contr_1 <- paste(c(group_1,"-",group_2), sep = "", collapse = " ")
contr_2 <- paste(c(group_1,"-",group_3), sep = "", collapse = " ")
contr_3 <- paste(c(group_2,"-",group_3), sep = "", collapse = " ")

contr.matrix <- makeContrasts(
  contr_1,
  contr_2, 
  contr_3, 
  levels = colnames(design))
contr.matrix

# delete if everything works
#contr.matrix <- makeContrasts(
#  Uninf_vs_H_8hrs = Uninfected - Hendra_8hrs, 
#  Uninf_vs_H_24hrs = Uninfected - Hendra_24hrs, 
#  H_8hrs_vs_H_24hrs = Hendra_8hrs - Hendra_24hrs, 
#  levels = colnames(design))
#contr.matrix

##### we will start here swapping v for our fpkm table
fit <- lmFit(x$counts, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit, trend=TRUE)    # we will use trend=TRUE here
efit
#summary(decideTests(efit))    # this shows if genes are up/downegulated or unchanged
dt <- decideTests(efit)   

#### this selects those genes that are differentially expressed from each comparrision
#DE_UvsT8 <- rownames(as.matrix(which(dt[,1]!=0)))
#DE_UvsT24 <- rownames(as.matrix(which(dt[,2]!=0)))
#DE_T8vsT24 <- rownames(as.matrix(which(dt[,3]!=0)))
#DE_all <- c(DE_UvsT8,DE_UvsT24,DE_T8vsT24)
#DE_all <- unique(DE_all)
#length(DE_all)

DE_grp1_vs_grp2 <- rownames(as.matrix(which(dt[,1]!=0)))
DE_grp1_vs_grp3 <- rownames(as.matrix(which(dt[,2]!=0)))
DE_grp2_vs_grp3 <- rownames(as.matrix(which(dt[,3]!=0)))
DE_all <- c(DE_grp1_vs_grp2 , DE_grp1_vs_grp3, DE_grp2_vs_grp3)
DE_all <- unique(DE_all)


#### this looks up their gene names
de_grp1_vs_grp2 <- genes[DE_grp1_vs_grp2 ,]
de_grp1_vs_grp3 <- genes[DE_grp1_vs_grp3 ,]
de_grp2_vs_grp3 <- genes[DE_grp2_vs_grp3 ,]
de_all <- genes[DE_all,]

### here we are forming a merged de_all + evalue table for visualisations
e_values <- efit$p.value[DE_all,]
merged <- merge(de_all,e_values, by.x=1, by.y=0)


### this writes out the list of DE genes for each comparision
write.table(de_grp1_vs_grp2, file="static/results/tables/2.de_grp1_vs_grp2.tsv", quote=FALSE, sep='\t', row.names=FALSE)
write.table(de_grp1_vs_grp3, file="static/results/tables/3.de_grp1_vs_grp3.tsv", quote=FALSE, sep='\t', row.names=FALSE)
write.table(de_grp2_vs_grp3, file="static/results/tables/4.de_grp2_vs_grp3.tsv", quote=FALSE, sep='\t', row.names=FALSE)
write.table(de_all, file="static/results/tables/5.de_all.tsv", quote=FALSE, sep='\t', row.names=FALSE)
write.table(merged, file="static/results/tables/6.de_all_with_evalues.tsv", quote=FALSE, sep='\t', row.names=FALSE)

######################### Volcano Plots ######################
# 1) make the plot of logFC vs -log10(P.Value)
# 2) highlight those genes with adj.P.Val < 0.05 & |logFC|>2) in red
# 3) label significant genes using the genes lookup table


### Examining individual DE genes - topTreat displays they in descending order of significance
topT_grp1_vs_grp2 <- topTreat(efit, coef=1, n=Inf)
topT_grp1_vs_grp3 <- topTreat(efit, coef=2, n=Inf)
topT_grp2_vs_grp3 <- topTreat(efit, coef=3, n=Inf)


# comparison names for graph titles
comp_1 <- paste(c(group_1,"vs",group_2), sep = "", collapse = " ")
comp_2 <- paste(c(group_1,"vs",group_3), sep = "", collapse = " ")
comp_3 <- paste(c(group_2,"vs",group_3), sep = "", collapse = " ")

#### Uninfected vs Hendra_8hrs
png('static/results/graphs/3.grp1_vs_grp2_Volcano.png')
with(topT_grp1_vs_grp2, plot(logFC, -log10(P.Value), pch=20, main=comp_1,
                                xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
with(subset(topT_grp1_vs_grp2, adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                                pch=20, , cex =0.5,col="red"))
with(subset(topT_grp1_vs_grp2, adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                                labs=genes[rownames(subset(topT_grp1_vs_grp2,
                                adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.8))
dev.off()

#### Uninfected vs Hendra_24hrs
png('static/results/graphs/4.grp1_vs_grp3_Volcano.png')
with(topT_grp1_vs_grp3, plot(logFC, -log10(P.Value), pch=20, main=comp_2,
                                xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
with(subset(topT_grp1_vs_grp3, adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                                pch=20, , cex =0.5,col="red"))
with(subset(topT_grp1_vs_grp3, adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                                labs=genes[rownames(subset(topT_grp1_vs_grp3,
                                adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.4))
dev.off()

#### Hendra_8hrs vs Hendra_24hrs
png('static/results/graphs/5.grp2_vs_grp3_Volcano.png')
with(topT_grp2_vs_grp3 , plot(logFC, -log10(P.Value), pch=20, main=comp_3,
                                     xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
with(subset(topT_grp2_vs_grp3 , adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                                     pch=20, , cex =0.5,col="red"))
with(subset(topT_grp2_vs_grp3 , adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                                     labs=genes[rownames(subset(topT_grp2_vs_grp3 ,
                                     adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.4))
dev.off()

############ Plotting a Heatmap of DE_Genes ###########

# the x$count data matrix and genes lookup table are subset by the DE_all genes list
png('static/results/graphs/6.Topgenes_Heatmap.png')
heatmap.2(x$counts[DE_all[1:60],], col=redgreen(30), scale="row", 
          labRow=genes[DE_all[1:60],3], labCol=group,
          key=T, keysize=1, density.info="none",
          trace= "none", cexCol=0.9, cexRow = 0.5, 
          margin=c(6,9))
dev.off()

################## Unsupervised exploratory data analysis ##################
############### PCA & HCA on genes identified by DGE-analysis ##############

### PCA plot
png('static/results/graphs/7.PCA_topgenes.png')
lfpkm <- t(x$counts[DE_all,])
Xpca <-prcomp(lfpkm, scale= TRUE)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)
dev.off()

### HCA dendrogram
png('static/results/graphs/8.HC_topgenes.png')
hc <-hclust(dist(lfpkm))
dend <-as.dendrogram(hc)
colorCodes <-c(group_1 = "black", group_2  = "red", group_3 = "green")
labels(dend) <- group[order.dendrogram(dend)]
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend)
dev.off()

rm (list =ls()) 
z= "analysis complete :) "
z
#sessionInfo()
