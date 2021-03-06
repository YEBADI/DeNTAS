##### RNA-seq analysis of FPKM data using limma & edgeR

rm (list =ls())  # clear the workkspace



######### loading the required packages ###############



# installing the required packages if needed
#install.packages('limma')
#install.packages('edgeR')
#install.packages('dendextend')
#install.packages('gplots')
#install.packages('calibrate')

library(limma)
library(edgeR)
library(dendextend)
library(gplots)
library(calibrate)



########## Reading in the variables set by python#######


myArgs <- commandArgs(trailingOnly = TRUE) 
myArgs       

genes  <- read.table(file=myArgs[1])   # reads in 1_gene_list.txt used to look up gene names by reference ID
myArgs <- myArgs[2:length(myArgs)]     # remove 1_gene_list.txt from the mrArgs list

x <- length(myArgs)                          # the first half of the myArgs represent the sample file names              
files  <- c(myArgs[1:(x/2)])                 # the second half of the myArgs represent their experimental groups
group  <- as.factor(c(myArgs[((x/2)+1):x]))

group_number<- length(unique(group))   # can be 2, 3 or 4 used downstream to dictate analysis options

group_IDs <- unique(group)             # extracting the name of each group
group_IDs
group_1 <- toString(group_IDs[1])
group_2 <- toString(group_IDs[2])
group_3 <- toString(group_IDs[3])
group_4 <- toString(group_IDs[4])



################ loading in the data ####################



### readDGE function creates a DGEList-object containing 1) sample info  & 2) expression matrix 
x <- readDGE(files, columns=c(1,3))
# class(x) # inspecting the data
# dim(x) # inspecting the data

### adding additional experimetal group information
samplenames <- colnames(x)
x$samples$group <- group

### adding in titles to the gene reference table
colnames(genes) <- c('ref_seqID', 'gene', 'symbol')
rownames(genes) <- genes$ref_seqID
#head(genes) # inspecting the data



## writing out a data matrix for the user to download ##

df <- data.frame(x$counts)
GeneID <- genes[rownames(x),3]
df <- cbind(GeneID, df)
write.table(df, file = "static/results/tables/1.data_matrix.tsv",  quote=FALSE, sep='\t', row.names=FALSE)



#################### data filtering #######################



### changing our FPKM values to Log2(fpkm) and setting log2(0) to 0
x$counts <- log2(x$counts)
x$counts[is.infinite(x$counts)] <- 0
head(x$counts) # inspecting the data


# Genes must be expressed in at least one group or in >3 samples to be kept for downstream analysis.
#dim(x) # before filtering (~18907)
lfpkm <- x$counts
keep.exprs <- rowSums(lfpkm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
#dim(x) # after filtering (~14021)



########## Unsupervised exploratory data analysis #############
####PCA & HCA on total data set (ie BEFORE DGE-analysis)#######



### PCA plot

pdf('static/results/graphs/1.PCA_all_data.pdf')
lfpkm <- t(x$counts)
Xpca <-prcomp(lfpkm)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)
dev.off()

### HCA dendrogram
pdf('static/results/graphs/2.HC_all_data.pdf')
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


### The following code builds a contrast matrix with 1, 3 or 6
### comaprisions depending on the number of experimental groups (group_number)

if (group_number ==2){
  contr_1 <- paste(c(group_1,"-",group_2), sep = "", collapse = " ") 
  
  contr.matrix <- makeContrasts(
    contr_1,
    levels = colnames(design))
  contr.matrix
  
} else if (group_number == 3){
  contr_1 <- paste(c(group_1,"-",group_2), sep = "", collapse = " ")
  contr_2 <- paste(c(group_1,"-",group_3), sep = "", collapse = " ")
  contr_3 <- paste(c(group_2,"-",group_3), sep = "", collapse = " ")
  
  contr.matrix <- makeContrasts(
    contr_1, contr_2, contr_3, 
    levels = colnames(design))
  contr.matrix
  
} else if (group_number == 4){
  contr_1 <- paste(c(group_1,"-",group_2), sep = "", collapse = " ")
  contr_2 <- paste(c(group_1,"-",group_3), sep = "", collapse = " ")
  contr_3 <- paste(c(group_1,"-",group_4), sep = "", collapse = " ")
  contr_4 <- paste(c(group_2,"-",group_3), sep = "", collapse = " ")
  contr_5 <- paste(c(group_2,"-",group_4), sep = "", collapse = " ")
  contr_6 <- paste(c(group_3,"-",group_4), sep = "", collapse = " ")
  
  contr.matrix <- makeContrasts(
    contr_1, contr_2, contr_3, contr_4, contr_5, contr_6,
    levels = colnames(design))
  contr.matrix
}

##### Using limma to perform DGE analysis according to the ######
############## set design and contrast matrix ###################

fit <- lmFit(x$counts, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit, trend=TRUE)    

#summary(decideTests(efit))    # this shows if genes are up/downegulated or unchanged
dt <- decideTests(efit)   

### The following code pulls out the refID of those gense determined to
### be differentially expressed by the limma model

if (group_number ==2){
  DE_grp1_vs_grp2 <- rownames(as.matrix(which(dt[,1]!=0)))
  DE_all <- DE_grp1_vs_grp2
  
} else if (group_number == 3){
  DE_grp1_vs_grp2 <- rownames(as.matrix(which(dt[,1]!=0)))
  DE_grp1_vs_grp3 <- rownames(as.matrix(which(dt[,2]!=0)))
  DE_grp2_vs_grp3 <- rownames(as.matrix(which(dt[,3]!=0)))
  DE_all <- c(DE_grp1_vs_grp2 , DE_grp1_vs_grp3, DE_grp2_vs_grp3)
  DE_all <- unique(DE_all)
  
} else if (group_number == 4){
  DE_grp1_vs_grp2 <- rownames(as.matrix(which(dt[,1]!=0)))
  DE_grp1_vs_grp3 <- rownames(as.matrix(which(dt[,2]!=0)))
  DE_grp1_vs_grp4 <- rownames(as.matrix(which(dt[,3]!=0)))
  DE_grp2_vs_grp3 <- rownames(as.matrix(which(dt[,4]!=0)))
  DE_grp2_vs_grp4 <- rownames(as.matrix(which(dt[,5]!=0)))
  DE_grp3_vs_grp4 <- rownames(as.matrix(which(dt[,6]!=0)))
  DE_all <- c(DE_grp1_vs_grp2 , DE_grp1_vs_grp3, DE_grp1_vs_grp4, DE_grp2_vs_grp3, DE_grp2_vs_grp4,DE_grp3_vs_grp4 )
  DE_all <- unique(DE_all)
}

de_all <- genes[DE_all,] # this looks up their gene names from the refID


### here we are forming a merged de_all + evalue table for visualisations
e_values <- efit$p.value[DE_all,]
merged <- merge(de_all,e_values, by.x=1, by.y=0)
write.table(merged, file="static/results/tables/6.de_all_with_evalues.tsv", quote=FALSE, sep='\t', row.names=FALSE)



######################### Volcano Plots ######################
# 1) make the plot of logFC vs -log10(P.Value)
# 2) highlight those genes with adj.P.Val < 0.05 & |logFC|>2) in red
# 3) label significant genes using the genes lookup table


### Examining individual DE genes - topTreat displays them in descending order of significance
### here we evaluate the different coeffiects of the topTreat on efit depending on group_number

if (group_number ==2){
  topT_grp1_vs_grp2 <- topTreat(efit, coef=1, n=Inf)
  topT_list <- list(topT_grp1_vs_grp2)
  
} else if (group_number == 3){
  topT_grp1_vs_grp2 <- topTreat(efit, coef=1, n=Inf)  
  topT_grp1_vs_grp3 <- topTreat(efit, coef=2, n=Inf)
  topT_grp2_vs_grp3 <- topTreat(efit, coef=3, n=Inf)
  topT_list <- list(topT_grp1_vs_grp2, topT_grp1_vs_grp3, topT_grp2_vs_grp3)
  
} else if (group_number == 4){
  topT_grp1_vs_grp2 <- topTreat(efit, coef=1, n=Inf)  
  topT_grp1_vs_grp3 <- topTreat(efit, coef=2, n=Inf)
  topT_grp1_vs_grp4 <- topTreat(efit, coef=3, n=Inf)
  topT_grp2_vs_grp3 <- topTreat(efit, coef=4, n=Inf)  
  topT_grp2_vs_grp4 <- topTreat(efit, coef=5, n=Inf)
  topT_grp3_vs_grp4 <- topTreat(efit, coef=6, n=Inf)
  topT_list <- list(topT_grp1_vs_grp2, topT_grp1_vs_grp3, topT_grp1_vs_grp4, topT_grp2_vs_grp3, topT_grp2_vs_grp4, topT_grp3_vs_grp4)

}

# Generatng the comparison names for Volcano plot titles
if (group_number ==2){
  comp_1 <- paste(c(group_1,"vs",group_2), sep = "", collapse = " ")
  comp_list <- list(comp_1)
  
} else if (group_number == 3){
  comp_1 <- paste(c(group_1,"vs",group_2), sep = "", collapse = " ")
  comp_2 <- paste(c(group_1,"vs",group_3), sep = "", collapse = " ")
  comp_3 <- paste(c(group_2,"vs",group_3), sep = "", collapse = " ")
  comp_list <- list(comp_1, comp_2, comp_3)

} else if (group_number == 4){
  comp_1 <- paste(c(group_1,"vs",group_2), sep = "", collapse = " ")
  comp_2 <- paste(c(group_1,"vs",group_3), sep = "", collapse = " ")
  comp_3 <- paste(c(group_1,"vs",group_4), sep = "", collapse = " ")
  comp_4 <- paste(c(group_2,"vs",group_3), sep = "", collapse = " ")
  comp_5 <- paste(c(group_2,"vs",group_4), sep = "", collapse = " ")
  comp_6 <- paste(c(group_3,"vs",group_4), sep = "", collapse = " ")
  comp_list <- list(comp_1, comp_2, comp_3, comp_4, comp_5, comp_6)
}


graph_names=list()
for (i in 1:length(topT_list)){
  graph_names <- c(graph_names, (paste(c('static/results/graphs/Volcano_',i,'.pdf'), sep = "", collapse = "")))
}

#For loop used to plot a titled Volcano plot for each experimental comparison

for (i in 1:length(topT_list)){
  pdf(graph_names[[i]])
  with(topT_list[[i]], plot(logFC, -log10(P.Value), pch=20, main=comp_list[[i]],
                            xlim=c(-8,8), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
  text(-6.5,-log10(0.05),col="red","P.Value = 0.05", pos = 3, cex=0.8)
  with(subset(topT_list[[i]], adj.P.Val  < 0.05 & abs(logFC)>2), points(logFC, -log10(P.Value),
                            pch=20, , cex =0.5,col="red"))
  with(subset(topT_list[[i]], adj.P.Val  < 0.05 & abs(logFC)>2), textxy(logFC, -log10(P.Value), 
                            labs=genes[rownames(subset(topT_list[[i]],
                            adj.P.Val  < 0.05 & abs(logFC)>2)),3], cex=.4))
  dev.off()
}



############ Plotting a Heatmap of DE_Genes ###########



# the x$count data matrix and genes lookup table are subset by the DE_all genes list
pdf('static/results/graphs/6.Topgenes_Heatmap.pdf')
heatmap.2(x$counts[DE_all[1:100],], col=redgreen(30), scale="row", 
          labRow=genes[DE_all[1:100],3], labCol=group,
          key=T, keysize=1, density.info="none",
          trace= "none", cexCol=0.9, cexRow = 0.5, 
          margin=c(6,9))
dev.off()



################## Unsupervised exploratory data analysis ##################
############### PCA & HCA on genes identified by DGE-analysis ##############



### PCA plot
pdf('static/results/graphs/7.PCA_topgenes.pdf')
lfpkm <- t(x$counts[DE_all,])
Xpca <-prcomp(lfpkm, scale= TRUE)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)
dev.off()



### HCA dendrogram
pdf('static/results/graphs/8.HC_topgenes.pdf')
hc <-hclust(dist(lfpkm))
dend <-as.dendrogram(hc)
colorCodes <-c(group_1 = "black", group_2  = "red", group_3 = "green")
labels(dend) <- group[order.dendrogram(dend)]
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend)
dev.off()

rm (list =ls())              # clearing the workspace
z= "analysis complete :) "
z
#sessionInfo()