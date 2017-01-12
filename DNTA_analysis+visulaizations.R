#### RNA-seq analysis using limma & edgeR
#### NB this analysis starts from raw counts & library sizes - therfore it will need to be adapted
#### to enable the input of our FPKM data - however DGE analysis and data visulaization will remian
#### approcimately constant

rm (list =ls())  # clear the workkspace

#### installing the required packages
#install.packages('limma')
#install.packages('edgeR')
#install.packages('dendextend')
#install.packages('gplots')

### loading the required packages
library(limma)
library(edgeR)
library(dendextend)
library(gplots)
library(RColorBrewer)

### downloading & unzipping the input files 
### Each of these text files contains the raw gene-level counts for a given sample.
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"

utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
R.utils::gunzip(i, overwrite=TRUE)

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")

read.delim(files[1], nrow=5) # this just checks what is the files

### readDGE function creates a DGEList-object containing 1) sample info  & 2) expression matrix 
### NB in this sample data it is Entrez ID, GeneLength & raw count

x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)

### adding additional information, ie experimetal group, batch number/date etc
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
x$samples$group

### here we are finding the gene names from their ENTREZID compared against the installed
### mouse annotations
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")

head(genes)
genes <- genes[!duplicated(genes$ENTREZID),] # removing duplicates
x$genes <- genes        # adding the gene information to the DGElist object

x  # check what the final DGElist object contains


#### Data pre-processing
# now converting reads to CPM (in our case we will just start with FPKM)
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

table(rowSums(x$counts==0)==9)   # this shows us the genes with no expression in any samples
                                 # we shouldn't have any of these


# Genes must be expressed in at least one group (or in at least three samples across ...
# ... the entire experiment) to be kept for downstream analysis.

dim(x) # before normalization (27179 rows)
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) # after normalization (14165 rows)

### graphically demonstrating the filtering using density plots
### - this would be a nice feature to return to our user 
par(mfrow = c(1, 2))
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", cex=0.75)

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", cex=0.75)


### normalizing the gene expression distributions
### it is assumed that all samples should have a similar range and distribution of expression values
### this can be altered by batch effects etc
### Normalisation ensures expression distributions are similar across the entire experiment.

par(mfrow = c(1, 1))
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
boxplot(lcpm, las=2, col=col, main="")    
title(main="Normalised expression distributions",ylab="Log-cpm")


### Unsupervised exploratory data analysis. PCA & HCA on total data set (ie BEFORE DGE-analysis)
### PCA plot
lcpm <- t(cpm(x, log=TRUE))
Xpca <-prcomp(lcpm, scale= TRUE)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)

### HCA dendrogram
hc <-hclust(dist(lcpm))
dend <-as.dendrogram(hc)
colorCodes <-c(Basal = "black", LP = "red", ML = "green")
labels(dend) <- group[order.dendrogram(dend)]
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend)

### differntial gen expression analysis
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

###### Here I use Voom (which starts from raw counts) - we will skip this bit
v <- voom(x, design)
v
##### we will start here swapping v for our fpkm table
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)    # we will use trend=TRUE here



summary(decideTests(efit))    # this shows if genes are up/downegulated or unchanged
tfit <- treat(vfit, lfc=1)    # The treat method is used to calculate p-values from empirical Bayes 
dt <- decideTests(tfit)       # moderated t-statistics with a minimum log-FC requirement.
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)  # this identifies those gene that are different between both basal & LP and basal & MP
length(de.common)

head(tfit$genes$SYMBOL[de.common], n=20)
write((tfit$genes$SYMBOL[de.common]), file="de_common.txt") # this writes out a list of the diffectially expressed genes 
write.fit(tfit, dt, file="results.txt")   # this saves a tsv txt file of the DGE reults

### Examining individual DE genes - topTreat displays they in descending order of significance
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)

### volcano plots
par(mfrow = c(1, 2))
with(basal.vs.lp, plot(basal.vs.lp$logFC, -log10(basal.vs.lp$P.Value), pch=20, 
          main="Basal vs LP", xlim=c(-10,10), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))
with(basal.vs.ml, plot(basal.vs.ml$logFC, -log10(basal.vs.lp$P.Value), pch=20, 
          main="Basal vs ml", xlim=c(-10,10), cex =0.25, abline(h=-log10(0.05),col="red",lty="44")))


### now ploting a heatmap

library(gplots)
topgenes <- c(basal.vs.lp$ENTREZID[1:40], basal.vs.ml$ENTREZID[1:40])
#length(topgenes)
topgenes <- unique(topgenes)
#length(topgenes)

i <- which(x$genes$ENTREZID %in% topgenes)
par(mfrow = c(1, 1))
heatmap.2(v$E[i,], col=redgreen(30), scale="row", 
          labRow=v$genes$SYMBOL[i], labCol=group,
          key=T, keysize=1, density.info="none",
          trace= "none", cexCol=0.9, dendrogram="column",
          margin=c(2,14))


### Unsupervised exploratory data analysis. PCA & HCA on genes identified by DGE-analysis
### PCA plot
lcpm <- t(cpm(x, log=TRUE))
lcpm <- lcpm[,i]
Xpca <-prcomp(lcpm, scale= TRUE)
Xscores   <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,
     cex.lab=0.7, cex.axis = 0.7, col=group,
     main="Principal components analysis")
text(Xscores, labels=group, cex= 1, pos=1)




### HCA dendrogram
hc <-hclust(dist(lcpm[,i]))
dend <-as.dendrogram(hc)
colorCodes <-c(Basal = "black", LP = "red", ML = "green")
labels(dend) <- group[order.dendrogram(dend)]
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend)


sessionInfo()

