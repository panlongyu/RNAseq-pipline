# plotting
library(DESeq2) # rna-seq
library(stringr)
library(rafalib) # nice plot arrangement
rafalib::mypar(mar=c(6,2.5,2.5,1)) #sets nice arrangement for the whole document

co <- read.delim("allSamples.featureCounts.tsv",sep="\t",header=TRUE,stringsAsFactors=F,comment.char="#")
cr <- co[,7:24]
rownames(cr) <- co$Geneid
colnames(cr) <- substr(colnames(cr),43,47)
group=c(rep("Cont",9),rep("NEC",9))
mr=as.data.frame(group)
row.names(mr)=colnames(cr)
all.equal(colnames(cr),rownames(mr))
rafalib::mypar(1,2,mar=c(6,3,3,2))
boxplot(log2(as.matrix(cr)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Raw data")
hist(log2(as.matrix(cr)+1),ylab="",las=2,main="Raw data")
par(mfrow=c(1,1))

barplot(colSums(cr>3),ylab="Number of detected genes",las=2) 
abline(h=median(colSums(cr>3)))
barplot(rowSums(cr>3),xlab="Genes",ylab="Number of samples",names.arg="")
abline(h=median(rowSums(cr>3)),col="red"

hist(rowSums(cr>3))
# remove genes with low counts
keep_genes <- rowSums( cr > 5 ) >= 3
cf <- cr[keep_genes,]
boxplot(log2(as.matrix(cf)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Filtered data")
hist(rowSums(cf>3))
all.equal(colnames(cf),rownames(mr))
write.csv(cf,"./counts_filtered.csv",quote=F)

#均一化
cc <- t( t(cf) / colSums(cf) * 1e6 )
cc <- log2( cc + 1 )
boxplot(cc,ylab=expression('Log'[2]~'Read counts'),las=2,main="CPM")

#TPM
#' @title Compute TPM from a read count matrix
#' @param counts A numeric data.frame of read counts with samples (columns) and genes (rows).
#' @param len A vector of gene cds length equal to number of rows of dfr.
#'
#' https://support.bioconductor.org/p/91218/
#'
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

g <- data.frame( ensembl_gene_id = co$Geneid , 
                 transcript_length = co$Length,
                 stringsAsFactors = F, row.names = co$Geneid)

g <- g[!duplicated(g$ensembl_gene_id),]

igenes <- intersect(rownames(cf),g$ensembl_gene_id)
g1 <- g[igenes,]
cf1 <- cf[igenes,]
all.equal(rownames(cf1),g1$ensembl_gene_id)

ct <- tpm(cf1,g1$transcript_length)
ct <- log2( ct + 1 )
boxplot(ct,ylab=expression('Log'[2]~'Read counts'),las=2,main="TPM")
write.csv(ct,"./counts_tpm.csv",quote=F)

library(DESeq2)
mr$Group <- factor(mr$Group)
d <- DESeqDataSetFromMatrix(countData=cf,colData=mr,design=~Group)
d <- DESeq2::estimateSizeFactors(d,type="ratio")
cd <- log2( counts(d,normalized=TRUE) + 1 ) 
saveRDS(cd,"gene_counts_normalised_deseq2.Rds")
boxplot(cd,ylab=expression('Log'[2]~'Read counts'),las=2,main="DESeq2")
library(DESeq2)
mr$Group <- factor(mr$Group)
d <- DESeqDataSetFromMatrix(countData=cf,colData=mr,design=~Group)
d <- DESeq2::estimateSizeFactors(d,type="ratio")
d <- DESeq2::estimateDispersions(d)
cv <- as.data.frame(assay(varianceStabilizingTransformation(d,blind=T)),check.names=F)
write.csv(cv,"./gene_counts_vst.csv",quote=FALSE)
boxplot(cv,ylab=expression('Log'[2]~'Read counts'),las=2,main="VST")
rowVar <- function(x) apply(x,1,var)
rafalib::mypar(mfrow=c(2,2))
plot(rowMeans(cc),rowVar(cc),xlab=expression('Log'[2]~'Mean count'),ylab=expression('Log'[2]~'Variance'),main="CPM",cex=.1)
plot(rowMeans(ct),rowVar(ct),xlab=expression('Log'[2]~'Mean count'),ylab=expression('Log'[2]~'Variance'),main="TPM",cex=.1)
plot(rowMeans(cd),rowVar(cd),xlab=expression('Log'[2]~'Mean count'),ylab=expression('Log'[2]~'Variance'),main="DESeq2",cex=.1)
plot(rowMeans(cv),rowVar(cv),xlab=expression('Log'[2]~'Mean count'),ylab=expression('Log'[2]~'Variance'),main="VST",cex=.1)
rafalib::mypar(mar=c(6,2.5,2.5,1))

write.csv(ct,"./counts_tpm.csv",quote=F)

library(biomaRt)
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
ds <- as.data.frame(listDatasets(mart=mart))

# find all rows in dataset 'ds' where column 'description' contains the string 'mouse'
ds %>% filter(grepl("human",tolower(description)))
#ds %>% filter(grepl("mouse",tolower(description)))

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset(mart=mart,dataset="hsapiens_gene_ensembl")
#mart <- useDataset(mart=mart,dataset="mmusculus_gene_ensembl")
la <- listAttributes(mart=mart)
head(la)
searchAttributes(mart=mart,pattern="entrez")
#创建感兴趣的list
myattributes <- c("ensembl_gene_id",
                  "entrezgene_id",
                  "external_gene_name",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "strand",
                  "gene_biotype",
                  "description")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset(mart=mart,dataset="hsapiens_gene_ensembl")
#mart <- useDataset(mart=mart,dataset="mmusculus_gene_ensembl")
bdata <- getBM(mart=mart,attributes=myattributes,uniqueRows=T,
               useCache=FALSE)
head(bdata)
sum(duplicated(bdata$ensembl_gene_id))
sum(duplicated(bdata$entrezgene_id))
sum(duplicated(bdata$external_gene_name))

# arrange table by chr name and start position
bdata <- dplyr::arrange(bdata,chromosome_name,start_position)
write.table(bdata,"./human_genes.txt",sep="\t",dec=".",row.names=FALSE,quote=FALSE)
#转录本
mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#mart <- useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl")
t2g <- getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","external_gene_name"),mart=mart,useCache=FALSE)
write.table(t2g,"./human_transcripts.txt",sep="\t",dec=".",row.names=F,quote=F)

#与GO的关系
mart <- biomaRt::useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#mart <- biomaRt::useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl")
la <- listAttributes(mart=mart)
# find all rows in dataset 'lf' where column 'name' contains the string 'go'
head(la[grepl("go",tolower(la$name)),])

mart <- biomaRt::useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#mart <- biomaRt::useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl")
bdata <- getBM(mart=mart,attributes=c("entrezgene_id","go_id","go_linkage_type"),uniqueRows=T,useCache=FALSE)
write.table(bdata,"./human_go.txt",sep="\t",dec=".",row.names=FALSE,quote=FALSE)
#ID转换(人鼠转化)
mouse_genes <- c("ENSMUSG00000035847","ENSMUSG00000000214")
mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
getLDS(attributes=c("ensembl_gene_id"),filters="ensembl_gene_id",values=mouse_genes,mart=mouse, attributesL=c("external_gene_name"),martL=human,valuesL="external_gene_name",uniqueRows=F)


library(pheatmap) #plot heatmap
library(DESeq2) #differential gene expression for RNA-seq
library(pvclust) #clustering bootstrapping
library(biomaRt) #gene annotation
library(rafalib) #nice plot arrangement
rafalib::mypar(mar=c(6,2.5,2.5,1)) #sets nice arrangement for the whole document

# Load data and metadata
download_data("data/counts_tpm.csv")
data <- read.csv("data/counts_tpm.csv",header=TRUE,stringsAsFactors=FALSE,row.names=1)

#Convert ensembl IDs to gene names

download_data("data/mouse_genes.txt")
#mouse <- biomaRt::useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
human <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
annot <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"),mart = human,useCache=FALSE)
gene_names <- as.character (annot[match(rownames(data),annot[,"ensembl_gene_id"]),"external_gene_name"] )
gene_names[is.na(gene_names) ] <- ""
#删除丢失或融合的重复注释，并将转录本汇总成基因
data <- rowsum(data, group = as.character(gene_names) ) #assign gene names and sum duplicates
data <- data[ rownames(data) != "" , ] #remove blank annotation
data <- data[ rownames(data) != "NA" , ] #remove blank annotation
data <- data[ order(rownames(data)),]  #order genes alphabetically
dim(data)

write.csv(data ,"./gene_tpm_counts.csv",row.names=T)

cf <- read.csv("./counts_filtered.csv",header=TRUE,stringsAsFactors=FALSE,row.names=1)
cf <- rowsum(cf, group = as.character(gene_names) ) #assign gene names and sum duplicates
cf <- cf[ rownames(cf) != "" , ] #remove blank annotation
cf <- cf[ rownames(cf) != "NA" , ] #remove blank annotation
cf <- cf[ order(rownames(cf)),]  #order genes alphabetically
write.csv(cf ,"./gene_counts.csv",row.names=T)


#PCA Z分数归一化

Znorm <- t(apply(data,1,function(x) scale(x,center=T,scale=T)))
colnames(Znorm) <- colnames(data)

{
  mypar(1,2,mar=c(5,3,2,1))
  boxplot(t(data[1:30,]),ylim=c(0,12),las=2,col="grey",main="data raw_data",cex=.2)
  boxplot(t(Znorm[1:30,]),ylim=c(-4,4),las=2,col="grey",main="Z-score data",cex=.2)
  abline(h=0,col="red",lty=2)
}

#选择基因进行PCA

#Compute the mean, variance and cv for each gene and sort them in decreasing order
gene_stats <- data.frame(row.names=rownames(data))
means <- apply(data,1,mean)
vars <- apply(data,1,var)

#Plot the row means versus row means
{
  mypar(1,2,mar=c(5,3,2,1))
  plot(means,vars,cex=.1)
}

#Sort and select the top 500 highly variable genes from the data
vars <- sort(vars,decreasing=T)
top_var <- names(vars)[1:100]
boxplot(t(data[top_var[1:15],]),ylim=c(0,12),las=2,col="grey", main="data (top var genes)",cex=.2)

{
  mypar(1,2,mar = c(3,3,2,1))
  PC <-  prcomp( t( Znorm[ top_var, ]) ) #Method1
  PC <-  prcomp( t( data[ top_var, ]), center = TRUE, scale. = TRUE) #Method2
  
  mypar(1,2)
  plot(PC$x[,1],PC$x[,2],cex=2,col=factor(mr$Group),xlab="PC1",ylab="PC2",pch=16,main="PCA",las=1)
  text(PC$x[,1],PC$x[,2],cex=.7,labels = paste0(rownames(mr)),pos=3)
  
  plot(PC$x[,3],PC$x[,4],cex=2,col=factor(mr$Group),xlab="PC3",ylab="PC4",pch=16,main="PCA",las=1)
  text(PC$x[,3],PC$x[,4],cex=.7,labels = paste0(rownames(mr)),pos=3)
}

#计算PC方差
PC_sd <- setNames(PC$sdev,paste0("PC",1:length(PC$sdev)))
PC_var_expl <- (PC_sd^2)/sum(PC_sd^2)*100

{
  mypar(1,1)
  barplot(PC_var_expl,las=2,ylab="% variance explained")
  abline(h=10,lty=2)
}

#导致差异的Leading genes

leading_genes <- PC$rotation
head(leading_genes)

leading_PC1 <- sort(leading_genes[,1],decreasing=T)
leading_PC2 <- sort(leading_genes[,2],decreasing=T)

{
  mypar(1,2,mar=c(3,4,2,2))
  barplot(leading_PC1[15:1],las=2,horiz=T,cex.names=.8,cex.axis=0.8,yaxs="i")
  abline(v=0,lwd=2)
  barplot(leading_PC2[15:1],las=2,horiz=T,cex.names=.8,cex.axis=0.8,yaxs="i")
  abline(v=0,lwd=2)
}

#层次聚类,Distance between samples,The base R stats package already contains a function dist() that calculates distances between all pairs of samples.
d <- dist( t(data) , method="euclidean")
d
#相关性，We can first compute sample correlations using the cor() function. As you might already know, correlation range from -1 to 1, where 1 indicates that two samples are closest, -1 indicates that two samples are the furthest and 0 is somewhat in between.
#Compute sample correlations
sample_cor <- cor( data )
round(sample_cor,4)
pheatmap(sample_cor)

#Transform the scale from correlations
cor_distance <- -(sample_cor-1)/2
round(cor_distance,4)
# 利用colorRampPalette生成梯度渐变颜色，参考?colorRampPalette
first <- colorRampPalette(c("red", "blue"))(10)
second <- colorRampPalette(c("blue", "black"))(40)
palette <- c(first, second)
pheatmap(cor_distance,color=palette)

#Convert it to a distance object
d2 <- as.dist(cor_distance)
d2

#Clustering samples,After having calculated the distances between samples calculated, we can now proceed with the hierarchical clustering per-se. We will use the function hclust() for this purpose, in which we can simply run it with the distance objects created above.

#Clustering using euclidean distance
{
  mypar(1,2,mar=c(6,4,2,1))
  h <- hclust(d,method="complete")
  plot( as.dendrogram(h),las=1,main="d=euclidean\nh=complete")
  points(1:ncol(data),rep(0,ncol(data)),pch=16,cex=2,col=mr$Group[h$order])
}


h2 <- hclust(d2,method="complete")
{
  plot( as.dendrogram(h2),las=1, main="d=correlation\nh=complete")
  points(1:ncol(data),rep(0,ncol(data)),pch=16,cex=2, col=mr$Group[h2$order])
}

#Defining clusters,ince it is time consuming to cluster on all genes, this step is usually done only on the differentially expressed gene list (say up to ~3000 genes). Here, we can for instance cluster on the genes with highest variance top_var.

gene_cor  <- cor(t(Znorm[top_var, ]))
gene_dist <- as.dist(-(gene_cor-1)/2)
gene_clus <- hclust(gene_dist,method="complete")

HEIGHT <- 0.9
gene_clusters <- cutree(gene_clus,h=HEIGHT)
gene_clusters

{
  mypar(1,1,mar=c(6,4,2,1))
  plot( as.dendrogram(gene_clus),las=1,main="d=correlation\nh=complete")
  
  rect.hclust(gene_clus,h=HEIGHT)
  abline(h=HEIGHT,col="red",lty=2)
  points(1:length(gene_clusters),rep(0,length(gene_clusters)),pch=16,cex=2, col=factor(gene_clusters)[gene_clus$order])
  legend("topright",levels(factor(gene_clusters)),pch=16,col=  factor(levels(factor(gene_clusters))))
}

pheatmap( data[top_var,] , scale="row" , color = colorRampPalette(c("navy","white","firebrick"))(90), border_color=NA, cluster_cols=F)

#3.5 Clustering bootstrapping (optional)

# Clustering on our dataset (this step is slow)
#One way to measure clustering robustness / accuracy is by selecting part of the data set (say 90% of the genes), performing the clustering and recording which samples fall together. 
#pvclust paper: https://academic.oup.com/bioinformatics/article/22/12/1540/207339
#AU (Approximately Unbiased) p-value and BP (BootstrapProbability) value. AU p-value, which is computed by multiscale bootstrap resampling, is a better approximation to unbiased p-value than BP value computed by normal bootstrap resampling
pvc <- pvclust( data = data[top_var[1:100],] , method.dist = "correlation", method.hclust = "complete")

{
  mypar()
  plot(pvc,las=2,hang = -0.5)
  pvrect(pvc, alpha = 0.9)
  
  points(1:ncol(data) ,rep(0,ncol(data)), pch= 16, cex=2, col=metadata$Group[pvc$hclust$order])
}

#mclust算法
library(mclust)
fit <- Mclust(data = data[top_var[1:100],])
plot(fit) # plot results

library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(edgeR) # rna-seq

mr$Group <- factor(mr$Group)
d <- DESeqDataSetFromMatrix(countData=cf,colData=mr,design=~Group)

d <- DESeq2::estimateSizeFactors(d,type="ratio")

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gmean <- apply(cf,1,gm_mean)
head(gmean)

ratio <- cf/gmean
head(ratio)[,1:5]
#每个样本（所有基因）的中位数比率被视为该样本的大小因子。
sf <- apply(ratio,2,median)
sf
#通过与DESeq2生成的大小因子进行比较，我们可以验证这些值是正确的。
sizeFactors(d)

plot(sizeFactors(d),colSums(cf),xlab="Size factors",ylab="Total counts")

# custom
head(t(t(cf)/sf))[,1:5]
# deseq2
head(counts(d,normalized=TRUE))[,1:5]

#Gene dispersion
#We can create a mean counts vs variance plot for all genes in our data set.
dm <- apply(cf,1,mean)
dv <- apply(cf,1,var)

ggplot(data.frame(mean=log10(dm+1),var=log10(dv+1)),
       aes(mean,var))+
  geom_point(alpha=0.2)+
  geom_smooth(method="lm")+
  labs(x=expression('Log'[10]~'Mean counts'),y=expression('Log'[10]~'Variance'))+
  theme_bw()
#一种选择是变异系数（CV）。让我们计算每个基因的CV并绘制CV与平均值的关系。

cva <- function(x) sd(x)/mean(x)
dc <- apply(cf,1,cva)

ggplot(data.frame(mean=log10(dm+1),var=dc),
       aes(mean,var))+
  geom_point(alpha=0.2)+
  geom_smooth()+
  labs(x=expression('Log'[10]~'Mean counts'),y="Coefficient of variation")+
  theme_bw()
#如果我们计算CV并求出样本组（时间）内重复样本的平均值，这一点就更加明显

dx1 <- data.frame(Cont=apply(cf[,1:9],1,cva),EPO=apply(cf[,10:18],1,cva))
dx1$gene <- rownames(dx1)
dx1 <- tidyr::gather(dx1,key=sample,value=cv,-gene)
rownames(dx1) <- paste0(dx1$gene,"-",dx1$sample)

dx2 <- data.frame(Cont=apply(cf[,1:9],1,mean),EPO=apply(cf[,10:18],1,mean))
dx2$gene <- rownames(dx2)
dx2 <- tidyr::gather(dx2,key=sample,value=mean,-gene)
rownames(dx2) <- paste0(dx2$gene,"-",dx2$sample)

dx3 <- merge(dx1,dx2,by=0)

ggplot(dx3,aes(x=log10(mean+1),y=cv))+
  geom_point(alpha=0.2)+
  geom_smooth()+
  facet_wrap(~sample.x)+
  labs(x=expression('Log'[10]~'Mean counts'),y="Coefficient of variation")+
  theme_bw()

d <- DESeq2::estimateDispersions(d)
plotDispEsts(d)
#FC = corrected counts group B / corrected counts group A
#logFC = log2(FC)
dg <- nbinomWaldTest(d)
resultsNames(dg)
res <- results(dg,name="Group_NEC_vs_Cont",alpha=0.05)
res$padj[is.na(res$padj)] <- 1
# res <- na.omit(res)
summary(res)
write.csv(res,"./dge_results.csv",row.names=T)