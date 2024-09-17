#loading the tables in ####
gene_counts = read.csv("C:\\Users\\2873826G\\Documents\\My Masters Degree UofG\\RNA-seq and NGS\\Assignment\\Part1 Bulk\\gene_count_matrix.csv",row.names = 1)
transcript_counts = read.csv("C:\\Users\\2873826G\\Documents\\My Masters Degree UofG\\RNA-seq and NGS\\Assignment\\Part1 Bulk\\transcript_count_matrix.csv",row.names = 1)
exp_design = read.csv("C:\\Users\\2873826G\\Documents\\My Masters Degree UofG\\RNA-seq and NGS\\Assignment\\Part1 Bulk\\Design-Exp.txt",header = TRUE,row.names = 1)

#loading the needed libraries ####
library(DESeq2)
library(ggrepel)
library(vsn)
library(ashr)
library(limma)

#using plotPCA function
plotPCAWithSampleNames = function(x, targets=targets, intgroup=colnames(targets)[1], ntop=500)
{
  library(RColorBrewer)
  library(genefilter)
  library(lattice)
  
  # pca
  #rv = rowVars(assay(x))
  rv = rowVars(x)
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  #pca = prcomp(t(assay(x)[select,]))
  pca = prcomp(t(x[select,]))
  
  # proportion of variance
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # sample names
  names = colnames(x)
  #names = as.character(x$sample)
  
  # factor of groups
  #fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  fac = factor(apply(as.data.frame(targets[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  # colors
  if( nlevels(fac) >= 10 )
    colors = rainbow(nlevels(fac))
  else if( nlevels(fac) >= 3 )
    colors = brewer.pal(nlevels(fac), "Set1")
  else
    colors = c( "dodgerblue3", "firebrick3" )
  
  # plot
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    aspect = "fill",
    col = colors,
    xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
    ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}

#reordering the columns in the gene and transcript counts tables to match the experiment design table (and make sense) ####
col_order = c(1,5,6,7,8,9,10,11,12,2,3,4)
gene_counts = gene_counts[,col_order]
transcript_counts = transcript_counts[,col_order]

#creating dds objects for both ####
dds_gene = DESeqDataSetFromMatrix(countData=gene_counts,colData=exp_design, design= ~ group)
dds_transcript = DESeqDataSetFromMatrix(countData=transcript_counts,colData=exp_design, design= ~ group)

#making sure that the dds objects were created correctly. ####
is.matrix(counts(dds_gene))
is.matrix(counts(dds_transcript))
#both are TRUE, moving on

#pre-filtering ####
#for the pre-filtering, keep the genes that have a minimum count of 10 reads in at least 3 samples
filtered_rows_gene = (rowSums(counts(dds_gene) >= 10)) >= 3
dds_gene = dds_gene[filtered_rows_gene,]
#for transcripts, we only filter out those that are all zeros
filtered_rows_transcript = (rowSums(counts(dds_transcript))>0)
dds_transcript = dds_transcript[filtered_rows_transcript,]

#removing filtered_rows variables because they are not needed now.
#same with col_order
rm(filtered_rows_gene)
rm(filtered_rows_transcript)
rm(col_order)

#running DE analysis for both ####
dds_gene = DESeq(dds_gene)
dds_transcript = DESeq(dds_transcript)

#making dispersion plots for both objects ####
plotDispEsts(dds_gene)
plotDispEsts(dds_transcript)

#performing rlog-based transformation for both objects, and then making PCA plots. #### 
rld_gene = rlog(dds_gene,blind = TRUE)
rld_transcript = rlog(dds_transcript,blind = TRUE)
# matrix of log2(raw counts)
lgc.raw.gene = log2(counts(dds_gene,normalized=FALSE)+0.000001)
lgc.raw.transcript = log2(counts(dds_transcript,normalized=FALSE)+0.000001)
# matrix of log2(normalized counts)
lgc.norm.gene = log2(counts(dds_gene,normalized=TRUE)+0.000001)
lgc.norm.transcript = log2(counts(dds_transcript,normalized=TRUE)+0.000001)
#creating the PCA plots now
plotPCA(rld_gene,intgroup=c("group")) + geom_label_repel(aes(label= rownames(exp_design))) 
plotPCA(rld_transcript,intgroup=c("group")) + geom_label_repel(aes(label= rownames(exp_design))) 

#at this point in the script, we work with the gene data only.####

#making SD vs mean plots ####
meanSdPlot(lgc.norm.gene, ranks = TRUE)
meanSdPlot((as.matrix(assay(rld_gene))), ranks = TRUE)

#extracting DE results for BvA and CvA with unshrunken LFC ####
#first for standard NULL hypothesis (LFC=0)
BvA.unshrunken.lfc0 = results(dds_gene,contrast = c("group","B","A"))
CvA.unshrunken.lfc0 = results(dds_gene,contrast = c("group","C","A"))
#for NULL hypothesis of LFC<1
BvA.unshrunken.lfc1 = results(dds_gene,contrast = c("group","B","A"),lfcThreshold = 1)
CvA.unshrunken.lfc1 = results(dds_gene,contrast = c("group","C","A"),lfcThreshold = 1)

#finding the number of significant genes in these DE comparisons
check_sig_genes = subset(BvA.unshrunken.lfc0,BvA.unshrunken.lfc0$padj<0.05)
summary(check_sig_genes) #547
check_sig_genes = subset(BvA.unshrunken.lfc1,BvA.unshrunken.lfc1$padj<0.05)
summary(check_sig_genes) #46
check_sig_genes = subset(CvA.unshrunken.lfc0,CvA.unshrunken.lfc0$padj<0.05)
summary(check_sig_genes) #35
check_sig_genes = subset(CvA.unshrunken.lfc1,CvA.unshrunken.lfc1$padj<0.05)
summary(check_sig_genes) #5

#creating MA plots ####
DESeq2::plotMA(BvA.unshrunken.lfc0,alpha=0.05)
DESeq2::plotMA(CvA.unshrunken.lfc0,alpha=0.05)
DESeq2::plotMA(BvA.unshrunken.lfc1,alpha=0.05)
DESeq2::plotMA(CvA.unshrunken.lfc1,alpha=0.05)

#extracting results and making plots same as above, but this time with shrunken data ####
BvA.shrunken.lfc0 = lfcShrink(dds_gene,contrast = c("group","B","A"),type = 'ashr')
CvA.shrunken.lfc0 = lfcShrink(dds_gene,contrast = c("group","C","A"),type = 'ashr')
BvA.shrunken.lfc1 = lfcShrink(dds_gene,contrast = c("group","B","A"),type = 'ashr',lfcThreshold = 1)
CvA.shrunken.lfc1 = lfcShrink(dds_gene,contrast = c("group","C","A"),type = 'ashr',lfcThreshold = 1)

#writing the Significant.Genes.BvsA.csv and Significant.Genes.CvsA.csv files ####
file_to_write = subset(BvA.shrunken.lfc1,BvA.shrunken.lfc1$padj<0.05 & abs(BvA.shrunken.lfc1$log2FoldChange > 1)) #creating a subset by p values and LFC
file_to_write = file_to_write[order(file_to_write$log2FoldChange),] #ordering by the log fold change
write.csv(file_to_write,file = "Significant.Genes.BvsA.csv") #writing the file
#repeating for CvA
file_to_write = subset(CvA.shrunken.lfc1,CvA.shrunken.lfc1$padj<0.05 & abs(CvA.shrunken.lfc1$log2FoldChange) > 1) 
file_to_write = file_to_write[order(file_to_write$log2FoldChange),] 
write.csv(file_to_write,file = "Significant.Genes.CvsA.csv") 
#removing unneeded variable after it has been written
rm(file_to_write)

#creating MA plots ####
DESeq2::plotMA(BvA.shrunken.lfc0,alpha=0.05)
DESeq2::plotMA(CvA.shrunken.lfc0,alpha=0.05)
DESeq2::plotMA(BvA.shrunken.lfc1,alpha=0.05)
DESeq2::plotMA(CvA.shrunken.lfc1,alpha=0.05)

#removing batch-effect using limma ####
exp_design_withbatch = read.csv('Design.Group.Batch.csv',row.names = 1)

dds.gb <- DESeqDataSetFromMatrix(countData=gene_counts,colData=exp_design_withbatch, design= ~ batch + group)
filtered_rows_gene = (rowSums(counts(dds.gb) >= 10)) >= 3
dds.gb = dds.gb[filtered_rows_gene,]
rm(filtered_rows_gene)
dds.gb <- DESeq(dds.gb)

mydesign <- model.matrix(design(dds_gene),colData(dds_gene))
b.corrected <- limma::removeBatchEffect(assay(rld_gene), batch=colData(dds.gb)$batch, design=mydesign)

#creating a PCA plot ####
plotPCAWithSampleNames(b.corrected,targets=exp_design_withbatch, intgroup='group')

#writing results to a file called BatchCorrected.Rlog.csv
write.csv(b.corrected,file = "BatchCorrected.Rlog.csv") 

#done