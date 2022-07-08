library(dplyr)
library(DESeq2)
library(annotables)
library(tibble)
wt_rawcounts<-DataFrame
column_to_rownames(wt_metadata,"sample")->wt_metadata
  #turns the first column into a rownames
column_to_rownames(wt_rawcounts,"X")
all(rownames(wt_metadata)==colnames(wt_rawcounts))
  #checks if the row and column names are the same
as.integer()
  #if things aren't an integer
dds_wt<-DESeqDataSetFromMatrix(countData = wt_rawcounts,colData=wt_metadata,design=~.condition)
dds_wt<-estimateSizeFactors(dds_wt)
sizeFactors(dds_wt)
normalized_wt_counts<-counts(dds_wt,normalized=TRUE)
#normalizes counts
normalized_wt_counts <- forceMatrixToInteger(normalized_wt_counts)
forceMatrixToInteger <- function(m){
+     apply (m, c (1, 2), function (x) {
  +         (as.integer(x))
  +     })
+ }
> 
  #makes normalized counts into an integer matrix
dds_wt<-DESeqDataSetFromMatrix(countData=normalized_wt_counts, colData=wt_metadata, design=~Condition)
  #data formulation
dds_wt<-DESeq(dds_wt)
  #runs DESeq analysis
plotDispEsts(dds_wt)
  #Plot dispersion estimates
results(dds_wt, alpha=0.05)
  #results of testing
dds_wt<-estimateSizeFactors(dds_wt)
  #size factors

wt_res<-results(dds,alpha=0.10, lfcThreshold = 0.32)
wt_res<-lfcShrink(dds_wt,res=wt_res)
  #DESeq2 contrasts
plotMA(wt_res,ylim=c(-30,30
                     ))
  #DESeq2 LFC shrinkage
mcols(wt_res)
  #DESeq2 results table
head(wt_res, n=10)
  #Identifying differentially expressed genes
summary(wt_res)
wt_res_all<-data.frame(wt_res) 
rownames(wt_res_all) <- gsub("\\..*","",rownames(wt_res_all))
rownames_to_column(wt_res_all,var="ensgene")->wt_res_all
#gets rid of the decimals
left_join(x=wt_res_all,y=grch38[,c("ensgene","symbol","description")])->wt_res_all
  #results
wt_res_sig<-subset(wt_res_all,padj<0.1)
wt_res_sig<-wt_res_sig %>% arrange(padj)
view(wt_res_sig)
#shows significant genes
view(wt_res_all)
sig_norm_counts_wt<-normalized_wt_counts[wt_res_sig$ensgene,]
