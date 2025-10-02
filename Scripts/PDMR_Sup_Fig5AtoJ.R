#Note: Paired seq data from TruSeq RNA Exome libraries are stranded (https://www.biostars.org/p/9484953/). This is what we use for our PDX data.
library(tidyverse)
library(SummarizedExperiment)  # Useful for handling the data
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

load('./data/PDX_TCGA_CCLE_expression.RData')
load('./data/PDX_TCGA_CCLE_expression_z_score.RData')

#load indices
load('./data/PDX_TCGA_CCLE_inds.RData')

#load meta data
load('./data/meta_in_pub.RData')


#load pca objects
load('./data/pca_PDX_TCGA_CCLE.RData')
load('./data/pca_PDX_TCGA_CCLE_z_score.RData')


#PCA for 3 cohort - before & after ComBat-Seq
pdf(file=paste0('sup_fig_6AB.pdf'),width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4,4,4),mgp=c(1, 0.5, 0))
types_=c('PCA - original score','PCA - z-score')
j=1
for(pca_ in list(pca_PDX_TCGA_CCLE,pca_PDX_TCGA_CCLE_z_score)){
  X=pca_$x
  
  sdev = pca_$sdev
  variances = sdev^2
  total_variance = sum(variances)
  explained_variance_percentage = (variances / total_variance*100) %>% signif(3)
  
  xlab_=paste0('PC1 (',explained_variance_percentage[1],'%)')
  ylab_=paste0('PC2 (',explained_variance_percentage[2],'%)')
  
  xmin=min(X[,'PC1']); xmax=max(X[,'PC1']); ymin=min(X[,'PC2']); ymax=max(X[,'PC2'])
  
  cex=1
  plot(X[,'PC1'],X[,'PC2'], col = "white", bg='white', pch = 21,cex=cex,cex.lab=1,
       main=types_[j],
       xlab=xlab_,
       ylab=ylab_,
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       xaxt='n',yaxt='n') 
  j=j+1
  
  cols_=c('#87cadeCC','#8a87deCC','#c92c6eCC') #CC at the end means 80% transparency 
  
  
  inds=list(unlist(PDX_TCGA_CCLE_inds$TCGA),unlist(PDX_TCGA_CCLE_inds$PDX),unlist(PDX_TCGA_CCLE_inds$CCLE))
  
  for(i in 1:length(cols_)){
    points(X[,'PC1'][inds[[i]]],X[,'PC2'][inds[[i]]], col = "#00000099", bg=cols_[i], pch = 21, cex = cex)
  }
}
par(mar = c(5, 4, 4, 2)) #set to default
dev.off()




#Spearman correlation matrix for PDX, TCGA and CCLE (Z-score)
#Note: The correlation matrix is computationally intensive. Consider running it from HPC with >=64g RAM.

# Initialize an empty list to collect rows
rows <- list()
# Loop through each cohort: PDX, TCGA, CCLE
for (cohort in names(PDX_TCGA_CCLE_inds)) {
  histology_list <- PDX_TCGA_CCLE_inds[[cohort]]
  
  # Loop through each histology category (e.g., BLCA_sample_ind)
  for (histology in names(histology_list)) {
    samples <- histology_list[[histology]]
    if (length(samples) > 0) {
      rows[[length(rows) + 1]] <- data.frame(
        sample = samples,
        cohort = cohort,
        histology = histology,
        stringsAsFactors = FALSE
      )
    }
  }
}
df_=do.call(rbind, rows)

df_$histology=sub('_ind','',df_$histology)
df_$histology=sub('TCGA|CCLE','',df_$histology)
df_$histology=sub('originator','',df_$histology)
df_$histology=sub('sample','',df_$histology)
df_$histology=sub('_','',df_$histology)
rownames(df_)=df_$sample
df_=df_[,-1]
rownames(df_)=colnames(PDX_TCGA_CCLE_expression)

#For simplicity, don't show READ and COADREAD; Combine SKCM and MEL
samples_used=rownames(df_)[df_$histology!='READ' & df_$histology!='COADREAD']
df_=df_[samples_used,]
df_$histology[df_$histology=='SKCM']='MEL'





legend_label='z-score'
mat_=PDX_TCGA_CCLE_expression_z_score
mat_=mat_[,ordered_samples]
cor_=cor(mat_,method='spearman',use='complete.obs')


pdf(file=paste0('sup_fig_6C.pdf'),paper = "a4r",width=11)
ComplexHeatmap::pheatmap(cor_,
                         annotation_row=df_[ordered_samples,],
                         annotation_col=df_[ordered_samples,],
                         annotation_colors = list(
                           cohort = c("PDX" = "#8a87de", "TCGA" = "#87cade","CCLE"="#c92c6e"),
                           histology = c("COAD"= "#e64b35", "HNSC" = "#0aa087","LUAD"="#3d5488","LUSC"="#749adb",
                                         "BLCA" = "#f49b7f", "PAAD" = "#8e0153", "MEL" = "#8391b4")
                         ),
                         cluster_rows = F,
                         cluster_cols = F,
                         show_colnames = F,show_rownames = F,
                         name=legend_label,
                         breaks = seq(-0.5, 0.5, length.out = 101),  # anchors color scale to [-0.5, 0.5]
                         use_raster=T)
dev.off()



heatmap_list=list()
i=0
histologies=c("BLCA", "COAD", "HNSC", "MEL", "LUAD", "LUSC", "PAAD")
for(histology in histologies){
  i=i+1
  

  df_histology=df_[df_$histology==histology,]
  
  ordered_samples_histology=rownames(df_histology)[order(df_histology$cohort)]
  mat_=PDX_TCGA_CCLE_expression_z_score[,ordered_samples_histology]
  cor_=cor(mat_,method='spearman',use='complete.obs')
  
  heatmap_list[[i]]=ComplexHeatmap::pheatmap(cor_,
                                              annotation_row=df_histology[ordered_samples_histology,],
                                              annotation_col=df_histology[ordered_samples_histology,],
                                              annotation_colors = list(
                                                cohort = c("PDX" = "#8a87de", "TCGA" = "#87cade","CCLE"="#c92c6e"),
                                                histology = c("COAD"= "#e64b35", "HNSC" = "#0aa087","LUAD"="#3d5488","LUSC"="#749adb",
                                                              "BLCA" = "#f49b7f", "PAAD" = "#8e0153", "MEL" = "#8391b4")
                                              ),
                                              cluster_rows = F,
                                              cluster_cols = F,
                                              show_colnames = F,show_rownames = F,
                                              name=histology,
                                              breaks = seq(-0.5, 0.5, length.out = 101),  # anchors color scale to [-0.5, 0.5]
                                              use_raster=T)
}


pdf(file=paste0('sup_fig_6DtoJ.pdf'),paper = "a4r",width=11)
heatmap_list[[1]]
heatmap_list[[2]]
heatmap_list[[3]]
heatmap_list[[4]]
heatmap_list[[5]]
heatmap_list[[6]]
heatmap_list[[7]]
dev.off()

