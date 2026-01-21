#load z-score expression
load('./data/PDX_TCGA_CCLE_expression_z_score.RData')

#load indices
load('./data/PDX_TCGA_CCLE_inds.RData')


#Spearman correlation matrix for PDX, TCGA and CCLE (Z-score)

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
rownames(df_)=colnames(PDX_TCGA_CCLE_expression_z_score)

#For simplicity, don't show READ and COADREAD; Combine SKCM and MEL
samples_used=rownames(df_)[df_$histology!='READ' & df_$histology!='COADREAD']
df_=df_[samples_used,]
df_$histology[df_$histology=='SKCM']='MEL'



# Heatmap for all histologies
#Note: The correlation matrix is computationally intensive. Consider running it from HPC with >=64g RAM.
legend_label='z-score'
mat_=PDX_TCGA_CCLE_expression_z_score
ordered_samples=rownames(df_)[order(df_$histology,df_$cohort)]
mat_=mat_[,ordered_samples]
cor_=cor(mat_,method='spearman',use='complete.obs')


ht=ComplexHeatmap::pheatmap(cor_,
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

pdf(file=paste0('Extended_Data_Fig_6A.pdf'),paper = "a4r",width=11)
ComplexHeatmap::draw(ht, merge_legend = TRUE) #Merge duplicated legend
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


pdf(file=paste0('Extended_Data_Fig_6BtoH.pdf'),paper = "a4r",width=11)
for(i in 1:7){
  ComplexHeatmap::draw(heatmap_list[[i]], merge_legend = TRUE) #Merge duplicated legend
}
dev.off()




