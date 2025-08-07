#Need to run this on HPC

#Load data before and after ComBat-Seq correction
load('../data/PDX_TCGA_CCLE_expression.RData')
load('../data/PDX_TCGA_CCLE_expression_z_score.RData')


pca_PDX_TCGA_CCLE=prcomp(t(log2(PDX_TCGA_CCLE_expression+1)))

PDX_TCGA_CCLE_expression_z_score=impute_expression(PDX_TCGA_CCLE_expression_z_score)
pca_PDX_TCGA_CCLE_z_score=prcomp(t(PDX_TCGA_CCLE_expression_z_score))

save(pca_PDX_TCGA_CCLE,file='../data/pca_PDX_TCGA_CCLE.RData')
save(pca_PDX_TCGA_CCLE_z_score,file='../data/pca_PDX_TCGA_CCLE_z_score.RData')
