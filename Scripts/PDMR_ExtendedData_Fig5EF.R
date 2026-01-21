library(estimate)
source('./sup_code/self_defined_functions.R')
load('./data/PDX_TCGA_CCLE_expression.RData')
load('./data/PDX_TCGA_CCLE_expression_z_score.RData')

#load indices
load('./data/PDX_TCGA_CCLE_inds.RData')

#load meta data
load('./data/meta_in_pub.RData')

#ESTIMATE for different passages
gene_expression=PDX_TCGA_CCLE_expression
pdf(file=paste0('Extended_Data_Fig_5E.pdf'),paper = "a4r",width=11)
par(mfrow=c(2,3))
par(mar = c(2,4,2,2))
for(histology in c('BLCA','COAD','READ','HNSC','MEL','LUAD','LUSC','PAAD')){ #COADREAD not used from PDX data
  
  histology_TCGA_CCLE=histology
  if(histology=='MEL') histology_TCGA_CCLE='SKCM'    
  
  originator_inds_=PDX_TCGA_CCLE_inds$PDX[!grepl('COADREAD',names(PDX_TCGA_CCLE_inds$PDX)) & grepl(histology,names(PDX_TCGA_CCLE_inds$PDX)) & grepl('originator',names(PDX_TCGA_CCLE_inds$PDX))][[1]]
  
  PDX_inds_=PDX_TCGA_CCLE_inds$PDX[grepl(histology,names(PDX_TCGA_CCLE_inds$PDX)) & grepl('sample',names(PDX_TCGA_CCLE_inds$PDX))][[1]]
  P0_inds=PDX_inds_[gsub('~','-',PDX_inds_) %in% meta_in_pub$sample[meta_in_pub$passage=='P0']]
  P1_inds=PDX_inds_[gsub('~','-',PDX_inds_) %in% meta_in_pub$sample[meta_in_pub$passage=='P1']]
  P2_onwards_inds=PDX_inds_[gsub('~','-',PDX_inds_) %in% meta_in_pub$sample[meta_in_pub$passage %in% c('P2','P3','P4','P5','P6')]]  
  
  Originator=my_ESTIMATE(gene_expression[,originator_inds_])[,'TumorPurity']
  PDX=my_ESTIMATE(gene_expression[,PDX_inds_])[,'TumorPurity']
  P0=my_ESTIMATE(gene_expression[,P0_inds])[,'TumorPurity']
  P1=my_ESTIMATE(gene_expression[,P1_inds])[,'TumorPurity']
  P2_onwards=my_ESTIMATE(gene_expression[,P2_onwards_inds])[,'TumorPurity']
  
  boxplot(Originator,P0,P1,P2_onwards,
          names=c('Originator',
                  'P0',
                  'P1',
                  '>= P2'
          ),
          ylim=c(0,1),
          main=histology,
          ylab='Tumor fraction',
          cex.axis=1)
}
dev.off()
par(mfrow=c(1,1)) #set to default
par(mar = c(5, 4, 4, 2)) #set to default



#Scatter plot for tumor fraction vs. correlation for originator-PDX
gene_expression=log2(PDX_TCGA_CCLE_expression+1)
pdf(file=paste0('Extended_Data_Fig_5F.pdf'))
par(mfrow=c(3,3),mar = c(2,2,2,1))
for(histology in c('BLCA','COAD','READ','HNSC','MEL','LUAD','LUSC','PAAD')){ #COADREAD not used from PDX data
  
  originator_ind_=names(PDX_TCGA_CCLE_inds$PDX)[!grepl('COADREAD',names(PDX_TCGA_CCLE_inds$PDX)) & grepl(histology,names(PDX_TCGA_CCLE_inds$PDX)) & grepl('originator',names(PDX_TCGA_CCLE_inds$PDX))]
  originator_tumor_fraction=my_ESTIMATE(gene_expression[,PDX_TCGA_CCLE_inds$PDX[[originator_ind_]]])[,'TumorPurity']
  
  mean_paired_originator_PDX_spearman=c()
  for(originator_ in names(originator_tumor_fraction)){
    originator_expression=gene_expression[,originator_,drop=F]
    ind_=names(PDX_TCGA_CCLE_inds$PDX)[!grepl('COADREAD',names(PDX_TCGA_CCLE_inds$PDX)) & grepl(histology,names(PDX_TCGA_CCLE_inds$PDX)) & grepl('sample',names(PDX_TCGA_CCLE_inds$PDX))]
    PDX_expression=gene_expression[,PDX_TCGA_CCLE_inds$PDX[[ind_]][grepl(sub('~ORIGINATOR','',originator_),PDX_TCGA_CCLE_inds$PDX[[ind_]])],drop=F]
    mean_paired_originator_PDX_spearman=c(mean_paired_originator_PDX_spearman,
                                          mean(cor(originator_expression,PDX_expression,method='spearman')))
  }
  
  plot(originator_tumor_fraction,
       mean_paired_originator_PDX_spearman,
       xlab='',ylab='',
       xlim=c(0.3,1),ylim=c(0.5,1),
       main=histology)
  
  lm_=lm(mean_paired_originator_PDX_spearman~originator_tumor_fraction)
  abline(a=lm_$coefficients[1],b=lm_$coefficients[2],lty='dotted',col = "blue")
  
  ct=cor.test(originator_tumor_fraction,
              mean_paired_originator_PDX_spearman,
              method = "pearson",
              use    = "complete.obs")
  
  pearson_cor=ct$estimate   # correlation (r)
  pearson_cor_pvalue=ct$p.value    # two-sided p-value
  
  text(x=0.6,y=0.6,labels=paste0('Pearson for (x,y)=',signif(pearson_cor,2)))
  text(x=0.6,y=0.55,labels=paste0('p-value=',signif(pearson_cor_pvalue,2)),font  = 3 ) #font  = 3: italic
}

#add legend for diagonal lines
par(mar = c(0, 0, 0, 0))          
plot.new()                        # empty canvas in the 4th cell
legend("center",
       legend = c("x: Originator tumor fraction -\n ESTIMATE","y: Mean Spearman correlation for\n paired originator-PDX\n gene expression","Regression line"), 
       lty=c(NA,NA,'dotted'),
       pt.cex = 1,
       col    = c('black','black',"blue"))
par(mar = c(5, 4, 4, 2)) #set to default

dev.off()
par(mfrow=c(1,1))
