library(dplyr)
library(Rtsne)
library(ggplot2)
setwd('./') #set this to the working directory
load('meta_in_pub.RData')
load('./model_dat_original.RData')
load('./model_dat.RData')
load('./n_samples_for_each_model.RData')
cancer_type_in_publication=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
passage_to_use=c('Originator','P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')

#########################
dat=model_dat_original
#dat=model_dat

method='pca'
#method='tsne'
#########################

#distribution for n samples of each model
hist(n_samples_for_each_model,main='Distribution for no. of samples per model',xlab='No. of PDX samples')

if(method=='pca'){
  start=Sys.time()
  pca=prcomp(t(log2(dat+1)))
  end=Sys.time()
  end-start #Time difference of 1.287885 mins
  X=pca$x
  xlab='PC1'
  ylab='PC2'
}

if(method=='tsne'){
  #do t-SNE
  set.seed(102)
  start=Sys.time()
  tsne = Rtsne(t(scale(log2(dat+1))),check_duplicates=F)
  end=Sys.time()
  end-start #Time difference of ~1min
  
  X=tsne$Y 
  colnames(X)=c('Dimension 1','Dimension 2') 
  rownames(X)=colnames(dat)
  xlab='Dimension 1'
  ylab='Dimension 2'
}

#Plot overall dimension reduction
cex=1
get_ind1=function(oncotree_code){
  originators=meta_in_pub$sample[grepl('ORIGINATOR',meta_in_pub$sample) & meta_in_pub$oncotree==oncotree_code]
  models=meta_in_pub$model[meta_in_pub$oncotree==oncotree_code] %>% unique
  
  return(c(originators,models))
}
ind_COADREAD=get_ind1('COADREAD'); ind_HNSC=get_ind1('HNSC'); ind_NSCLC=get_ind1('NSCLC'); ind_BLCA=get_ind1('BLCA'); ind_PAAD=get_ind1('PAAD'); ind_MEL=get_ind1('MEL'); ind_SARCOMA=get_ind1('SARCOMA')

#xmin, xmax, ymin, ymax
xmin=min(X[,xlab]); xmax=max(X[,xlab]); ymin=min(X[,ylab]); ymax=max(X[,ylab])


if(method=='pca'){
  sdev = pca$sdev
  variances = sdev^2
  total_variance = sum(variances)
  explained_variance_percentage = (variances / total_variance*100) %>% signif(3)
  xlab_=paste0('PC1 (',explained_variance_percentage[1],'%)')
  ylab_=paste0('PC2 (',explained_variance_percentage[2],'%)')
}

if(method=='tsne'){
  xlab_=''
  ylab_=''
} 

pdf('figures_pdf/sup_fig_7_pca_original.pdf',width=8*1.2,height=2.2*1.2)
#pdf('figures_pdf/sup_fig_7_pca_combat.pdf',width=8*1.2,height=2.2*1.2)
#pdf('figures_pdf/sup_fig_7_tsne_original.pdf',width=8*1.2,height=2.2*1.2)
#pdf('figures_pdf/sup_fig_7_tsne_combat.pdf',width=8*1.2,height=2.2*1.2)

cex=1
par(mfrow=c(1,3),mar=c(4,4,4,4),mgp=c(1, 0.5, 0))
plot(X[,xlab][ind_COADREAD],X[,ylab][ind_COADREAD], col = "#00000099", bg='#D5554099', pch = 21, cex = cex,main='All',xlab=xlab_,ylab=ylab_,
     xlim=c(xmin,xmax),ylim=c(ymin,ymax),xaxt='n',yaxt='n')

cancer_type_cols=c('#539E8799','#3E538399','#E79F8399','#821A5199','#8691B199','#79614A99')
inds=list(ind_HNSC,ind_NSCLC,ind_BLCA,ind_PAAD,ind_MEL,ind_SARCOMA)
for(i in 1:length(cancer_type_cols)){
  points(X[,xlab][inds[[i]]],X[,ylab][inds[[i]]], col = "#00000099", bg=cancer_type_cols[i], pch = 21, cex = cex)
}
#length(unique((meta_in_pub$model)))

#separate dim reduction plots for originators and non-originators
cex.axis=1.5
first_ind=meta_in_pub$sample[meta_in_pub$passage=='Originator' & (meta_in_pub$oncotree=='COADREAD')]
plot(X[,xlab][first_ind],X[,ylab][first_ind], col = "#00000099", bg='#D5554099', pch = 21, cex.axis = cex.axis,main='Originators',cex.main=1.5,xlab='',ylab='',
     xlim=c(xmin,xmax),ylim=c(ymin,ymax),xaxt='n',yaxt='n')
rest_cancer_type_in_publication=cancer_type_in_publication[-1]
for(i in 1:length(cancer_type_cols)){
  ind=meta_in_pub$sample[meta_in_pub$passage=='Originator' & (meta_in_pub$oncotree==rest_cancer_type_in_publication[i])]
  points(X[,xlab][ind],X[,ylab][ind], col = "#00000099", bg=cancer_type_cols[i], pch = 21, cex = cex)
}
#sum(meta_in_pub$'PDX Passage or PDC or PDOrg or Patient Orginator specimen'=='Originator')


get_ind2=function(oncotree_code){
  meta_in_pub$model[meta_in_pub$oncotree==oncotree_code & !grepl('ORIGINATOR',meta_in_pub$sample)] %>% unique
}
plot(X[,xlab][ind_COADREAD],X[,ylab][ind_COADREAD], col = "#00000099", bg='#D5554099', pch = 21, cex.axis = cex.axis,main='Models',cex.main=1.5,xlab='',ylab='',
     xlim=c(xmin,xmax),ylim=c(ymin,ymax),xaxt='n',yaxt='n')
for(i in 1:length(cancer_type_cols)){
  type=rest_cancer_type_in_publication[i]
  ind=get_ind2(type)
  points(X[,xlab][ind],X[,ylab][ind], col = "#00000099", bg=cancer_type_cols[i], pch = 21, cex = cex)
}


dev.off()


#statistical analysis for all cancer types using spearman
compute_spearmans=function(dat){
  ind_originators=meta_in_pub$sample[meta_in_pub$passage=='Originator']
  ind_models=meta_in_pub$model[!grepl('ORIGINATOR',meta_in_pub$sample)] %>% unique
  
  originator=log2(dat+1)[,ind_originators]
  models=log2(dat+1)[,ind_models]
  
  spearmans=c()
  names_spearmans=c()
  for(model_name in colnames(models)){
    model_expression=models[,model_name]
    originator_name=paste(model_name,'ORIGINATOR',sep='-')
    
    if(originator_name %in% colnames(originator)){
      originator_expression=originator[,originator_name]
      spearmans=c(spearmans,
                  cor(model_expression,originator_expression,method='spearman'))
      names_spearmans=c(names_spearmans,
                        paste(originator_name,model_name,sep=' - '))
    }
  }
  names(spearmans)=names_spearmans
  
  return(spearmans)
}

spearmans_dat_original=compute_spearmans(model_dat_original)
spearmans_dat=compute_spearmans(model_dat)

par(mfrow=c(1,1))
boxplot(spearmans_dat_original,spearmans_dat,ylim=c(0.6,1),cex.axis=2)

wilcox.test(spearmans_dat_original, spearmans_dat)


