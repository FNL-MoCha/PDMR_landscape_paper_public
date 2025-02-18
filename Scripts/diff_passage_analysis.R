#Note: Each time the samples/models change, I need to rerun: Rscript compute_dim_reduction_objs.R
library(dplyr)
library(Rtsne)
library(ggplot2)

setwd('./') #set this to the working directory

#load gene expression aggregated at model level (Currently they don't include PDC, PDOrg)
load('model_dat_original.RData')
load('model_dat.RData')

###############
method='pca'
#method='tsne'

#data_='original'
data_='combat'
###############

#Load RNA-seq data

##Without ComBat-Seq (original)
if(data_=='original'){
  load('./mutational_landscape_normalizedCount.RData') #loaded as: data.final
} 
##with ComBat-Seq
if(data_=='combat'){
  load('./batch_adjusted_mutational_landscape_normalizedCount.RData')
} 

if(identical(class(data.final),'data.frame')) data.final=as.matrix(data.final)
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])

#Need a df with: PDX sample ID - cancer type - passage
cancer_type_in_publication=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
passage_to_use=c('Originator','P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')

load('meta_in_pub.RData')
dat=data.final[,meta_in_pub$sample] 



if(method=='pca'){
if(data_=='original'){
  load('RNAseq_pca.RData')
  pca=RNAseq_pca
} 
if(data_=='combat'){
  load('batch_adjusted_RNAseq_pca.RData')
  pca=RNAseq_pca
} 
  X=pca$x
  main='PCA'
  xlab='PC1'
  ylab='PC2'
  
  sdev = pca$sdev
  variances = sdev^2
  total_variance = sum(variances)
  explained_variance_percentage = (variances / total_variance*100) %>% signif(3)
}



if(method=='tsne'){
if(data_=='original'){
  load('RNAseq_tsne.RData')
  tsne=RNAseq_tsne
} 
if(data_=='combat'){
  load('batch_adjusted_RNAseq_tsne.RData')
  tsne=RNAseq_tsne
} 
 X=tsne$Y 
 colnames(X)=c('Dimension 1','Dimension 2') 
 rownames(X)=colnames(dat)
 main='t-SNE'
 xlab='Dimension 1'
 ylab='Dimension 2'
}



#Plot overall 2D plot
get_ind1=function(oncotree_code){
  meta_in_pub$sample[meta_in_pub$oncotree==oncotree_code]
}

#total no. of originators, PDX, PDC, PDO
#nrow(meta_in_pub)

ind_COADREAD=get_ind1('COADREAD'); ind_HNSC=get_ind1('HNSC'); ind_NSCLC=get_ind1('NSCLC'); ind_BLCA=get_ind1('BLCA'); ind_PAAD=get_ind1('PAAD'); ind_MEL=get_ind1('MEL'); ind_SARCOMA=get_ind1('SARCOMA')

#xmin, xmax, ymin, ymax
xmin=min(X[,xlab]); xmax=max(X[,xlab]); ymin=min(X[,ylab]); ymax=max(X[,ylab])



#create figure legend for fig 3D

pdf('figures_pdf/fig_3C_legend.pdf')
data <- data.frame(
  Category = factor(c("COADREAD", "HNSC", "NSCLC", "BLCA", "PAAD", "MEL", "SARCOMA"),
                    levels = c("COADREAD", "HNSC", "NSCLC", "BLCA", "PAAD", "MEL", "SARCOMA")),
  Value = c(1, 1, 1, 1, 1, 1, 1)
)


colors <- c("COADREAD" = "#e64b35", "HNSC" = "#0aa087", "NSCLC" = "#3d5488", "BLCA" = "#f49b7f", "PAAD" = "#8e0153", "MEL" = "#8391b4", "SARCOMA" = "#7d6046")

ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 12))
dev.off()




pdf('figures_pdf/fig_3C_1.pdf')
#pdf('figures_pdf/sup_fig_6_1.pdf')

if(method=='pca'){
  xlab_=paste0('PC1 (',explained_variance_percentage[1],'%)')
  ylab_=paste0('PC2 (',explained_variance_percentage[2],'%)')
} 

if(method=='tsne'){
  xlab_=xlab
  ylab_=ylab
}

cex=1
par(mfrow=c(1,1),mar=c(4,4,4,4),mgp=c(1, 0.5, 0))
plot(X[,xlab][ind_COADREAD],X[,ylab][ind_COADREAD], col = "#00000099", bg='#D5554099', pch = 21,cex=cex,cex.lab=2,main='',
     xlab=xlab_,
     ylab=ylab_,
     xlim=c(xmin,xmax),ylim=c(ymin,ymax),
     xaxt='n',yaxt='n') 

cancer_type_cols=c('#539E8799','#3E538399','#E79F8399','#821A5199','#8691B199','#79614A99')
inds=list(ind_HNSC,ind_NSCLC,ind_BLCA,ind_PAAD,ind_MEL,ind_SARCOMA)

for(i in 1:length(cancer_type_cols)){
  points(X[,xlab][inds[[i]]],X[,ylab][inds[[i]]], col = "#00000099", bg=cancer_type_cols[i], pch = 21, cex = cex)
}

#find an COAD example to highlight the drift from originator to all passages
COAD_meta_in_pub=meta_in_pub[meta_in_pub$oncotree=='COADREAD',]
#table(COAD_meta_in_pub$model) %>% sort 
##=> 746757-062-R has the largest no. of samples, but they are scattered
##=> CN0446-F447 has no originator
##=> 983718-287-R has 1 outlier
##=> 947758-054-R has 2 outlier
##=> 944381-210-T has 1 outlier 
##=> 884544-143-R has 1 outlier 
##=> 764851-200-R has 1 outlier 
##=> 762968-020-R has 1 outlier 
##=> 746538-026-R has 1 outlier 
##=> 738494-208-R2 has 1 outlier 
##=> Looks like many of them have 1 outlier
## => Just randomly pick one that has fewer PDX samples: 145919-334-R => This one has no originator
## => Randomly picked 381356-305-R. This one has 1 outlier, but has >=P3 samples

#find an example programmatically
# good_examples=c()
# for(model in unique(COAD_meta_in_pub$model)){
#   model_tab=COAD_meta_in_pub[COAD_meta_in_pub$model==model,]
#   if(('originator' %in% model_tab$Passage) &  #if originator sample is available for that model
#      sum(c('0','1','2') %in% model_tab$Passage)==3 & #if passage 0,1,2 is available for that model
#      sum(model_tab$Passage %in% c("3","4","5","6"))>=1 & #if one of passage 3-6 is available for that model
#      max(X[,xlab][model_tab$SAMPLE_ID])-min(X[,xlab][model_tab$SAMPLE_ID])<=80 & #if the range of PC1 is smaller than the given value
#      max(X[,ylab][model_tab$SAMPLE_ID])-min(X[,ylab][model_tab$SAMPLE_ID])<=80 #if the range of PC2 is smaller than the given value
#      ){
#     good_examples=c(good_examples,
#                     model)
#   }
# }
# good_examples


drift_examples1=rownames(X)[grepl('431354-103-R',rownames(X))]
points(X[,xlab][drift_examples1],X[,ylab][drift_examples1], col = "blue", pch = 1, cex = 2,lwd=2)

drift_examples2=rownames(X)[grepl('551195-344-R',rownames(X))]
points(X[,xlab][drift_examples2],X[,ylab][drift_examples2], col = "green", pch = 1, cex = 2,lwd=2)

dev.off()

#show where 381576-328-R is (The model that has poor within-model correlation ) 
#ind_bad_correlation=rownames(X)[grepl('381576-328-R',rownames(X))]
#points(X[,xlab][ind_bad_correlation],X[,ylab][ind_bad_correlation], col = "red", pch = 1, cex = 2,lwd=2)



#different 2D plots with different passages
get_ind2=function(passage){
  meta_in_pub$sample[meta_in_pub$passage==passage]
}

ind_originator=get_ind2('Originator')
ind_P0=get_ind2('P0'); ind_P1=get_ind2('P1'); ind_P2=get_ind2('P2'); ind_P3=get_ind2('P3'); ind_P4=get_ind2('P4'); ind_P5=get_ind2('P5'); ind_P6=get_ind2('P6')
ind_PDC=get_ind2('PDC'); ind_PDO=get_ind2('PDOrg')

plot_passage_pca=function(ind_Px,main,...){
  cex.axis=1.5
  plot(X[,xlab][intersect(ind_COADREAD,ind_Px)],X[,ylab][intersect(ind_COADREAD,ind_Px)], col = "#00000099", bg='#D5554099', pch = 21, cex.axis = cex.axis,main=main,cex.main=1.5,xlab='',ylab='',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),...)
  
  for(i in 1:length(cancer_type_cols)){
    points(X[,xlab][intersect(inds[[i]],ind_Px)],X[,ylab][intersect(inds[[i]],ind_Px)], col = "#00000099", bg=cancer_type_cols[i], pch = 21, cex = 1)
  }
  
  #points(X[,xlab][intersect(ind_Px,drift_examples1)],X[,ylab][intersect(ind_Px,drift_examples1)], col = "blue", pch = 1, cex = 2,lwd=2) #add drift examples (COAD)
  #points(X[,xlab][intersect(ind_Px,drift_examples2)],X[,ylab][intersect(ind_Px,drift_examples2)], col = "green", pch = 1, cex = 2,lwd=2) #add drift examples (COAD)
  #points(X[,xlab][intersect(ind_Px,ind_bad_correlation)],X[,ylab][intersect(ind_Px,ind_bad_correlation)], col = "red", pch = 1, cex = 2,lwd=2)
}


pdf('figures_pdf/fig_3C_2.pdf')
#pdf('figures_pdf/sup_fig_6_2.pdf')

#plot originators
par(mfrow=c(3,3),mar=c(2,3,2,2))
for(i in 1:2) plot.new() #these 2 plots are place holders
plot_passage_pca(ind_originator,'Originators',xaxt='n',yaxt='n')
text(-20, 115, paste0("n=",length(ind_originator)),col='#00000099',cex=2)
#text(-34, 52, paste0("n=",length(ind_originator)),col='#00000099',cex=1.8)

#plot the PDX, PDC, PDO
main=c('P0','P1','P2','P3-P6','PDC','PDOrg'); i=1
for(data in list(ind_P0,ind_P1,ind_P2,c(ind_P3,ind_P4,ind_P5,ind_P6),ind_PDC,ind_PDO)){
  plot_passage_pca(data,main[i],xaxt='n',yaxt='n')
  text(-20, 115, paste0("n=",length(data)),col='#00000099',cex=2)
  #text(-34, 52, paste0("n=",length(data)),col='#00000099',cex=1.8)
  i=i+1
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


#statistical analysis for all cancer types

#correlation between:

##1. same model (excluding originators) 
spearmans_models=c()
outliers_by_arbitrary_cutoff=list()
for(model in unique(meta_in_pub$model)){
  sample_ind=meta_in_pub$sample[(meta_in_pub$model==model) &
                                     (meta_in_pub$passage %in% c('P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')) 
  ]
  if(length(sample_ind)>=2){
    cors=cor(dat[,sample_ind],method='spearman') %>% as.dist %>% as.numeric
    spearmans_models=c(spearmans_models,
                       cors)
    if(sum(cors<0.75)>=1){
      outliers_by_arbitrary_cutoff[[model]]=sample_ind
    } 
  }
}

##2. same histology (excluding originators)
spearmans_histology=c()
for(histology in unique(meta_in_pub$oncotree)){
  sample_ind=meta_in_pub$sample[(meta_in_pub$oncotree==histology) &
                                     (meta_in_pub$passage %in% c('P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')) 
  ]
  if(length(sample_ind)>=2){
    cors=cor(dat[,sample_ind],method='spearman') %>% as.dist %>% as.numeric
    
    #if(sum(cors>0.95)>=1) cat(sample_ind,'\n')
    
    spearmans_histology=c(spearmans_histology,
                          cors)
  }
}


##3. random pair (excluding originators)
sample_ind=meta_in_pub$sample[meta_in_pub$passage %in% c('P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')]
start=Sys.time()
spearmans_passages=cor(dat[,sample_ind],method='spearman') %>% as.dist %>% as.numeric
end=Sys.time()
end-start #Time difference of 3.05186 mins


#compare correlation distributions
pdf('figures_pdf/fig_3D.pdf')
par(mgp=c(1, 0.5, 0))
boxplot(list(spearmans_models,spearmans_histology,spearmans_passages),ylim=c(0.5,1),xaxt='n',cex=2,cex.axis=2)
dev.off()

wilcox.test(spearmans_models, spearmans_histology) #p-value < 2.2e-16
wilcox.test(spearmans_models, spearmans_passages) #p-value < 2.2e-16
wilcox.test(spearmans_histology, spearmans_passages) #p-value < 2.2e-16



#Using spearman and only consider corresponding originator-passage pairs (This way spearman for originators is simply 1)
start=Sys.time()
originator_dat=log2(dat+1)[,ind_originator]
spearmans_list=list()
for(passage in list('P0','P1','P2',c('P3','P4','P5','P6'),'PDC','PDOrg')){
  cat('doing ',passage)
  
  ind_Px=meta_in_pub$sample[meta_in_pub$passage %in% passage]
  
  p_dat=log2(dat+1)[,ind_Px,drop=F]
  
  
  spearmans=c()
  for(j in 1:ncol(p_dat)){
    model=meta_in_pub$model[meta_in_pub$sample==colnames(p_dat)[j]] 
      
    ind_=which(grepl(model,colnames(originator_dat)))
    
    
    if(!identical(ind_,integer(0))){
      originator_sample=originator_dat[,ind_]
      spearman=cor(p_dat[,j],originator_sample,method='spearman')
      spearmans=c(spearmans,spearman)
    }       
  }
  cat(',n=',length(spearmans),'\n')
  spearmans_list=c(spearmans_list,
                   list(spearmans))
}
end=Sys.time()
end-start #Time difference of 13.1203 secs

pdf('figures_pdf/fig_3D_ComBat.pdf',width = 4*2,height = 3*2)
#pdf('figures_pdf/fig_3D_no_ComBat.pdf',width = 4*2,height = 3*2)
par(mfrow=c(1,1))
boxplot(spearmans_list,ylim=c(0.7,1),xaxt='n',cex.axis=1.5)
dev.off()

print(kruskal.test(spearmans_list[1:4])) #test for H0: P0-P6 all have the same Spearman correlation to originators
sapply(spearmans_list,length) #n for each box




#generate a sample1 - sample 2 - correlation table for Bishu
outlier_cor_table=data.frame()
for(sample_inds_ in outliers_by_arbitrary_cutoff){
  outlier_cors=cor(dat[,sample_inds_],method='spearman')
  outlier_cors[upper.tri(outlier_cors, diag = TRUE)] = NA #remove the redundant half
  cor_table=reshape2::melt(outlier_cors,na.rm=T)
  outlier_cor_table=rbind(outlier_cor_table,
                          cor_table)
}

names(outlier_cor_table)=c('SAMPLE_ID_1','SAMPLE_ID_2','spearman_cor')
outlier_cor_table=outlier_cor_table[outlier_cor_table$spearman_cor<0.8,] #remove correlations>=0.8 from pairings that do not result in outliers
outlier_cor_table=outlier_cor_table[outlier_cor_table$SAMPLE_ID_1!=outlier_cor_table$SAMPLE_ID_2,] #remove self-correlation
outlier_cor_table=left_join(outlier_cor_table,meta_in_pub[,c('sample','oncotree')],by=c('SAMPLE_ID_1'='sample'))
outlier_cor_table=left_join(outlier_cor_table,meta_in_pub[,c('sample','oncotree')],by=c('SAMPLE_ID_2'='sample'))
outlier_cor_table=outlier_cor_table[,c('SAMPLE_ID_1','oncotree.x','SAMPLE_ID_2','oncotree.y','spearman_cor')]
colnames(outlier_cor_table)[2]='OncoTree_1'
colnames(outlier_cor_table)[4]='OncoTree_2'

write.table(outlier_cor_table,file='outlier_mutational_landscape_cor_table.tsv',row.names = F,col.names = T,quote=F,sep='\t')