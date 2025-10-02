source('self_defined_functions.R')

#load metadata
load('../data/meta_in_pub.RData')

#load RNA-seq data
load('../data/mutational_landscape_normalizedCount.RData') #loaded as: data.final
if(identical(class(data.final),'data.frame')) data.final=as.matrix(data.final)
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])
dat=data.final[,meta_in_pub$sample] 





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


##1. correlation within the same model (excluding originators) 
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

##2. correlation within the same histology (excluding originators)
spearmans_histology=c()
for(histology in unique(meta_in_pub$histology_in_figures)){
  sample_ind=meta_in_pub$sample[(meta_in_pub$histology_in_figures==histology) &
                                  (meta_in_pub$passage %in% c('P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')) 
  ]
  if(length(sample_ind)>=2){
    cors=cor(dat[,sample_ind],method='spearman') %>% as.dist %>% as.numeric
    
    #if(sum(cors>0.95)>=1) cat(sample_ind,'\n')
    
    spearmans_histology=c(spearmans_histology,
                          cors)
  }
}


##3. correlation for random pair (excluding originators)
sample_ind=meta_in_pub$sample[meta_in_pub$passage %in% c('P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')]
start=Sys.time()
spearmans_passages=cor(dat[,sample_ind],method='spearman') %>% as.dist %>% as.numeric
end=Sys.time()
end-start #Time difference of 3.05186 mins


#compare correlation distributions
pdf('sup_fig_5K.pdf')
par(mgp=c(1, 0.5, 0))
boxplot(list(spearmans_models,spearmans_histology,spearmans_passages),ylim=c(0,1),xaxt='n',cex=2,cex.axis=2)
dev.off()

wilcox.test(spearmans_models, spearmans_histology) #p-value < 2.2e-16
wilcox.test(spearmans_models, spearmans_passages) #p-value < 2.2e-16
wilcox.test(spearmans_histology, spearmans_passages) #p-value < 2.2e-16

