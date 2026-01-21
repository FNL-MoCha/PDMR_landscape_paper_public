library(dplyr)

source('./sup_code/self_defined_functions.R')

#load metadata
load('./data/meta_in_pub.RData')

#load gene expression aggregated at model level (Currently they don't include PDC, PDOrg)
load('./data/model_dat_original.RData') 


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


pdf('Extended_Data_Fig_5D.pdf')
boxplot(spearmans_dat_original,ylim=c(0.6,1),cex.axis=2,main='Model vs. originator Spearman correlations')
dev.off()

