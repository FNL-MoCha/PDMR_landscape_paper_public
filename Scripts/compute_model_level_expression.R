source('./self_defined_functions.R')

#########################
#Load meta data
load('./meta_in_pub.RData')

#Load RNA-seq data
load('./mutational_landscape_normalizedCount.RData') #loaded as: data.final
data.final_original=data.final
colnames(data.final_original)[grepl('_FFPE',colnames(data.final_original))]=sub('_FFPE$','',colnames(data.final_original)[grepl('_FFPE',colnames(data.final_original))]) ##correct some originator names
colnames(data.final_original)=gsub('~','-',colnames(data.final_original))

load('./batch_adjusted_mutational_landscape_normalizedCount.RData')
colnames(data.final)=gsub('~','-',colnames(data.final))
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))]) ##correct some originator names

cancer_type_in_publication=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
passage_to_use=c('originator','0','1','2','3','4','5','6')

dat_list=list()
for(data_ in list(data.final_original,data.final)){
  dat_sample=data_[,meta_in_pub$sample] 
  
  #aggregate the data to the model level
  n_samples_for_each_model=c()
  dat=create_new_tab(rownames(dat_sample),unique(meta_in_pub$model[meta_in_pub$passage!='Originator']))
  for(model in colnames(dat)){
    sample_in_model=meta_in_pub$sample[meta_in_pub$model==model & meta_in_pub$passage!='Originator']
    n_samples_for_each_model=c(n_samples_for_each_model,
                               length(sample_in_model))
    if(length(sample_in_model)==1){
      dat[,model]=dat_sample[,sample_in_model]
    }else{
      dat[,model]=rowMeans(dat_sample[,sample_in_model])
    }
  }
  dat=cbind(dat,dat_sample[,meta_in_pub$passage=='Originator'])
  
  dat_list=c(dat_list,
             list(dat))
}

model_dat_original=dat_list[[1]]
model_dat=dat_list[[2]]

save(model_dat_original,file='model_dat_original.RData')
save(model_dat,file='model_dat.RData')
#########################