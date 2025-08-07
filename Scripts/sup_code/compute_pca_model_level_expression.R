#compute necessary dim reduction objects

library(dplyr)
library(Rtsne)
#BiocManager::install("sva")
library(sva)

source('self_defined_functions.R')


#compute meta data
cancer_type_in_publication=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
passage_to_use=c('Originator','P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')

patientID_oncotree=readxl::read_excel('../data/Supp_tables_v31.xlsx',sheet='Supp_Table_1',skip=1)[,c("Patient ID","ONCOTREE CODE","Tumor histology code for manuscript figures")] %>% as.data.frame
patientID_oncotree=unique(patientID_oncotree)
patientID_modelID_sampleID_passage=readxl::read_excel('../data/Supp_tables_v31.xlsx',sheet='Supp_Table_2',skip=1)[,c("Patient ID","Model Set ID","Sample ID","PDX Passage (P#), PDC, PDOrg, or Patient Orginator specimen")] %>% as.data.frame
patientID_modelID_sampleID_passage=patientID_modelID_sampleID_passage[patientID_modelID_sampleID_passage$"PDX Passage (P#), PDC, PDOrg, or Patient Orginator specimen" %in% passage_to_use,]
meta_in_pub=inner_join(patientID_modelID_sampleID_passage,patientID_oncotree,by='Patient ID')
colnames(meta_in_pub)=c('patient','model','sample','passage','oncotree','histology_in_figures')
save(meta_in_pub,file='../data/meta_in_pub.RData')



#Load RNA-seq data
load('../data/mutational_landscape_normalizedCount.RData') #loaded as: data.final

if(identical(class(data.final),'data.frame')) data.final=as.matrix(data.final)
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])
dat=data.final[,meta_in_pub$sample] 
save(meta_in_pub,file='../data/meta_in_pub.RData')

##PCA
pca=prcomp(t(log2(dat+1)))
RNAseq_pca=pca
save(RNAseq_pca,file='../data/RNAseq_pca.RData')




#compute model-level expression
passage_to_use=c('originator','0','1','2','3','4','5','6')

dat_sample=data.final[,meta_in_pub$sample] 

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
model_dat_original=dat

#model_dat_original is the result
save(model_dat_original,file='../data/model_dat_original.RData')
