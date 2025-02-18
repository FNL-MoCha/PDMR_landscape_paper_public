library(dplyr)
library(Rtsne)
#BiocManager::install("sva")
library(sva)

setwd('./') #set this to the working directory
source('./self_defined_functions.R')

cancer_type_in_publication=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
passage_to_use=c('Originator','P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')

#compute necessary dim reduction objects
start=Sys.time()

for(method in c('pca','tsne')){
  for(data_ in c('original','combat')){
    
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
    
    patientID_oncotree=readxl::read_excel('Supp_tables_v27.xlsx',sheet='Supp_Table_1',skip=1)[,c("Patient ID","Tumor histology code for manuscript figures")] %>% as.data.frame
    patientID_oncotree=patientID_oncotree[patientID_oncotree$`Tumor histology code for manuscript figures` %in% cancer_type_in_publication,]
    patientID_oncotree=unique(patientID_oncotree)
    
    patientID_modelID_sampleID_passage=readxl::read_excel('Supp_tables_v27.xlsx',sheet='Supp_Table_2',skip=1)[,c("Patient ID","Model Set ID","Sample ID","PDX Passage (P#), PDC, PDOrg, or Patient Orginator specimen")] %>% as.data.frame
    patientID_modelID_sampleID_passage=patientID_modelID_sampleID_passage[patientID_modelID_sampleID_passage$'PDX Passage (P#), PDC, PDOrg, or Patient Orginator specimen' %in% passage_to_use,]
    
    meta_in_pub=inner_join(patientID_modelID_sampleID_passage,patientID_oncotree,by='Patient ID')
    colnames(meta_in_pub)[3]=c('SAMPLE_ID')  
    
    ##old approach to get sample Id, Oncotree code and passage utilizes these 2 files: 
    ##/Users/wui2/Documents/MoCha/PDX/data_samples.txt
    ##/Users/wui2/Documents/MoCha/PDX/Master.txt
    
    dat=data.final[,meta_in_pub$SAMPLE_ID] 
    
    if(method=='pca'){
    pca=prcomp(t(log2(dat+1)))
    
    RNAseq_pca=pca
    if(data_=='original') save(RNAseq_pca,file='RNAseq_pca.RData')
    if(data_=='combat') save(RNAseq_pca,file='batch_adjusted_RNAseq_pca.RData')
    }

    if(method=='tsne'){
      #do t-SNE equivalence
      set.seed(102)
      tsne = Rtsne(t(scale(log2(dat+1))),check_duplicates=F)
      RNAseq_tsne=tsne
    
      if(data_=='original') save(RNAseq_tsne,file='RNAseq_tsne.RData')
      if(data_=='combat') save(RNAseq_tsne,file='batch_adjusted_RNAseq_tsne.RData')
    }
  }
}


colnames(meta_in_pub)=c('patient','model','sample','passage','oncotree')
save(meta_in_pub,file='meta_in_pub.RData')


#compute model-level dim reduction objs
#########################
passage_to_use=c('originator','0','1','2','3','4','5','6')

#Load RNA-seq data
load('./mutational_landscape_normalizedCount.RData') #loaded as: data.final
data.final_original=data.final
colnames(data.final_original)[grepl('_FFPE',colnames(data.final_original))]=sub('_FFPE$','',colnames(data.final_original)[grepl('_FFPE',colnames(data.final_original))]) ##correct some originator names
colnames(data.final_original)=gsub('~','-',colnames(data.final_original))

load('./batch_adjusted_mutational_landscape_normalizedCount.RData')
colnames(data.final)=gsub('~','-',colnames(data.final))
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))]) ##correct some originator names


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
save(n_samples_for_each_model,file='n_samples_for_each_model.RData')
#########################



#remove batch effect by oncotree codes
#########################
load('./mutational_landscape_normalizedCount.RData') #loaded as: data.final
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])
load('meta_in_pub.RData')
new_data.final=data.final[,meta_in_pub$SAMPLE_ID] 

oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list=list()
OncoTree_code=unique(meta_in_pub$`ONCOTREE CODE used in FIGURES`)
for(code in OncoTree_code){
  cat('doing ',code,'\n')
  
  ind_=meta_in_pub$SAMPLE_ID[meta_in_pub$`ONCOTREE CODE used in FIGURES`==code]
  dat=new_data.final[,ind_]
  
  batch=ifelse(grepl('ORIGINATOR',ind_),'originator','sample')
  
  adjusted = ComBat_seq(as.matrix(dat), batch=batch)
  oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list[[code]]=adjusted
}

save(oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list,file='oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list.RData')
#########################

end=Sys.time()
end-start





