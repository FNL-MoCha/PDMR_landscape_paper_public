library(R.matlab) #BiocManager::install("R.matlab")
library(dplyr)
library(estimate)

source('./sup_code/self_defined_functions.R')

#load the PDX data
load('./data/expression_no_organoid_no_PDC.RData')
load('./data/PDX_TCGA_CCLE_inds.RData')
COAD_PDX_expression=expression_no_organoid_no_PDC[,PDX_TCGA_CCLE_inds$PDX$COAD_sample_ind]
PAAD_PDX_expression=expression_no_organoid_no_PDC[,PDX_TCGA_CCLE_inds$PDX$PAAD_sample_ind]
HNSC_PDX_expression=expression_no_organoid_no_PDC[,PDX_TCGA_CCLE_inds$PDX$HNSC_sample_ind]

#For final K selection, correlate all major components to PDX (Considering all originator-PDX pairs)



#COAD
################################################
################################################

#load the deconvoluted matrices
res=readMat("./data/decoder/DECODER_COAD/K4_res.mat")
genesigs=res$genesigs          # genes × K matrix
samplesigs=res$samplesigs          # K × samples

res_K2=readMat("./data/decoder/DECODER_COAD/K2_res.mat")
genesigs_K2=res_K2$genesigs          # genes × K matrix
samplesigs_K2=res_K2$samplesigs          # K × samples



#load the original DECODER input matrix to get gene names
COAD_originator_expression=read.table('./sup_code/COAD_originator_expression.tsv',sep='\t',header = T,check.names=F,row.names=1)
rownames(genesigs)=rownames(genesigs_K2)=rownames(COAD_originator_expression)


#reconstruct expression matrices of optimal K
component_expression_matrix_list=list()
for(i in 1:ncol(genesigs)){
  component_expression_matrix_list[[i]]=genesigs[,i,drop=F] %*% samplesigs[i,,drop=F]
  colnames(component_expression_matrix_list[[i]])=colnames(COAD_originator_expression)
  rownames(component_expression_matrix_list[[i]])=rownames(COAD_originator_expression)
}
names(component_expression_matrix_list)=colnames(component_expression_matrix_list)

major_component_expression=component_expression_matrix_list[[1]]+component_expression_matrix_list[[2]]+component_expression_matrix_list[[3]]
rest_component=component_expression_matrix_list[[4]]


filler=rep(NA,ncol(COAD_originator_expression))
mean_cor_list_K2=list(mean_cors_originator_PDX=filler,
                      mean_cors_DeconOriginatorComponent1_PDX=filler,
                      mean_cors_DeconOriginatorComponent2_PDX=filler,
                      major_component_expression=filler,
                      rest_component=filler)
i=1
for(originator_expression_mat in list(log2(COAD_originator_expression+1),
                                      K2_component_expression_matrix_list[[1]],
                                      K2_component_expression_matrix_list[[2]],
                                      major_component_expression,
                                      rest_component
)){
  mean_cors_=c()
  for(originator in colnames(originator_expression_mat)){
    model_=sub('~ORIGINATOR.*','',originator)
    PDX_=colnames(COAD_PDX_expression)[grepl(model_,colnames(COAD_PDX_expression))]
    
    mean_cors_=c(mean_cors_,
                 cor(originator_expression_mat[,originator,drop=F],COAD_PDX_expression[,PDX_] ,method='spearman') %>% mean)
  }
  mean_cor_list_K2[[i]]=mean_cors_
  i=i+1
}


#Calculate tumor fraction for each major component and aggregate into 2: high tumor fraction vs. low tumor fraction
pdf(file='./fig_3D_1.pdf',paper="letter")

op <- par(no.readonly = TRUE) # save everything to restore later

par(mar = c(10, 4, 4, 2) + 0.1,  # bottom margin large enough for long labels
    pin = c(2.5, 5), # Makes the boxplot thinner. Will be reeset everytime a new device is called
    las = 2)                     # vertical tick labels


boxK1=my_ESTIMATE(component_expression_matrix_list[[1]])[,'TumorPurity']
boxK2=my_ESTIMATE(component_expression_matrix_list[[2]])[,'TumorPurity']
boxK3=my_ESTIMATE(component_expression_matrix_list[[3]])[,'TumorPurity']
boxK4=my_ESTIMATE(component_expression_matrix_list[[4]])[,'TumorPurity']
# boxplot(list(boxK1,boxK2,boxK3,boxK4),
#         names=c('Component 1','Component 2','Component 3','Component 4'),
#         main='ESTIMATE tumor fraction')
#=> Choose component 1,2,4

high_TF_component=component_expression_matrix_list[[1]]+component_expression_matrix_list[[2]]+component_expression_matrix_list[[4]]
low_TF_component=component_expression_matrix_list[[3]]


filler=rep(NA,ncol(COAD_originator_expression))
mean_TF_components_cor_list=list(mean_cors_originator_PDX=filler,
                                 mean_cors_high_TF_component_PDX=filler,
                                 mean_cors_low_TF_component_PDX=filler)
i=1
for(originator_expression_mat in list(log2(COAD_originator_expression+1),
                                      high_TF_component,
                                      low_TF_component)){
  mean_cors_=c()
  for(originator in colnames(originator_expression_mat)){
    model_=sub('~ORIGINATOR.*','',originator)
    PDX_=colnames(COAD_PDX_expression)[grepl(model_,colnames(COAD_PDX_expression))]
    
    mean_cors_=c(mean_cors_,
                 cor(originator_expression_mat[,originator,drop=F],COAD_PDX_expression[,PDX_] ,method='spearman') %>% mean)
  }
  mean_TF_components_cor_list[[i]]=mean_cors_
  i=i+1
}



boxplot(mean_TF_components_cor_list,names=c('Paired originator vs. PDX',
                                            'Paired high TF component vs. PDX',
                                            'Paired low TF component vs. PDX'),
        main='Mean Spearman correlation distribution - COAD',
        ylim=c(0.6,1),
        border=c('#4323ad','#bf2453','#1f6682'),
        las=2
)
par(op) # restore previous settings
dev.off()


#statistical test
wilcox.test(mean_TF_components_cor_list[[1]], mean_TF_components_cor_list[[2]]) #p-value = 0.05449
################################################
################################################


#PAAD
################################################
################################################
#load the deconvoluted matrices
res=readMat("./data/decoder/DECODER_PAAD/K6_res.mat")
genesigs=res$genesigs          # genes × K matrix
samplesigs=res$samplesigs          # K × samples

res_K2=readMat("./data/decoder/DECODER_PAAD/K2_res.mat")
genesigs_K2=res_K2$genesigs          # genes × K matrix
samplesigs_K2=res_K2$samplesigs          # K × samples


#load the original DECODER input matrix to get gene names
PAAD_originator_expression=read.table('./sup_code/PAAD_originator_expression.tsv',sep='\t',header = T,check.names=F,row.names=1)
rownames(genesigs)=rownames(genesigs_K2)=rownames(PAAD_originator_expression)


#reconstruct expression matrices of optimal K
component_expression_matrix_list=list()
for(i in 1:ncol(genesigs)){
  component_expression_matrix_list[[i]]=genesigs[,i,drop=F] %*% samplesigs[i,,drop=F]
  colnames(component_expression_matrix_list[[i]])=colnames(PAAD_originator_expression)
  rownames(component_expression_matrix_list[[i]])=rownames(PAAD_originator_expression)
}
names(component_expression_matrix_list)=colnames(component_expression_matrix_list)

major_component_expression=component_expression_matrix_list[[1]]+component_expression_matrix_list[[2]]+component_expression_matrix_list[[3]]+component_expression_matrix_list[[4]]
rest_component=component_expression_matrix_list[[5]]+component_expression_matrix_list[[6]]



filler=rep(NA,ncol(PAAD_originator_expression))
mean_cor_list_K2=list(mean_cors_originator_PDX=filler,
                      mean_cors_DeconOriginatorComponent1_PDX=filler,
                      mean_cors_DeconOriginatorComponent2_PDX=filler,
                      major_component_expression=filler,
                      rest_component=filler)
i=1
for(originator_expression_mat in list(log2(PAAD_originator_expression+1),
                                      K2_component_expression_matrix_list[[1]],
                                      K2_component_expression_matrix_list[[2]],
                                      major_component_expression,
                                      rest_component
)){
  mean_cors_=c()
  for(originator in colnames(originator_expression_mat)){
    model_=sub('~ORIGINATOR.*','',originator)
    PDX_=colnames(PAAD_PDX_expression)[grepl(model_,colnames(PAAD_PDX_expression))]
    
    mean_cors_=c(mean_cors_,
                 cor(originator_expression_mat[,originator,drop=F],PAAD_PDX_expression[,PDX_] ,method='spearman') %>% mean)
  }
  mean_cor_list_K2[[i]]=mean_cors_
  i=i+1
}


#Calculate tumor fraction for each major component and aggregate into 2: high tumor fraction vs. low tumor fraction
pdf(file='./fig_3D_2.pdf',paper="letter")

op <- par(no.readonly = TRUE) # save everything to restore later

par(mar = c(10, 4, 4, 2) + 0.1,  # bottom margin large enough for long labels
    pin = c(2.5, 5), # Makes the boxplot thinner. Will be reeset everytime a new device is called
    las = 2)                     # vertical tick labels


boxK1=my_ESTIMATE(component_expression_matrix_list[[1]])[,'TumorPurity']
boxK2=my_ESTIMATE(component_expression_matrix_list[[2]])[,'TumorPurity']
boxK3=my_ESTIMATE(component_expression_matrix_list[[3]])[,'TumorPurity']
boxK4=my_ESTIMATE(component_expression_matrix_list[[4]])[,'TumorPurity']
boxK5=my_ESTIMATE(component_expression_matrix_list[[5]])[,'TumorPurity']
boxK6=my_ESTIMATE(component_expression_matrix_list[[6]])[,'TumorPurity']
# boxplot(list(boxK1,boxK2,boxK3,boxK4,boxK5,boxK6),
#         names=c('Component 1','Component 2','Component 3','Component 4','Component 5','Component 6'),
#         main='ESTIMATE tumor fraction')
#=> Choose component 1,2,4,6

high_TF_component=component_expression_matrix_list[[1]]+component_expression_matrix_list[[2]]+component_expression_matrix_list[[4]]+component_expression_matrix_list[[6]]
low_TF_component=component_expression_matrix_list[[3]]+component_expression_matrix_list[[5]]


filler=rep(NA,ncol(PAAD_originator_expression))
mean_TF_components_cor_list=list(mean_cors_originator_PDX=filler,
                                 mean_cors_high_TF_component_PDX=filler,
                                 mean_cors_low_TF_component_PDX=filler)
i=1
for(originator_expression_mat in list(log2(PAAD_originator_expression+1),
                                      high_TF_component,
                                      low_TF_component)){
  mean_cors_=c()
  for(originator in colnames(originator_expression_mat)){
    model_=sub('~ORIGINATOR.*','',originator)
    PDX_=colnames(PAAD_PDX_expression)[grepl(model_,colnames(PAAD_PDX_expression))]
    
    mean_cors_=c(mean_cors_,
                 cor(originator_expression_mat[,originator,drop=F],PAAD_PDX_expression[,PDX_] ,method='spearman') %>% mean)
  }
  mean_TF_components_cor_list[[i]]=mean_cors_
  i=i+1
}



boxplot(mean_TF_components_cor_list,names=c('Paired originator vs. PDX',
                                            'Paired high TF component vs. PDX',
                                            'Paired low TF component vs. PDX'),
        main='Mean Spearman correlation distribution - PAAD',
        ylim=c(0.6,1),
        border=c('#4323ad','#bf2453','#1f6682'),
        las=2
)
par(op) # restore previous settings
dev.off()

#statistical test
wilcox.test(mean_TF_components_cor_list[[1]], mean_TF_components_cor_list[[2]]) #p-value = 0.227

################################################
################################################

#HNSC
################################################
################################################
#load the deconvoluted matrices
res=readMat("./data/decoder/DECODER_HNSC/K6_res.mat")
genesigs=res$genesigs          # genes × K matrix
samplesigs=res$samplesigs          # K × samples

res_K2=readMat("./data/decoder/DECODER_HNSC/K2_res.mat")
genesigs_K2=res_K2$genesigs          # genes × K matrix
samplesigs_K2=res_K2$samplesigs          # K × samples



#load the original DECODER input matrix to get gene names
HNSC_originator_expression=read.table('./sup_code/HNSC_originator_expression.tsv',sep='\t',header = T,check.names=F,row.names=1)
rownames(genesigs)=rownames(genesigs_K2)=rownames(HNSC_originator_expression)


#reconstruct expression matrices of optimal K
component_expression_matrix_list=list()
for(i in 1:ncol(genesigs)){
  component_expression_matrix_list[[i]]=genesigs[,i,drop=F] %*% samplesigs[i,,drop=F]
  colnames(component_expression_matrix_list[[i]])=colnames(HNSC_originator_expression)
  rownames(component_expression_matrix_list[[i]])=rownames(HNSC_originator_expression)
}
names(component_expression_matrix_list)=colnames(component_expression_matrix_list)

major_component_expression=component_expression_matrix_list[[1]]+component_expression_matrix_list[[2]]+component_expression_matrix_list[[3]]+
  component_expression_matrix_list[[4]]+component_expression_matrix_list[[5]]
rest_component=component_expression_matrix_list[[6]]


filler=rep(NA,ncol(HNSC_originator_expression))
mean_cor_list_K2=list(mean_cors_originator_PDX=filler,
                      mean_cors_DeconOriginatorComponent1_PDX=filler,
                      mean_cors_DeconOriginatorComponent2_PDX=filler,
                      major_component_expression=filler,
                      rest_component=filler)
i=1
for(originator_expression_mat in list(log2(HNSC_originator_expression+1),
                                      K2_component_expression_matrix_list[[1]],
                                      K2_component_expression_matrix_list[[2]],
                                      major_component_expression,
                                      rest_component
)){
  mean_cors_=c()
  for(originator in colnames(originator_expression_mat)){
    model_=sub('~ORIGINATOR.*','',originator)
    PDX_=colnames(HNSC_PDX_expression)[grepl(model_,colnames(HNSC_PDX_expression))]
    
    mean_cors_=c(mean_cors_,
                 cor(originator_expression_mat[,originator,drop=F],HNSC_PDX_expression[,PDX_] ,method='spearman') %>% mean)
  }
  mean_cor_list_K2[[i]]=mean_cors_
  i=i+1
}


#Calculate tumor fraction for each major component and aggregate into 2: high tumor fraction vs. low tumor fraction
pdf(file='./fig_3D_3.pdf',paper="letter")

op <- par(no.readonly = TRUE) # save everything to restore later

par(mar = c(10, 4, 4, 2) + 0.1,  # bottom margin large enough for long labels
    pin = c(2.5, 5), # Makes the boxplot thinner. Will be reeset everytime a new device is called
    las = 2)                     # vertical tick labels

boxK1=my_ESTIMATE(component_expression_matrix_list[[1]])[,'TumorPurity']
boxK2=my_ESTIMATE(component_expression_matrix_list[[2]])[,'TumorPurity']
boxK3=my_ESTIMATE(component_expression_matrix_list[[3]])[,'TumorPurity']
boxK4=my_ESTIMATE(component_expression_matrix_list[[4]])[,'TumorPurity']
boxK5=my_ESTIMATE(component_expression_matrix_list[[5]])[,'TumorPurity']
boxK6=my_ESTIMATE(component_expression_matrix_list[[6]])[,'TumorPurity']
# boxplot(list(boxK1,boxK2,boxK3,boxK4,boxK5,boxK6),
#         names=c('Component 1','Component 2','Component 3','Component 4','Component 5','Component 6'),
#         main='ESTIMATE tumor fraction')
#=> Choose component 1,3,4,5,6

high_TF_component=component_expression_matrix_list[[1]]+component_expression_matrix_list[[3]]+component_expression_matrix_list[[4]]+
  component_expression_matrix_list[[5]]+component_expression_matrix_list[[6]]
low_TF_component=component_expression_matrix_list[[2]]


filler=rep(NA,ncol(HNSC_originator_expression))
mean_TF_components_cor_list=list(mean_cors_originator_PDX=filler,
                                 mean_cors_high_TF_component_PDX=filler,
                                 mean_cors_low_TF_component_PDX=filler)
i=1
for(originator_expression_mat in list(log2(HNSC_originator_expression+1),
                                      high_TF_component,
                                      low_TF_component)){
  mean_cors_=c()
  for(originator in colnames(originator_expression_mat)){
    model_=sub('~ORIGINATOR.*','',originator)
    PDX_=colnames(HNSC_PDX_expression)[grepl(model_,colnames(HNSC_PDX_expression))]
    
    mean_cors_=c(mean_cors_,
                 cor(originator_expression_mat[,originator,drop=F],HNSC_PDX_expression[,PDX_] ,method='spearman') %>% mean)
  }
  mean_TF_components_cor_list[[i]]=mean_cors_
  i=i+1
}



boxplot(mean_TF_components_cor_list,names=c('Paired originator vs. PDX',
                                            'Paired high TF component vs. PDX',
                                            'Paired low TF component vs. PDX'),
        main='Mean Spearman correlation distribution - HNSC',
        ylim=c(0.6,1),
        border=c('#4323ad','#bf2453','#1f6682'),
        las=2
)
par(op) # restore previous settings
dev.off()

#statistical test
wilcox.test(mean_TF_components_cor_list[[1]], mean_TF_components_cor_list[[2]]) #p-value = 0.4721
################################################
################################################
