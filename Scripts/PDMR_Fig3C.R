library(dplyr)
library(Rtsne)
library(ggplot2)

source('self_defined_functions.R')



#Load RNA-seq data
load('./data/mutational_landscape_normalizedCount.RData') #loaded as: data.final
if(identical(class(data.final),'data.frame')) data.final=as.matrix(data.final)
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])

#Need a df with: PDX sample ID - cancer type - passage
cancer_type_in_publication=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
passage_to_use=c('Originator','P0','P1','P2','P3','P4','P5','P6','PDC','PDOrg')

#load metadata
load('./data/meta_in_pub.RData')

dat=data.final[,meta_in_pub$sample] 



load('./data/RNAseq_pca.RData')
pca=RNAseq_pca

X=pca$x
main='PCA'
xlab='PC1'
ylab='PC2'

sdev = pca$sdev
variances = sdev^2
total_variance = sum(variances)
explained_variance_percentage = (variances / total_variance*100) %>% signif(3)


#Plot overall 2D plot
get_ind1=function(oncotree_code){
  meta_in_pub$sample[meta_in_pub$histology_in_figures==oncotree_code]
}

#total no. of originators, PDX, PDC, PDO
#nrow(meta_in_pub)

ind_COADREAD=get_ind1('COADREAD'); ind_HNSC=get_ind1('HNSC'); ind_NSCLC=get_ind1('NSCLC'); ind_BLCA=get_ind1('BLCA'); ind_PAAD=get_ind1('PAAD'); ind_MEL=get_ind1('MEL'); ind_SARCOMA=get_ind1('SARCOMA')

#xmin, xmax, ymin, ymax
xmin=min(X[,xlab]); xmax=max(X[,xlab]); ymin=min(X[,ylab]); ymax=max(X[,ylab])



#create figure legend for fig 3D

pdf('fig_3C_legend.pdf')
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




pdf('fig_3C_1.pdf')
xlab_=paste0('PC1 (',explained_variance_percentage[1],'%)')
ylab_=paste0('PC2 (',explained_variance_percentage[2],'%)')

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
COAD_meta_in_pub=meta_in_pub[meta_in_pub$histology_in_figures=='COADREAD',]

drift_examples1=rownames(X)[grepl('431354-103-R',rownames(X))]
points(X[,xlab][drift_examples1],X[,ylab][drift_examples1], col = "blue", pch = 1, cex = 2,lwd=2)

drift_examples2=rownames(X)[grepl('551195-344-R',rownames(X))]
points(X[,xlab][drift_examples2],X[,ylab][drift_examples2], col = "green", pch = 1, cex = 2,lwd=2)

dev.off()



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
}


pdf('fig_3C_2.pdf')


#plot originators
par(mfrow=c(3,3),mar=c(2,3,2,2))
for(i in 1:2) plot.new() #these 2 plots are place holders
plot_passage_pca(ind_originator,'Originators',xaxt='n',yaxt='n')
text(130, 220, paste0("n=",length(ind_originator)),col='#00000099',cex=1.6)
#text(-34, 52, paste0("n=",length(ind_originator)),col='#00000099',cex=1.6)

#plot the PDX, PDC, PDO
main=c('P0','P1','P2','P3-P6','PDC','PDOrg'); i=1
for(data in list(ind_P0,ind_P1,ind_P2,c(ind_P3,ind_P4,ind_P5,ind_P6),ind_PDC,ind_PDO)){
  plot_passage_pca(data,main[i],xaxt='n',yaxt='n')
  text(130, 220, paste0("n=",length(data)),col='#00000099',cex=1.6)
  #text(-34, 52, paste0("n=",length(data)),col='#00000099',cex=1.6)
  i=i+1
}

dev.off()
par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))



