#need to run this script first: compute_dim_reduction_objs.R

library(dplyr)
library(Rtsne)
setwd('./') #set this to the working directory

method='pca'
#method='tsne'

if(method=='pca'){
  xlab='PC1'
  ylab='PC2'
}

if(method=='tsne'){
  xlab=''
  ylab=''
}


load('./oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list.RData')
load('./meta_in_pub.RData')
load('./mutational_landscape_normalizedCount.RData') #loaded as: data.final
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])

new_data.final=data.final[,meta_in_pub$sample]

cancer_type=c('COADREAD','HNSC','NSCLC','BLCA','PAAD','MEL','SARCOMA')
cancer_type_cols=c('#D5554099','#539E8799','#3E538399','#E79F8399','#821A5199','#8691B199','#79614A99')



start=Sys.time()
dim_reduction_list=list()
for(i in seq(oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list)){
  dat=oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list[[i]]
  if(method=='pca') dim_reduction_list[[i]]=prcomp(t(log2(dat+1)))
  if(method=='tsne'){
    set.seed(102)
    dim_reduction_list[[i]] = Rtsne(t(scale(log2(dat+1))),check_duplicates=F)
  } 
  } 

names(dim_reduction_list)=names(oncotreeCode_batch_adjusted_mutational_landscape_normalizedCount_list)

end=Sys.time()
end-start #Time difference of 1.380356 mins


start=Sys.time()
no_combat_dim_reduction_list=list()
for(i in seq(cancer_type)){
  dat=new_data.final[,meta_in_pub$oncotree==cancer_type[[i]]]
  if(method=='pca') no_combat_dim_reduction_list[[i]]=prcomp(t(log2(dat+1)))
  if(method=='tsne'){
    set.seed(102)
    no_combat_dim_reduction_list[[i]] = Rtsne(t(scale(log2(dat+1))),check_duplicates=F)
  } 
}
names(no_combat_dim_reduction_list)=cancer_type

end=Sys.time()
end-start #Time difference of 1.456895 mins


plot_all_dim_reduction=function(dim_reduction_list,...){
  cex=1
  cex.lab = 1.5
  par(mfrow=c(3,3),mar=c(2,3,2,2),mgp=c(1, 0.5, 0))
  
  for(i in seq(cancer_type)){
    dim_reduction=dim_reduction_list[[cancer_type[i]]]
    if(method=='pca'){
      X=dim_reduction$x
      
      sdev = dim_reduction$sdev
      variances = sdev^2
      total_variance = sum(variances)
      explained_variance_percentage = (variances / total_variance*100) %>% signif(3)
      xlab_=paste('PC1 (',explained_variance_percentage[1],'%)')
      ylab_=paste('PC2 (',explained_variance_percentage[2],'%)')
    } 
    if(method=='tsne'){
      X=dim_reduction$Y
      colnames(X)=c('Dimension 1','Dimension 2') 
      rownames(X)=colnames(new_data.final)[meta_in_pub$oncotree==cancer_type[[i]]]
      xlab='Dimension 1'; ylab='Dimension 2'
      xlab_=''; ylab_=''
    } 
    originator_ind=meta_in_pub$sample[meta_in_pub$passage=='Originator' &
                                           meta_in_pub$oncotree==cancer_type[i]]
    sample_ind=meta_in_pub$sample[meta_in_pub$passage!='Originator' &
                                       meta_in_pub$oncotree==cancer_type[i]]
    
    plot(X[,xlab],X[,ylab], col = "#FFFFFF",pch = 21, 
         cex = cex,main=cancer_type[i],
         cex.lab = cex.lab,
         xlab=xlab_,
         ylab=ylab_,
         xaxt='n',yaxt='n')
    
    points(X[sample_ind,xlab],X[sample_ind,ylab], col = "#00000099", bg=cancer_type_cols[[i]], pch = 21, cex = cex)
    points(X[originator_ind,xlab],X[originator_ind,ylab], col = "blue", bg=cancer_type_cols[[i]], pch = 21, cex = cex+1)
  }
}


pdf('figures_pdf/fig_3E.pdf')
plot_all_dim_reduction(dim_reduction_list) 
dev.off()

pdf('figures_pdf/fig_3E_without_ComBat.pdf')
plot_all_dim_reduction(no_combat_dim_reduction_list)
dev.off()

#==============================================================================================================
library(RColorBrewer)

generate_colors = function(n) {
  # Get all palettes from RColorBrewer
  palettes <- brewer.pal.info
  
  # Filter palettes that are qualitative
  qualitative_palettes <- rownames(palettes[palettes$category == "qual", ])
  
  # Extract colors from these palettes
  brewer_colors <- unlist(lapply(qualitative_palettes, function(pal) {
    brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
  }))
  
  # If n is less than or equal to the length of brewer_colors, return the first n colors
  if (n <= length(brewer_colors)) {
    return(brewer_colors[1:n])
  } else {
    # Generate random colors for the remainder
    random_colors_count <- n - length(brewer_colors)
    random_colors <- replicate(random_colors_count, {
      rgb(runif(1), runif(1), runif(1))
    })
    
    return(c(brewer_colors, random_colors))
  }
}


plot_all_dim_reduction_diff_models=function(dim_reduction_list){
  cex=1
  cex.axis = 1.1
  par(mfrow=c(3,3),mar=c(2,2,2,2))
  
  for(i in seq(cancer_type)){
    dim_reduction=dim_reduction_list[[cancer_type[i]]]
    if(method=='pca'){
      X=dim_reduction$x
      
      sdev = dim_reduction$sdev
      variances = sdev^2
      total_variance = sum(variances)
      explained_variance_percentage = (variances / total_variance*100) %>% signif(3)
      xlab_=paste('PC1 (',explained_variance_percentage[1],'%)')
      ylab_=paste('PC2 (',explained_variance_percentage[2],'%)')
    } 
    if(method=='tsne'){
      X=dim_reduction$Y
      colnames(X)=c('Dimension 1','Dimension 2') 
      rownames(X)=colnames(new_data.final)[meta_in_pub$oncotree==cancer_type[i]]
      xlab='Dimension 1'
      ylab='Dimension 2'
      xlab_=''
      ylab_=''
    } 
    
    model_ind=unique(meta_in_pub$model[meta_in_pub$oncotree==cancer_type[i]])
    plot(X[,xlab],X[,ylab], col = "#FFFFFF",pch = 21, cex = cex,main=paste(cancer_type[i],'(n model=',length(model_ind),')'),
         xlab=xlab_,ylab=ylab_,
         cex.axis = cex.axis,
         xaxt='n',yaxt='n')
    colors=generate_colors(length(model_ind)) #try to color each model
    j=1
    for(model in model_ind){
      model_tab=meta_in_pub[meta_in_pub$model==model,]
      
      if('Originator' %in% model_tab$passage){
        originator_ind=model_tab$sample[model_tab$passage=='Originator']
        points(X[originator_ind,xlab],X[originator_ind,ylab], col = "blue", bg=colors[j], pch = 21, cex = cex+1)
      } 
      
      sample_ind=model_tab$sample[model_tab$passage!='Originator']
      points(X[sample_ind,xlab],X[sample_ind,ylab], col = "#00000099", bg=colors[j], pch = 21, cex = cex)
      j=j+1
    }
  }
}


pdf('figures_pdf/sup_fig_8_ComBat.pdf')
set.seed(100)
plot_all_dim_reduction_diff_models(dim_reduction_list)
dev.off()
set.seed(100)
plot_all_dim_reduction_diff_models(no_combat_dim_reduction_list)
#==============================================================================================================



#Boxplots and statistical analysis
get_ind2=function(passage){
  meta_in_pub$sample[meta_in_pub$passage==passage]
}

ind_originator=get_ind2('Originator')
ind_P0=get_ind2('P0'); ind_P1=get_ind2('P1'); ind_P2=get_ind2('P2'); ind_P3=get_ind2('P3'); ind_P4=get_ind2('P4'); ind_P5=get_ind2('P5'); ind_P6=get_ind2('P6')


#pairwise euclidean distance for hypothsis testing
calculate_pairwise_distance=function(originator_x_y,Px_x_y){
  P0_to_originator_dist_num=c()
  for(i in 1:nrow(originator_x_y)){
    if(nrow(Px_x_y)==0){ #if no samples in the passage
      break
    }  
    
    for(j in 1:nrow(Px_x_y)){
      P0_to_originator_dist_num=c(P0_to_originator_dist_num,
                                  sqrt(sum((Px_x_y[j,]-originator_x_y[i,])^2)))
    }
    return(P0_to_originator_dist_num)
  }
}


par(mfrow=c(7,1),mar = c(0.3, 2, 0.3, 2))
for(cancer_type_ in cancer_type){
  ind_cancer_type=meta_in_pub$sample[meta_in_pub$oncotree==cancer_type_]
  
  ##originator itself (negative control)
  originator_x_y=dim_reduction_list[[cancer_type_]]$x[intersect(ind_cancer_type,ind_originator),c(xlab,ylab)]
  originator_dist_num=originator_x_y %>% dist(method='euclidean') %>% as.numeric
  
  Px_dist_num_list=list(originator_dist_num)
  for(passage in list('P0','P1','P2',paste0('P',as.character(3:6)))){
    ind_Px=meta_in_pub$sample[meta_in_pub$passage %in% passage]
    Px_x_y=dim_reduction_list[[cancer_type_]]$x[intersect(ind_cancer_type,ind_Px),c(xlab,ylab),drop=F]
    Px_dist_num=calculate_pairwise_distance(originator_x_y,Px_x_y)
    
    Px_dist_num_list=c(Px_dist_num_list,list(Px_dist_num))
  }
  
  boxplot(Px_dist_num_list,names=rep('',length(Px_dist_num_list))) #8 empty names
  
  list_for_test=Px_dist_num_list[sapply(Px_dist_num_list,FUN = function(element) !is.null(element)) ]
  #print(kruskal.test(list_for_test[2:length(list_for_test)])) #test for H0: P0-P6 all have the same median
  #print(kruskal.test(list_for_test)) #test for H0: originator, P0-P6 all have the same median
  
  print(wilcox.test(originator_dist_num, unlist(list_for_test[2:length(list_for_test)]))) #test for H0: originator vs. (P0-P6 as a whole) are the same
}






