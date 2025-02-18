#create an empty dataframe/matrix
create_new_tab=function(rownames_=NULL,colnames_=NULL,type='data.frame',fill=NA){
  tab=matrix(fill,nrow=length(rownames_),ncol=length(colnames_),dimnames=list(rownames_,colnames_))
  if(type=='data.frame') tab=as.data.frame(tab)
  return(tab)
}


#plot heatmap by pheatmap or ComplexHeatmap
plot_my_heatmap=function(X,
                         my_annotation_columns=NA, #a column name-named vector corresponding to colnames(X). A dataframe can also be used to add more than 1 annotation.
                         my_annotation_rows=NA, #a row name-named vector corresponding to rownames(X). A dataframe can also be used to add more than 1 annotation.
                         color=colorRampPalette(c("#4575B4", "#BACCE3","white","#F0B2AF","#D73027"))( 100 ), #if given NA, this will automatically be assigned
                         border_color='grey60',
                         cluster_rows=F, #replace default value
                         cluster_cols=F, #replace default value
                         center_value=NA, #make this the center of the color scale
                         unequal_portion=F, # if a center value is given and this is T, assign full spectrum of color but assign center color for the center value
                         name=' ', #this removes the undesired annotation on top of the heatmap color bar. It doesn't take empty string, so I give a blank character
                         ...){
  
           # Code to execute if module is 'ComplexHeatmap::pheatmap'
           print("- ComplexHeatmap::pheatmap is called - ")
           
           #column annotations and colors
           #--------------------------------------------------------
           if(!identical(NA,my_annotation_columns)){
             annotation_col=data.frame(my_annotation_columns)
             if(!is.null(names(my_annotation_columns)) & !is.vector(my_annotation_columns)){
               colnames(annotation_col)=names(my_annotation_columns)
             }else{
               colnames(annotation_col)='group'
             }
             rownames(annotation_col)=colnames(X)
           }else{
             annotation_col=my_annotation_columns
           }
           
           #--------------------------------------------------------
           
           #row annotations and colors
           #--------------------------------------------------------
           if(!identical(NA,my_annotation_rows)){
             annotation_row=data.frame(my_annotation_rows)
             if(!is.null(names(my_annotation_rows)) & !is.vector(my_annotation_rows)){
               colnames(annotation_row)=names(my_annotation_rows)
             }else{
               colnames(annotation_row)='group2'
             }
             rownames(annotation_row)=rownames(X)
           }else{
             annotation_row=my_annotation_rows
           }
           
           #--------------------------------------------------------
           
           if(identical(NA,color)){
             if(!is.na(center_value)){
               n_color_density=1001
               data_range=range(as.numeric(as.matrix(X)))
               
               p_greater_than_value=(data_range[2]-center_value)/(data_range[2]-data_range[1])
               p_less_equal_than_value=1-p_greater_than_value
               
               half=(n_color_density-1)/2
               color=colorRampPalette(c(rev(brewer.pal(n = 7, name ="RdYlBu"))))(n_color_density)
               
               if(unequal_portion){
                 if(p_greater_than_value>p_less_equal_than_value){
                   color=color[c(seq(from=1,to=half,by=p_greater_than_value/p_less_equal_than_value),
                                 half+1,
                                 (half+2):n_color_density)]
                 }else{
                   color=color[c(1:half,
                                 half+1,
                                 seq(from=half+2,to=n_color_density,by=p_less_equal_than_value/p_greater_than_value))]
                 }
               }else{
                 if(p_greater_than_value>p_less_equal_than_value){
                   
                   
                   color=color[c((half-round(half*(p_less_equal_than_value/p_greater_than_value))):half,
                                 half+1,
                                 (half+2):n_color_density)]
                 }else{
                   color=color[c(1:half,
                                 half+1,
                                 (half+2):(half+2+round(half*(p_greater_than_value/p_less_equal_than_value))))]
                 }
               }
             }else{ #the default values
               color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
               breaks=seq(100)
             }
           }
           
           if('data.frame' %in% class(X)) X=as.matrix(X)
           
           ComplexHeatmap::pheatmap(X,
                                    cluster_rows=cluster_rows,
                                    cluster_cols=cluster_cols,
                                    annotation_col=annotation_col,
                                    annotation_row=annotation_row,
                                    color=color,
                                    name=name,
                                    border_color=border_color,
                                    ...)
         
}



#For "condition", the first entry will be treated as the baseline (eg. If c( "ctr", "ctr", "trt", "trt"), "ctr" is treated as the baseline). 
#(This is coded in the middle of this function)
perform_my_DE=function(gene_expression_matrix=NULL,txi.rsem=NULL,condition,baseline,padj=0.05,fold_change_cutoff=2,reads_cutoff=10,input_type='df'){
  
  if(input_type=='tximport'){
    
    coldata=data.frame(condition = as.factor(condition))
    rownames(coldata)=colnames(txi.rsem$counts)
    
    dds=DESeqDataSetFromTximport(txi.rsem, coldata, ~condition)
    
  }else if(input_type=='df'){
    
    cts=as.matrix(gene_expression_matrix)
    
    coldata=data.frame(condition = as.factor(condition))
    rownames(coldata)=colnames(cts)
    
    
    dds=DESeqDataSetFromMatrix(countData = cts,
                               colData = coldata,
                               design= ~ condition)
  }
  
  
  
  #(Different between studies) remove the genes having < [reads_cutoff] reads 
  dds=dds[rowSums(counts(dds)) >= reads_cutoff,]
  
  # set the baseline as the reference
  dds$condition=relevel(dds$condition, ref = baseline) 
  
  
  # Perform differential gene expression analysis
  dds=DESeq(dds)
  
  res=results(dds)
  res=as.data.frame(res[order(res$log2FoldChange),] )
  res=res[!is.na(res$padj) & (res$padj<=padj) & (abs(res$log2FoldChange)>=log2(fold_change_cutoff)),]
  
  DE_result=res
  
  return(DE_result)
}


#The "DE_result" input should be one of the DE result matrix or a stats (eg. T-test) vector whose names are genes 
perform_my_GSEA=function(DE_result, pathway_or_GO_gmt_file_path,pval=0.05,seed=1230,n_up=NA,n_down=NA) { 
  ##If DE_result is a DE table, ranking_method can be based on stat, log2FoldChange...etc. They must be columns in the DE object. 
  #n_up, n_down subselect the top up- or down-regulated pathways
  
  set.seed(seed)
  my_gmt = fgsea::gmtPathways(pathway_or_GO_gmt_file_path)
  
  if((class(DE_result) %in% c('matrix','data.frame')) && ncol(DE_result)>=2 ){ #if this check failed, the remaining possibility must be that DE_result is a ranked gene vector
    gene_ranking=DE_result[,'stat'] #this calculates the ranking
    names(gene_ranking)=rownames(DE_result)
    gene_ranking = sort(gene_ranking, decreasing = TRUE)
    
  }else{
    gene_ranking=DE_result
  }
  
  
  fgRes = fgsea::fgsea(pathways = my_gmt,
                       stats = gene_ranking,
                       scoreType='std') %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES)) 
  
  #(!) Some warning message
  
  "
    Warning message:
    In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
      There are ties in the preranked stats (1.42% of the list).
    The order of those tied genes will be arbitrary, which may produce unexpected results.
    
    => Seems I can ignore the error: https://github.com/YuLab-SMU/clusterProfiler/issues/214
    "
  
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),pathways = my_gmt,stats = gene_ranking) #this collapses pathways having too many overlapped genes
  
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  if(!is.na(n_up) & !is.na(n_down)){
    fgRes = rbind(head(fgRes, n = n_up),
                  tail(fgRes, n = n_down ))
  }
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("(Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(fgRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) 
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

