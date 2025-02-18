library(oncoEnrichR)
library(circlize)
library(gridExtra)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)
library(DESeq2)
library(fgsea)

setwd('./') #set this to the working directory
source('./self_defined_functions.R')



sup_tab=readxl::read_excel('Supp_tables_v27.xlsx',sheet='Supp_Table_18',col_names = T,skip = 1)


sup_tab$`EFSx4 for Erlotinib`[sup_tab$`EFSx4 for Erlotinib`=='>1.53']=1.53 #put a star in the ppt indicating this is not actually 1.53. Rather, it's >1.53
sup_tab=sup_tab[order(as.numeric(sup_tab$`EFSx4 for Erlotinib`),decreasing=T),]
EFSx4=ifelse(as.numeric(sup_tab$`EFSx4 for Erlotinib`)>1.5,'EFSx4 > 1.5','EFSx4 <= 1.5')

#Load RNA-seq data
load('./mutational_landscape_normalizedCount.RData') #loaded as: data.final
data.final_original=data.final
colnames(data.final_original)[grepl('_FFPE',colnames(data.final_original))]=sub('_FFPE$','',colnames(data.final_original)[grepl('_FFPE',colnames(data.final_original))]) ##correct some originator names
colnames(data.final_original)=gsub('~','-',colnames(data.final_original))

load('./batch_adjusted_mutational_landscape_normalizedCount.RData')
colnames(data.final)=gsub('~','-',colnames(data.final))
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))]) ##correct some originator names

#compute the meta table
##### 
if(identical(class(data.final),'data.frame')) data.final=as.matrix(data.final)
colnames(data.final)=gsub('~','-',colnames(data.final))
##correct some originator names
colnames(data.final)[grepl('_FFPE',colnames(data.final))]=sub('_FFPE$','',colnames(data.final)[grepl('_FFPE',colnames(data.final))])

patientID_modelID_sampleID_passage=readxl::read_excel('Supp_tables_v21.xlsx',sheet='Supp_table_3',skip=1)[,c("Patient ID","Specimen/Parent model ID","Sample ID","PDX Passage or PDC or PDOrg or Patient Orginator specimen")] %>% as.data.frame
patientID_modelID_sampleID_passage=patientID_modelID_sampleID_passage #I don't limit the passage no. here. Otherwise, some models will have no samples at all.

meta=patientID_modelID_sampleID_passage
colnames(meta)=c('patient','model','sample','passage')  
#####


dat_list=list()
for(data_ in list(data.final_original,data.final)){
  dat_sample=data_[,meta$sample] 
  
  #aggregate the data to the model level
  n_samples_for_each_model=c()
  dat=create_new_tab(rownames(dat_sample),unique(meta$model[meta$passage!='Originator']))
  dat=dat[,colnames(dat) %in% sup_tab$`PDX Model`]
  
  for(model in colnames(dat)){
    sample_in_model=meta$sample[meta$model==model & meta$passage!='Originator']
    n_samples_for_each_model=c(n_samples_for_each_model,
                               length(sample_in_model))
    if(length(sample_in_model)==1){
      dat[,model]=dat_sample[,sample_in_model]
    }else{
      dat[,model]=rowMeans(dat_sample[,sample_in_model])
    }
  }
  
  #this line includes the originator data:
  #dat=cbind(dat,dat_sample[,(meta$model %in% sup_tab$`PDX Model`) & (meta$passage=='Originator')])
  
  dat_list=c(dat_list,
             list(dat))
}

model_dat_original=dat_list[[1]] #originator data are not combined into the model data
model_dat=dat_list[[2]] #originator data are not combined into the model data





#setdiff(sup_tab$`PDX Model`,colnames(model_dat)) #all models are present

models_have_RNAseq=sup_tab$`PDX Model`
EFSx4=sapply(models_have_RNAseq,function(model) EFSx4[sup_tab$`PDX Model`==model])
histology=sapply(models_have_RNAseq,function(model) sup_tab$`Tumor histology (Oncotree code)`[sup_tab$`PDX Model`==model])
EGFR_amplification=ifelse(sup_tab$`EGFR Amplification (Copy number >=7)`=='Yes','Yes','No') #Use the last column
condition_df=data.frame(EFSx4,'EGFR amplification'=EGFR_amplification,'Tumor histology'=histology,check.names=F)
  
DE=perform_my_DE(round(model_dat[,models_have_RNAseq]),condition=EFSx4,baseline='EFSx4 <= 1.5')

#all DE genes

#plot EFSx4 barplot + gene expression heatmap (Z-score after log transformation)
quantitative_EFSx4=as.numeric(sup_tab$`EFSx4 for Erlotinib`)


pdf('figures_pdf/fig_6A.pdf',width = 4.3*2,height = 3*2)
DE_expression_Z_score=log2(model_dat[rownames(DE),models_have_RNAseq]+1) %>% apply(MARGIN=1,FUN=scale) %>% t
colnames(DE_expression_Z_score)=models_have_RNAseq; rownames(DE_expression_Z_score)=rownames(DE)
plot_my_heatmap(DE_expression_Z_score[rownames(DE),],
                my_annotation_columns = condition_df,
                annotation_colors=list('EGFR amplification'=c(Yes='#c40a35',No='#03a1fc'),
                                       EFSx4=c('EFSx4 > 1.5'='red','EFSx4 <= 1.5'='blue'),
                                        'Tumor histology'=c(BLCA='#f49b7f',COADREAD='#e64b35',HNSC='#0aa087',NSCLC='#3d5488',MEL='#8391b4',PAAD='#8e0153',SARCOMA='#7d6046',Other='#b09c85')),
                cluster_rows = T,fontsize=6,
                color=NA)
dev.off()





#GSEA for all the genes, not just DE genes
DE_complete=perform_my_DE(round(model_dat[,models_have_RNAseq]),condition=condition_df$EFSx4,baseline='EFSx4 <= 1.5',padj=1,fold_change_cutoff=1)

pdf('figures_pdf/sup_fig_10_A.pdf',width=10,height=10)
GSEA_hallmark=perform_my_GSEA(DE_complete,'./h.all.v2022.1.Hs.symbols.gmt') #this file can be downloaded from GSEA website
GSEA_hallmark
dev.off()

# Plot enrichment for the selected pathways

gene_stats=DE_complete$stat #this calculates the ranking
names(gene_stats)=rownames(DE_complete)
gene_stats = sort(gene_stats, decreasing = TRUE)

GSEA_gene_list=fgsea::gmtPathways('./h.all.v2022.1.Hs.symbols.gmt') 
GSEA_gene_list=GSEA_gene_list[GSEA_hallmark$Results$pathway]

pdf('figures_pdf/sup_fig_10B.pdf',width = 9,height = 12)
plots <- list()
for(pathway in names(GSEA_gene_list)){
  plot_=plotEnrichment(pathway = GSEA_gene_list[[pathway]], stats = gene_stats) + labs(title = pathway)
  plot_$layers[[1]]$aes_params$colour='black'
  plot_$layers[[3]]$aes_params$colour='red'
  plot_$layers[[4]]$aes_params$colour='blue'
  plot_$labels$x=''
  plot_$labels$y=''
  
  plots[[pathway]]=plot_
}
do.call(grid.arrange, c(plots, ncol=2))
dev.off()


#GSEA for EGFR and AKT-related pathways
DE_result=DE_complete
myGO=c(fgsea::gmtPathways('./biocarta_EGFR_related_genesets.v2023.2.Hs.gmt'),
       fgsea::gmtPathways('./biocarta_AKT_related_genesets.v2023.2.Hs.gmt'),
       fgsea::gmtPathways('./biocarta_pathways_for_low_RMEFS.gmt')
       )
myGO=myGO[names(myGO) %>% unique] #remove duplicated pathways
pval=0.05

set.seed(1230)
gene_stats=DE_result$stat #this calculates the ranking
names(gene_stats)=rownames(DE_result)
gene_stats = sort(gene_stats, decreasing = TRUE)

fgRes = fgsea::fgsea(pathways = myGO,
                     stats = gene_stats) %>% 
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

#message("Collapsing Pathways -----")
concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                    pathways = myGO,
                                    stats = gene_stats) #ref: https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")

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
  theme(axis.text.y = element_text(size = 14))+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=header) 

pdf('figures_pdf/fig_6B.pdf',width = 3*2.5,height = 4*2)
g1
dev.off()

output = list("Results" = fgRes, "Plot" = g1)

output
GSEA_EGFR_related=output

GSEA_gene_list=myGO[output$Results$pathway]

genes=unique(unlist(GSEA_gene_list))
pathways=names(GSEA_gene_list)
gene_union_df=create_new_tab(rownames_=genes,colnames_=pathways)
for(i in 1:nrow(gene_union_df)){
  for(j in 1:ncol(gene_union_df)){
    gene_union_df[i,j]=ifelse(genes[i] %in% GSEA_gene_list[[j]],1,0)
  }
}

color_mapping <- colorRamp2(c(0, 1), c("white", "blue"))
ComplexHeatmap::pheatmap(as.matrix(gene_union_df),fontsize_row = 5,fontsize_col=8,cluster_rows = F,cluster_cols = F,color=color_mapping,
                         border_color='grey60',
                         name=' ')


# Plot each enriched pathway
pdf('figures_pdf/fig_6C.pdf',width = 3*2,height = 4*2)
plots <- list()
for(pathway in names(GSEA_gene_list)){
  plot_=plotEnrichment(pathway = GSEA_gene_list[[pathway]], stats = gene_stats) + labs(title = pathway)
  plot_$layers[[1]]$aes_params$colour='black'
  plot_$layers[[3]]$aes_params$colour='red'
  plot_$layers[[4]]$aes_params$colour='blue'
  plot_$labels$x=''
  plot_$labels$y=''
  
  plots[[pathway]]=plot_
}
do.call(grid.arrange, c(plots, ncol=1))
dev.off()



#GO enrichment
ego_list=list()
up_regulated_gene=rownames(DE)[DE$log2FoldChange>0]
down_regulated_gene=rownames(DE)[DE$log2FoldChange<0]
for(ontology in c("BP",'MF','CC')){
  for(gene in list(up_regulated_gene,down_regulated_gene)){

    ego <- enrichGO(gene=gene,
                    OrgDb='org.Hs.eg.db',  # Choose the appropriate OrgDb for your organism
                    keyType="SYMBOL",      # Change if using a different identifier type
                    ont=ontology,        
                    pAdjustMethod="fdr",         # Method for adjusting p values
                    readable=TRUE)          # Convert Entrez IDs to gene names in the output
    
    ego_list=c(ego_list,
               list(ego))
  }
}


#BP for up- and down-regulated genes
dotplot(ego_list[[1]])
dotplot(ego_list[[2]])

#MF for up- and down-regulated genes
#dotplot(ego_list[[3]])
#dotplot(ego_list[[4]]) #no enriched term

#CC for up- and down-regulated genes
#dotplot(ego_list[[5]]) #not very useful
#dotplot(ego_list[[6]]) #no enriched term


write.csv(DE_complete,file='sup_table_16_DE_all_genes.csv',quote=F)

temp=GSEA_hallmark$Results
temp$leadingEdge=sapply(temp$leadingEdge,function(gene_vector) paste(gene_vector,collapse=';') )
write.csv(temp,file='sup_table_17_GSEA_hallmark.csv',quote=F,row.names=F)

temp=GSEA_EGFR_related$Results
temp$leadingEdge=sapply(temp$leadingEdge,function(gene_vector) paste(gene_vector,collapse=';') )
write.csv(temp,file='sup_table_18_GSEA_EGFR_related.csv',quote=F,row.names=F)