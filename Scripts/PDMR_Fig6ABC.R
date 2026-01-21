library(circlize)
library(gridExtra)
library(ggplot2)
library(clusterProfiler)
library(readxl)
library(ComplexHeatmap)
library(DESeq2)

source('./sup_code/self_defined_functions.R')

sup_tab=readxl::read_excel('./data/erlotinib_sensitivity_models.xlsx',,col_names = T)
sup_tab$`EFSx4 for Erlotinib`[sup_tab$`EFSx4 for Erlotinib`=='>1.53']=1.53 #put a star in the ppt indicating this is not actually 1.53. Rather, it's >1.53
sup_tab=sup_tab[order(as.numeric(sup_tab$`EFSx4 for Erlotinib`),decreasing=T),]
EFSx4=ifelse(as.numeric(sup_tab$`EFSx4 for Erlotinib`)>1.5,'EFSx4 > 1.5','EFSx4 <= 1.5')

#Load RNA-seq data
load('./data/DESeq2_mutational_landscape.RData') #load this for DESeq2. Loaded as "ddsColl2".
raw_counts=counts(ddsColl2,normalized=F)
colnames(raw_counts)=sub('~RNASEQ','',colnames(raw_counts))
colnames(raw_counts)[grepl('_FFPE',colnames(raw_counts))]=sub('_FFPE$','',colnames(raw_counts)[grepl('_FFPE',colnames(raw_counts))]) ##correct some originator names
colnames(raw_counts)=gsub('~','-',colnames(raw_counts))

#load the meta table
load('./data/meta_in_pub.RData')

data_=raw_counts
dat_sample=data_[,meta_in_pub$sample] 

#aggregate the data to the model level using summation
n_samples_for_each_model=c()
dat=create_new_tab(rownames(dat_sample),unique(meta_in_pub$model[meta_in_pub$passage!='Originator']))
dat=dat[,colnames(dat) %in% sup_tab$`PDX Model`]

for(model in colnames(dat)){
  sample_in_model=meta_in_pub$sample[meta_in_pub$model==model & meta_in_pub$passage!='Originator']
  n_samples_for_each_model=c(n_samples_for_each_model,
                             length(sample_in_model))
  if(length(sample_in_model)==1){
    dat[,model]=dat_sample[,sample_in_model]
  }else{
    dat[,model]=rowSums(dat_sample[,sample_in_model]) #take sum like what DESeq2::collapseReplicates() would do 
  }
}

model_dat_original=dat




#setdiff(sup_tab$`PDX Model`,colnames(model_dat_original)) #all models are present

models_have_RNAseq=sup_tab$`PDX Model`
EFSx4=sapply(models_have_RNAseq,function(model) EFSx4[sup_tab$`PDX Model`==model])
names(EFSx4)=models_have_RNAseq
histology=sapply(models_have_RNAseq,function(model) sup_tab$`Tumor histology (Oncotree code)`[sup_tab$`PDX Model`==model])
EGFR_amplification=ifelse(sup_tab$`EGFR Amplification (Copy number >=7)`=='Yes','Yes','No') #Use the last column
condition_df=data.frame(EFSx4,'EGFR amplification'=EGFR_amplification,'Tumor histology'=histology,check.names=F)
condition_df$EFSx4=as.factor(condition_df$EFSx4)
condition_df$EFSx4=relevel(condition_df$EFSx4,ref='EFSx4 <= 1.5')
save(condition_df,file='./data/erlotinib_condition_df.RData')

#perform DE
cts=as.matrix(round(model_dat_original[,models_have_RNAseq]))
save(cts,file='./data/erlotinib_model_expression.RData')
condition=EFSx4
baseline='EFSx4 <= 1.5'
reads_cutoff=10
padj=0.05
fold_change_cutoff=2

coldata=data.frame(condition = as.factor(condition))
rownames(coldata)=colnames(cts)


dds=DESeqDataSetFromMatrix(countData = cts,
                           colData = coldata,
                           design= ~ condition)


#(Different between studies) remove the genes having < [reads_cutoff] reads 
dds=dds[rowSums(counts(dds)) >= reads_cutoff,]

# set the baseline as the reference
dds$condition=relevel(dds$condition, ref = baseline) 


# Perform differential gene expression analysis
dds=DESeq(dds)

res=results(dds)
res=as.data.frame(res[order(res$log2FoldChange),] )
res=res[!is.na(res$padj) & (res$padj<=padj) & (abs(res$log2FoldChange)>=log2(fold_change_cutoff)),]

DE=res
write.csv(DE,file='./data/erlotinib_response_DE_genes.csv',row.names = T)



#all DE genes

#plot EFSx4 barplot + gene expression heatmap (Z-score after log transformation)
quantitative_EFSx4=as.numeric(sup_tab$`EFSx4 for Erlotinib`)



pdf('fig_6A.pdf',width = 4.3*2,height = 3*2)
DE_expression_Z_score=log2(counts(dds,normalized=T)[rownames(DE),]+1) %>% apply(MARGIN=1,FUN=scale) %>% t
colnames(DE_expression_Z_score)=models_have_RNAseq; rownames(DE_expression_Z_score)=rownames(DE)
order_=order(condition_df$EFSx4,condition_df$`EGFR amplification`,condition_df$`Tumor histology`,decreasing=T)
ComplexHeatmap::pheatmap(DE_expression_Z_score[rownames(DE),order_],
                annotation_col = condition_df[order_,],
                annotation_colors=list('EGFR amplification'=c(Yes='#c40a35',No='#03a1fc'),
                                       EFSx4=c('EFSx4 > 1.5'='red','EFSx4 <= 1.5'='blue'),
                                       'Tumor histology'=c(BLCA='#f49b7f',COADREAD='#e64b35',HNSC='#0aa087',NSCLC='#3d5488',MEL='#8391b4',PAAD='#8e0153',SARCOMA='#7d6046',Other='#b09c85')),
                cluster_rows = T,
                cluster_cols = F,
                border_color='grey60',
                fontsize=6)
dev.off()




#GSEA for EGFR and AKT-related pathways

#GSEA for all the genes, not just DE genes
baseline='EFSx4 <= 1.5'
reads_cutoff=10
padj=1
fold_change_cutoff=1
coldata=data.frame(condition = as.factor(condition_df$EFSx4))
rownames(coldata)=colnames(cts)

dds=DESeqDataSetFromMatrix(countData = cts,
                           colData = coldata,
                           design= ~ condition)

#(Different between studies) remove the genes having < [reads_cutoff] reads 
dds=dds[rowSums(counts(dds)) >= reads_cutoff,]

# set the baseline as the reference
dds$condition=relevel(dds$condition, ref = baseline) 


# Perform differential gene expression analysis
dds=DESeq(dds)

res=results(dds)
res=as.data.frame(res[order(res$log2FoldChange),] )
res=res[!is.na(res$padj) & (res$padj<=padj) & (abs(res$log2FoldChange)>=log2(fold_change_cutoff)),]

DE_complete=res


DE_result=DE_complete
myGO=c(fgsea::gmtPathways('./data/biocarta_EGFR_related_genesets.v2023.2.Hs.gmt'),
       fgsea::gmtPathways('./data/biocarta_AKT_related_genesets.v2023.2.Hs.gmt'),
       fgsea::gmtPathways('./data/biocarta_pathways_for_low_RMEFS.gmt')
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

GSEA_EGFR_related=output
GSEA_gene_list=myGO[output$Results$pathway]



# Plot each enriched pathway
pdf('figures_pdf/fig_6C.pdf',width = 3*2.5,height = 4*2)
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
