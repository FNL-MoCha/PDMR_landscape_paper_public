library(circlize)
library(gridExtra)
library(ggplot2)
library(clusterProfiler)
library(readxl)
library(ComplexHeatmap)
library(DESeq2)

source('./sup_code/self_defined_functions.R')

load('./data/erlotinib_model_expression.RData') #loaded as cts
load('./data/erlotinib_condition_df.RData') #loaded as condition_df 


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


pdf('Extended_Data_Fig_9A.pdf',width=10,height=10)
GSEA_hallmark=perform_my_GSEA(DE_complete,'~/Documents/MoCha/util/h.all.v2022.1.Hs.symbols.gmt')
GSEA_hallmark
dev.off()



# Plot enrichment for the selected pathways
gene_stats=DE_complete$stat #this calculates the ranking
names(gene_stats)=rownames(DE_complete)
gene_stats = sort(gene_stats, decreasing = TRUE)

GSEA_gene_list=fgsea::gmtPathways('~/Documents/MoCha/util/h.all.v2022.1.Hs.symbols.gmt')
GSEA_gene_list=GSEA_gene_list[GSEA_hallmark$Results$pathway]

pdf('Extended_Data_Fig_9B.pdf',width = 10,height = 12)
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
