DE=read.csv('./data/erlotinib_response_DE_genes.csv',row.names=1)
DE_genes=rownames(DE)
DE_genes

#Manually paste the above gene list to OncoEnrich Web-based version. Supplementary figure 9 is based on the String interaction network within the Oncoenrich result.
#Got to: https://oncotools.elixir.no
#Click 'oncoEnrichR'. Then paste the genes, separated by new lines.


