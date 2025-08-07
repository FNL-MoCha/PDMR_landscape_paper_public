load('../data/expression_no_organoid_no_PDC.RData')
load('../data/PDX_TCGA_CCLE_inds.RData')

gene_id <- rownames(expression_no_organoid_no_PDC)
hists   <- c("BLCA","COAD","READ","HNSC","MEL","LUAD","LUSC","PAAD")


for (h in hists) {
  idx <- PDX_TCGA_CCLE_inds$PDX[[paste0(h, "_originator_ind")]]
  if (is.null(idx) || length(idx) == 0) next  # nothing to write
  
  out <- data.frame(
    gene_id = gene_id,
    expression_no_organoid_no_PDC[, idx, drop = FALSE],
    check.names = FALSE
  )
  
  write.table(
    out,
    file = sprintf("%s_originator_expression.tsv", h),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
