library(fgsea)
library(estimate) #https://bioinformatics.mdanderson.org/public-software/estimate/

#create an empty dataframe/matrix
create_new_tab=function(rownames_=NULL,colnames_=NULL,type='data.frame',fill=NA){
  tab=matrix(fill,nrow=length(rownames_),ncol=length(colnames_),dimnames=list(rownames_,colnames_))
  if(type=='data.frame') tab=as.data.frame(tab)
  return(tab)
}


#return the replicated name(s) from a vector and their frequencies
check_replicates=function(x){
  out=table(x)
  out=out[out>1]
  
  return(out)
}

#convert all values in a dataframe to the same type (eg. character, numeric)
#ref: https://stackoverflow.com/questions/26391921/how-to-convert-entire-dataframe-to-numeric-while-preserving-decimals
convert_val_type=function(X,type){
  
  X=as.data.frame(X)
  
  if(type=='numeric') func=as.numeric
  if(type=='character') func=as.character
  
  new_X=data.frame(lapply(X, func), stringsAsFactors = FALSE)
  rownames(new_X)=rownames(X) #resolve the disappearing rowname issue
  colnames(new_X)=colnames(X) #resolve the altering colname issue
  
  return(new_X)
}


##impute expression by medians
impute_expression <- function(mat_) {
  mat_ <- t(apply(mat_, 1, function(row) {
    row[is.na(row)] <- median(row, na.rm = TRUE) #impute by median
    return(row)
  }))
  return(mat_)
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



#The following is directly copied and modified from filterCommonGenes() and estimateScore() (did this so the input can accept a df rather than a file, and I can get the actual tumor fraction table)
my_ESTIMATE=function(gene_expression_matrix){
  id = c("GeneSymbol")
  input.df <- as.data.frame(gene_expression_matrix)
  merged.df <- merge(common_genes, input.df, by.x = id, by.y = "row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).", 
                nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
  
  #let this return the score df #(!)not sure if platform='illumina' stands for RNAseq, microarray, or both
  my_estimateScore=function (merged.df, platform = c("affymetrix", "agilent","illumina")){
    platform <- match.arg(platform)
    ds=cbind(Description=rownames(merged.df),merged.df)
    descs <- ds[, 1]
    ds <- ds[-1]
    row.names <- row.names(ds)
    names <- names(ds)
    dataset <- list(ds = ds, row.names = row.names, descs = descs, 
                    names = names)
    m <- data.matrix(dataset$ds)
    gene.names <- dataset$row.names
    sample.names <- dataset$names
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    
    for (j in 1:Ns) {
      m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m/Ng
    gs <- as.matrix(SI_geneset[, -1], dimnames = NULL)
    N.gs <- 2
    gs.names <- row.names(SI_geneset)
    score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
    for (gs.i in 1:N.gs) {
      gene.set <- gs[gs.i, ]
      gene.overlap <- intersect(gene.set, gene.names)
      print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", 
                  length(gene.overlap)))
      if (length(gene.overlap) == 0) {
        score.matrix[gs.i, ] <- rep(NA, Ns)
        next
      }
      else {
        ES.vector <- vector(length = Ns)
        for (S.index in 1:Ns) {
          gene.list <- order(m[, S.index], decreasing = TRUE)
          gene.set2 <- match(gene.overlap, gene.names)
          correl.vector <- m[gene.list, S.index]
          TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
          no.TAG <- 1 - TAG
          N <- length(gene.list)
          Nh <- length(gene.set2)
          Nm <- N - Nh
          correl.vector <- abs(correl.vector)^0.25
          sum.correl <- sum(correl.vector[TAG == 1])
          P0 <- no.TAG/Nm
          F0 <- cumsum(P0)
          Pn <- TAG * correl.vector/sum.correl
          Fn <- cumsum(Pn)
          RES <- Fn - F0
          max.ES <- max(RES)
          min.ES <- min(RES)
          if (max.ES > -min.ES) {
            arg.ES <- which.max(RES)
          }
          else {
            arg.ES <- which.min(RES)
          }
          ES <- sum(RES)
          EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, 
                                  RES = RES, indicator = TAG)
          ES.vector[S.index] <- EnrichmentScore$ES
        }
        score.matrix[gs.i, ] <- ES.vector
      }
    }
    score.data <- data.frame(score.matrix)
    names(score.data) <- sample.names
    row.names(score.data) <- gs.names
    estimate.score <- apply(score.data, 2, sum)
    if (platform != "affymetrix") {
      score.data <- rbind(score.data, estimate.score)
      rownames(score.data) <- c("StromalScore", "ImmuneScore", 
                                "ESTIMATEScore")
    }
    else {
      convert_row_estimate_score_to_tumor_purity <- function(x) {
        stopifnot(is.numeric(x))
        cos(0.6049872018 + 0.0001467884 * x)
      }
      est.new <- NULL
      for (i in 1:length(estimate.score)) {
        est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
        est.new <- rbind(est.new, est_i)
        if (est_i >= 0) {
          next
        }
        else {
          message(paste(sample.names[i], ": out of bounds", 
                        sep = ""))
        }
      }
      colnames(est.new) <- c("TumorPurity")
      estimate.t1 <- cbind(estimate.score, est.new)
      x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 
        0
      estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
      score.data <- rbind(score.data, t(estimate.t1))
      rownames(score.data) <- c("StromalScore", "ImmuneScore", 
                                "ESTIMATEScore", "TumorPurity")
    }
    
    return(t(score.data))
  }
  tumor_frac_matrix=my_estimateScore(merged.df) #(!)not sure if platform='illumina' stands for RNAseq, microarray, or both
  
  return(tumor_frac_matrix)
}


