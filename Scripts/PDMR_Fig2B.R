library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

###############################################################################
#Three input files were compiled from Supp_Table7 (PDX), Supp_Table8 (PDC) 
#and Supp_Table9 (PDOrg). Fusion data was obtained from Supp_Table10
###############################################################################

#pdx
pdx <- read.table("PDMR_Fig2B_PDX.txt",sep = "\t",header = TRUE)
df <- as.data.frame(pdx)

#pdc
pdc <- read.table("PDMR_Fig2B_PDC.txt",sep = "\t",header = TRUE)
df <- as.data.frame(pdc)

#pdorg
pdorg <- read.table("PDMR_Fig2B_PDOrg.txt",sep = "\t",header = TRUE)
df <- as.data.frame(pdorg)


mutations <- df[1:18, ]
amplifications <- df[19:21, ]
deletions <- df[22, ]
fusion <- df[23:29,]
sig = df[30:31,]
pdx_no = df[32,]

###############################################################################
#mutation
rownames(mutations) = mutations$Gene
mut_mat = as.matrix(mutations[,2:9])

cell_fun1 = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "gray", fill = fill))
  if(mut_mat[i, j] != 0) {
    if(mut_mat[i, j] > 9) {
      grid.text(mut_mat[i, j], x, y, gp = gpar(col = "white"))
    } else {
      grid.text(mut_mat[i, j], x, y, gp = gpar(col = "black"))
    }
  }
}

col_fun1 = colorRamp2(c(0, 5, 10, 20, 40), c("white","#bae4b3","#74c476", "#31a354", "#006d2c"))

hm1 <- Heatmap(mut_mat, 
               name = "#_of_LOE_mutations",
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE,
               show_column_names = TRUE,
               column_names_side = "top",
               row_names_side = "left",
               column_names_gp = gpar(fontsize=13,fontface="bold"),
               row_names_gp = gpar(fontsize=10,fontface="bold"),
               cell_fun = cell_fun1,
               col = col_fun1,
               show_heatmap_legend = FALSE,
               width = unit(6, "cm"),
               height = unit(13, "cm")
)
hm1

###############################################################################
#Amp
rownames(amplifications) = amplifications$Gene
amp_mat = as.matrix(amplifications[,2:9])

cell_fun2 = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "gray", fill = fill))
  if(amp_mat[i, j] != 0) {
    if(amp_mat[i, j] > 9) {
      grid.text(amp_mat[i, j], x, y, gp = gpar(col = "white"))
    } else {
      grid.text(amp_mat[i, j], x, y, gp = gpar(col = "black"))
    }
  }
}

col_fun2 = colorRamp2(c(0, 2, 5, 10, 20), c("white","#fee5d9", "#fcae91", "#fb6a4a", "#a50f15"))

hm2 <- Heatmap(amp_mat, 
               name = "#_of_LOE_CNV_Amp",
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE,
               show_column_names = FALSE,
               column_names_side = "top",
               row_names_side = "left",
               column_names_gp = gpar(fontsize=13,fontface="bold"),
               row_names_gp = gpar(fontsize=10,fontface="bold"),
               cell_fun = cell_fun2,
               col = col_fun2,
               show_heatmap_legend = FALSE,
               width = unit(7, "cm"),
               height = unit(2, "cm")
)
hm2

###############################################################################
#Del
rownames(deletions) = deletions$Gene
del_mat = as.matrix(deletions[,2:9])

cell_fun3 = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "gray", fill = fill))
  if(del_mat[i, j] != 0) {
    if(del_mat[i, j] > 9) {
      grid.text(del_mat[i, j], x, y, gp = gpar(col = "white"))
    } else {
      grid.text(del_mat[i, j], x, y, gp = gpar(col = "black"))
    }
  }
}

col_fun3 = colorRamp2(c(0, 2, 5, 10, 20), c("white","#E0FFFF", "#ADD8E6", "#4682B4", "#000080"))

hm3 <- Heatmap(del_mat, 
               name = "#_of_LOE_CNV_Del",
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE,
               show_column_names = FALSE,
               column_names_side = "top",
               row_names_side = "left",
               column_names_gp = gpar(fontsize=13,fontface="bold"),
               row_names_gp = gpar(fontsize=10,fontface="bold"),
               cell_fun = cell_fun3,
               col = col_fun3,
               show_heatmap_legend = FALSE,
               width = unit(7, "cm"),
               height = unit(1, "cm")
)
hm3

###############################################################################
#fusion
rownames(fusion) = fusion$Gene
del_mat = as.matrix(fusion[,2:9])

cell_fun4 = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "gray", fill = fill))
  if(del_mat[i, j] != 0) {
    if(del_mat[i, j] > 9) {
      grid.text(del_mat[i, j], x, y, gp = gpar(col = "white"))
    } else {
      grid.text(del_mat[i, j], x, y, gp = gpar(col = "black"))
    }
  }
}

col_fun4 = colorRamp2(c(0, 2, 5, 10, 20), c("white","#E6E6FA", "#CEA2FD", "#9370DB", "#9932CC"))

hm4 <- Heatmap(del_mat, 
               name = "#_of_LOE_Fusion",
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE,
               show_column_names = FALSE,
               column_names_side = "top",
               row_names_side = "left",
               column_names_gp = gpar(fontsize=13,fontface="bold"),
               row_names_gp = gpar(fontsize=10,fontface="bold"),
               cell_fun = cell_fun4,
               col = col_fun4,
               show_heatmap_legend = FALSE,
               width = unit(7, "cm"),
               height = unit(3, "cm")
)
hm4

###############################################################################
#sig
rownames(sig) = sig$Gene
sig_mat = as.matrix(sig[,2:9])
pdx = as.numeric(as.matrix(pdx_no[,2:9]))


cell_fun5 = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "gray", fill = fill))
  if(sig_mat[i, j] != 0) {
    if(sig_mat[i, j] > 9) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(col = "white"))
    } else {
      grid.text(sig_mat[i, j], x, y, gp = gpar(col = "black"))
    }
  }
}

col_fun5 = colorRamp2(c(0, 2, 5, 10, 20), c("#FFFFFF", "#D3D3D3", "#808080", "#A9A9A9", "#000000"))

col_anno <- columnAnnotation(Pct_model = anno_barplot(pdx,
                                                      bar_width = 0.3,
                                                      gp = gpar(fill= c("#F39B7F","#E64B35","#00A087","#8491B4",
                                                                        "#3C5488","#8e0152",
                                                                        "#7E6148","#B09C85")),
                                                      add_numbers = TRUE,
                                                      height = unit(2, "cm"),
                                                      numbers_rot = 0,
                                                      numbers_gp = gpar(fontsize = 12),
                                                      border = TRUE,
                                                      name=NULL
))

hm5 <- Heatmap(sig_mat, 
               name = "#_of_LOE_Sig",
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE,
               show_column_names = FALSE,
               column_names_side = "top",
               row_names_side = "left",
               column_names_gp = gpar(fontsize=13,fontface="bold"),
               row_names_gp = gpar(fontsize=10,fontface="bold"),
               cell_fun = cell_fun5,
               col = col_fun5,
               bottom_annotation = col_anno,
               show_heatmap_legend = FALSE,
               width = unit(7, "cm"),
               height = unit(1, "cm")
)
hm5

###############################################################################
ht_list = hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5

draw(ht_list, show_annotation_legend = FALSE)

