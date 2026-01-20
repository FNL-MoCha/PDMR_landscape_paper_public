# my script

suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
library("maftools")
library("ComplexHeatmap")
library(circlize)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(ggsci)


mat <- read.csv("../DATA/oncoplot_data_pdc.csv", check.names = F, stringsAsFactors = F)

# Function to set the size of the rectangle to be plotted as well as the color based on variat type
alter_fun = list(
  background = function(x, y, w, h) {grid.rect(x, y, w-unit(0.05, "mm"), h-unit(0.05, "mm"), gp = gpar(fill = "white", col = NA))},
  AMP        = function(x, y, w, h) {grid.rect(x, y, w-unit(0.05, "mm"), 0.6*h-unit(0.05, "mm"), gp = gpar(fill = "#de2d26", col = NA))},
  DEL        = function(x, y, w, h) {grid.rect(x, y, w-unit(0.05, "mm"), 0.6*h-unit(0.05, "mm"), gp = gpar(fill = "#3182bd", col = NA))},
  MUT        = function(x, y, w, h) {grid.rect(x, y, w-unit(0.05, "mm"), h*0.6,            gp = gpar(fill = "#008000", col = NA))},
  IND        = function(x, y, w, h) {grid.rect(x, y, w-unit(0.05, "mm"), h*0.3,            gp = gpar(fill = "#E69F00", col = NA))}
)

col = c(AMP = "#de2d26", DEL = "#3182bd", 
        MUT = "#008000", IND ="#E69F00")


#################################################################################################

#Make the top annotations
#need to sort the annotation file based on the matrix created by oncoplot
clinical_input <- read.csv("../DATA/Supp_tables_pdc_metadata.csv", check.names = F, stringsAsFactors = F)


mat2 <- mat[,2:dim(mat)[2]]
rownames(mat2) <- mat[,1]

idx <- which(colnames(mat2) %in% clinical_input$Tumor_Sample_Barcode)
mat2 <- mat2[,idx]

mat <- mat2


# idx <- which(clinical_input$)
clinical_input$Tumor_Sample_Barcode <- factor(clinical_input$Tumor_Sample_Barcode, levels=colnames(mat))
clinical_input <- clinical_input[order(clinical_input$Tumor_Sample_Barcode),]
clinical_input <- clinical_input[1:dim(mat)[2],]

for (i in 1:length(colnames(mat))){
  print(grep(colnames(mat)[i], clinical_input$Tumor_Sample_Barcode))
}

temp <- clinical_input$`Median TMB(mut/Mb)`
temp[temp>40] <- 40
clinical_input$`Median TMB(mut/Mb)` <- temp


##################
# start here

LOH <- clinical_input$`Median LOH Percent`
LOH2 <- clinical_input$`Median LOH Percent`
LOH[LOH2<=15] <- "0-15"
LOH[LOH2>15 & LOH2<=30] <- "16-30"
LOH[LOH2>30] <- ">30"

Aneuploidy <- clinical_input$Anueploidy
Aneuploidy2 <- clinical_input$Anueploidy
Aneuploidy[Aneuploidy2 <= 7] <- "0-7"
Aneuploidy[Aneuploidy2 > 7 & Aneuploidy2 <= 15] <- "8-15"
Aneuploidy[Aneuploidy2 > 15] <- ">15"



top_ann <-HeatmapAnnotation(
  show_annotation_name = F,
  barplot = anno_barplot(
    clinical_input$`Median TMB(mut/Mb)`,
    baseline=0,
    axis = TRUE,
    border=FALSE,
    ylim=c(0,40), # This is set to 80 max to show the variation in sample which are between 20 and 80, there are samples with TMB >200 and that suppresses the whole scale 
    gp=gpar(border =NA,fill="black",lty="blank")
  ),
  Signature=clinical_input$`Interpreted mutational Signature (V2)`,
  WGD=clinical_input$`Whole genome doubling (WGD)`,
  Aneuploidy=Aneuploidy,  ### changed here
  MSI=clinical_input$`MSI Status (model)`,
  LOH=LOH,    ### changed here
  ONCOTREE=clinical_input$ONCOTREE_CODE_in_FIGURES,
  col = list(
    Signature=c(
      "MMR"="#4C4C4C",
      "UV"="#E69F00",
      "TMZ"="#56B4E9",
      "BRCA1/2"="#009E73",
      "Smoking"="#F0E442",
      "POLE"="#0072B2",
      "AID/APOBEC"="#D55E00",
      "TMZ"="#DC0000FF",
      "Other"="#808180FF"
    ),
    WGD=c(
      "N"="#008B45FF",
      "Y"="#BB0021FF",
      "Y;N"="#631879FF"
    ),
    ONCOTREE=c(
      "COADREAD"="#E64B35FF",
      "HNSC"="#00A087",
      "NSCLC"="#3C5488",
      "BLCA"="#F39B7F",
      "MEL"="#8491B4",
      "PAAD"="#8e0152",
      "SARCOMA"="#7E6148",
      "Other"="#B09C85"
    ),
    MSI = c("MSI-H" =  "#EE0000FF",
            "MSI-S" = "light gray"),
    LOH = c("0-15" = "white",
            "16-30" = "#E8EAF6",
            ">30" = "#3F51B5"),
    Aneuploidy = c("0-7" ="white",
                   "8-15" = "#E0F2F1",
                   ">15" = "#4DB6AC")
  ),
  na_col = "white",
  annotation_legend_param = list(
    Signature = list(
      title = "Signature",
      at = c("MMR", "UV", "BRCA1/2", "Smoking", "POLE", "AID/APOBEC","TMZ","Other"),
      labels = c("MMR", "UV", "BRCA1/2", "Smoking", "POLE", "AID/APOBEC","TMZ","Other"),
      title_gp = gpar(fontsize = 10, 
                      fontface = "bold"), 
      labels_gp = gpar(fontsize = 8)
    ),
    MSI = list(
      title = "MSI",
      at = c("MSI-H", "MSI-S"),
      labels = c("MSI-H", "MSI-S"),
      title_gp = gpar(fontsize = 10, 
                      fontface = "bold"), 
      labels_gp = gpar(fontsize = 8)
    ),
    LOH=list(
      title="LOH",
      at=c("0-15","16-30",">30"),
      labels=c("0-15","16-30",">30"),
      title_gp = gpar(fontsize = 10, 
                      fontface = "bold"), 
      labels_gp = gpar(fontsize = 8)
    ),
    Aneuploidy=list(
      title="Aneuploidy",
      at=c("0-7","8-15",">15"),
      labels=c("0-7","8-15",">15"),
      title_gp = gpar(fontsize = 10, 
                      fontface = "bold"), 
      labels_gp = gpar(fontsize = 8)
    ),
    ONCOTREE=list(
      title="Tumor Histology",
      at=c("COADREAD","HNSC","NSCLC",
             "BLCA","MEL","PAAD",
             "SARCOMA","Other"),
      labels=c("COADREAD","HNSC","NSCLC",
               "BLCA","MEL","PAAD",
               "SARCOMA","Other"),
      title_gp = gpar(fontsize = 10, 
                      fontface = "bold"), 
      labels_gp = gpar(fontsize = 8)
    )
  ),
  annotation_height = unit(c(2.5, 0.4, 0.4,0.4,0.4, 0.4,0.4), "cm")
)



#################################################################################################

# bottom_anno<-HeatmapAnnotation(
# )
#################################################################################################
#oncoPrint using complex heatmap

ht=oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
             alter_fun = alter_fun, col = col, 
             alter_fun_is_vectorized = F,
             remove_empty_columns = F,
             # show_row_barplot = F,
             top_annotation=top_ann,
             # top_annotation_height = unit(2.5, "cm"),
             # bottom_annotation = bottom_anno,
             # bottom_annotation_height = unit(1, "cm"),
             show_pct = F,
             #pct_gp = gpar(fontsize = 2),
             #pct_digits = 0,
             #pct_side = "left",
             remove_empty_rows = TRUE,
             row_order=c(1:dim(mat)[1]),
             right_annotation = NULL,
             #column_title = "OncoPrint for PDMR\n oncoKB oncogenic variants",
             column_order=c(1:dim(mat)[2]),
             row_names_gp = gpar(fontsize = 6, fontface=2),
             row_gap = unit(0.01, "mm"),
             height = unit(12,"cm"),
             heatmap_legend_param = list(
               title = "Alterations", 
               at = c("AMP", "DEL", "MUT","IND"), 
               labels = c("Amplification", "Deep deletion", "Mutations","Indel"),
               nrow = 4
             )
)

pdf("Figure_oncoplot_PDC.pdf",width = 10,height = 10)
draw(ht, padding = unit(c(.5,2.5,3,.5), "cm"), heatmap_legend_side = "right",
     annotation_legend_side = "bottom")

decorate_annotation(
  "barplot", {
    grid.text("TMB\n(Mut/mb)", unit(-10, "mm"), just = "bottom", rot = 90,gp=gpar(fontsize=8))
  }
)
decorate_annotation(
  "MSI",{
    grid.text("MSI Status", unit(-10, "mm"), just = "bottom",gp=gpar(fontsize=8))
  }
)
decorate_annotation(
  "Signature",{
    grid.text("Mutational Signature", unit(-15, "mm"), just = "bottom",gp=gpar(fontsize=8))
  }
)
decorate_annotation(
  "WGD",{
    grid.text("WGD", unit(-10, "mm"), just = "bottom",gp=gpar(fontsize=8))
  }
)
decorate_annotation(
  "Aneuploidy",{
    grid.text("Aneuploidy", unit(-10, "mm"), just = "bottom",gp=gpar(fontsize=8))
  }
)
decorate_annotation(
  "LOH",{
    grid.text("% LOH", unit(-10, "mm"), just = "bottom",gp=gpar(fontsize=8))
  }
)
decorate_annotation(
  "ONCOTREE",{
    grid.text("Tumor Histology", unit(-10, "mm"), just = "bottom",gp=gpar(fontsize=8))
  }
)

dev.off()

