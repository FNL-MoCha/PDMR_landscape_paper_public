# PDMR Figure 1 A Script

library(dplyr)
library(ComplexUpset)
library(ggplot2)

#Read the PDMR Supp Table 1
data <- read.table("PDMR_Supp_Table1.txt", sep = "\t", header = TRUE)

new_data <- data %>%
  select(
    Model             = `Model.Set.ID`,
    Oncotree          = `Tumor.histology.code.for.manuscript.figures`,
    PDC               = `Number.of.PDC`,
    PDOrg             = `Number.of.PDOrg`,
    PDX               = `Number.of.PDX.Samples`
  ) %>%
  mutate(
    PDC   = if_else(PDC != 0, 1, 0),
    PDOrg = if_else(PDOrg != 0, 1, 0),
    PDX   = if_else(PDX != 0, 1, 0)
  ) %>%
  # Keep only the columns we need in the final table
  select(
    Model,
    Oncotree,
    PDC,
    PDOrg,
    PDX
  )


df = as.data.frame(new_data)
pdx_source = colnames(df)[3:5]

upset(
  df,
  pdx_source,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill=factor(Oncotree, 
                              levels=c("BLCA","MEL","NSCLC","HNSC",
                                       "PAAD","SARCOMA","COADREAD","Other")))) +
      scale_fill_manual(
        values = c("COADREAD"="#E64B35","HNSC"="#00A087","NSCLC"="#3C5488",
                   "BLCA"="#F39B7F","MEL"="#8491B4","PAAD"="#8e0152",
                   "SARCOMA"="#7E6148","Other"="#B09C85"),
        guide="none"
      ) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            strip.clip = "off"
      )
  ),
  width_ratio=0.1,
  set_sizes=FALSE
)




