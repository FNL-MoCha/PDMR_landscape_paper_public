# PDMR Figure 1 B Script

library(dplyr)
library(ggplot2)
library(ggrepel)

# A helper function to generate a pie chart
make_pie <- function(data, column_name) {
  # Summarize the data
  df_count <- data %>%
    count(.data[[column_name]]) %>% 
    mutate(
      pct = n / sum(n),
      label_factor = as.factor(.data[[column_name]])
    )
  
  # Create the pie chart
  p <- ggplot(df_count, aes(x = "", y = pct, fill = label_factor)) +
    geom_bar(stat = "identity", color = "white") +
    geom_text(
      aes(label = label_factor), 
      position = position_stack(vjust = 0.5), 
      size = 6, color = "white"
    ) +
    coord_polar("y") +
    scale_fill_manual(values=c("#E64B35", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#8e0152", "#7E6148", "#B09C85","red","blue")) +
    theme_void() +
    theme(legend.position="none")
  
  return(p)
}

#Read the PDMR Supp Table 1 for five pie charts
  #Rare/Common
  #Inferred ancestry
  #Metastatic Disease
  #Number of specimens/model
  #age at diagnosis
  #germline data available
data <- read.table("PDMR_Supp_Table1.txt", sep = "\t",header=TRUE)

new_data <- data %>%
  select(
    Model         = `Model.Set.ID`,
    common_rare   = `Rare..Recalcitrant..Common.cancer`,
    ancestry      = `Inferred.Ancestry`,
    metastatic    = `Has.known.Metastatic.Disease`,
    PDX_sam       = `Number.of.PDX.Samples`,
    age           = `Age.at.Diagnosis`,
    germline      = `Germline.Available`
  )

df <-data.frame(new_data)

# List of columns for which we want pie charts
pie_vars <- c("common_rare", "ancestry", "metastatic", 
              "PDX_sam", "age", "germline")

# Generate and print (or save) each pie chart
for (col in pie_vars) {
  pie_plot <- make_pie(df, col)
  print(pie_plot)
}


#Read the PDMR Supp Table 3 for two additional pie charts
  #Culture of origin for PDC
  #Culture of origin for PDOrg

#Read the PDMR Supp Table 3
data1 <- read.table("PDMR_Supp_Table3.txt", sep = "\t",header = TRUE)

new_data1 <- data1 %>%
  select(
    source         = `PDC.PDOrg`,
    origin         = `Culture.Origin`
  )

df1 <-data.frame(new_data1)


pdc_table   <- subset(df1, source == "PDC")
pdorg_table <- subset(df1, source == "PDOrg")

pie_pdc   <- make_pie(pdc_table,   "origin")
pie_pdorg <- make_pie(pdorg_table, "origin")

pie_pdc
pie_pdorg
