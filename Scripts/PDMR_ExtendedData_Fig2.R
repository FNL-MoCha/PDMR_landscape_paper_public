# PDMR Supp Figure 2 Script
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(dplyr)
library(readr)

#Read the PDMR Supp Table 1
df <- read.table("PDMR_Supp_Fig2_table1.txt", sep = "\t",header = TRUE)

#IMPACT/TCGA-(PDX, PDC, PDOrg)
df$label <- ifelse(df$Hugo %in% c("TP53", "APC", "KRAS", "PIK3CA", "ARID1A"), df$Hugo, NA)

df_pdx <- df[df$Group == "PDX",]
df_pdc <- df[df$Group == "PDC",]
df_pdorg <- df[df$Group == "PDOrg",]

ggplot(data = df_pdx, aes(x=PDMR, y=IMPACT,label=label)) +
  geom_point(size = 5,show.legend = F,color="red") +
  geom_text_repel(aes(label=label), box.padding = 0.5, point.padding = 0.5, 
                  nudge_y = 0.5, nudge_x = 0.5, color="black",size=5) +
  scale_x_continuous(limits=c(0,70), expand = c(0, 0),breaks = scales::pretty_breaks(n = 5)) + 
  scale_y_continuous(limits=c(0,50),expand = c(0, 0),breaks = scales::pretty_breaks(n = 6))+
  theme_classic()+
  xlab("\nFrequency in PDMR Cohort")+
  ylab("Frequency in MSK-IMPACT Cohort\n")+
  theme(
    axis.text.x = element_text(colour = "black", size = rel(1.5), hjust = 1),
    axis.text.y = element_text(colour = "black", size = rel(1.5), hjust = 1),
    legend.text=element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title = element_text(colour = "black", 
                              size = rel(1.5)
    )
  )


#BLCA/COAD/MEL/HNSC/NSCLC/PAAD/SARCOMA
supp_fig2 <- read.table("PDMR_Supp_Fig2_table2.txt", 
                        sep = "\t",header = TRUE)

input <- as.data.frame(supp_fig2)

df_blca <- input[input$Oncotree == "BLCA",]
df_coad <- input[input$Oncotree == "COAD",]
df_mel <- input[input$Oncotree == "MEL",]
df_hnsc <- input[input$Oncotree == "HNSC",]
df_nsclc <- input[input$Oncotree == "NSCLC",]
df_padd <- input[input$Oncotree == "PAAD",]
df_sarcnos <- input[input$Oncotree == "SARCOMA",]

ggplot(data = df_blca, aes(x=PDMR, y=IMPACT, label = Hugo_Symbol,color="red")) +
  geom_point(size = 3,show.legend = F) +
  scale_x_continuous(limits=c(0,100), expand = c(0, 0),breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(limits=c(0,100),expand = c(0, 0),breaks = scales::pretty_breaks(n = 6))+
  theme_classic()+
  stat_dens2d_filter(geom = "text_repel", keep.fraction = 0.1,col="black")+
  theme(
    axis.text.x = element_text(colour = "black", size = rel(1.5), hjust = 1),
    axis.text.y = element_text(colour = "black", size = rel(1.5), hjust = 1),
    legend.text=element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title = element_text(colour = "black", 
                              size = rel(1.5)
    )
  )


