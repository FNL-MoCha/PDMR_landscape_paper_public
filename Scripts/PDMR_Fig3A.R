# PDMR Figure 3 A Script

library(ggplot2)
library(dplyr)
library(ggtext)
library(patchwork)

fig3a_c <- read.table("PDMR_Fig3A_table.txt",sep = "\t",header = TRUE)

df_c <-data.frame(fig3a_c)

hist_color = c("BLCA"="#F39B7F","COADREAD"="#E64B35","HNSC"="#00A087","MEL"="#8491B4",
               "NSCLC"="#3C5488","PAAD"="#8e0152","Other"="#B09C85")
hist_color <- data.frame(Hist = names(hist_color), color = hist_color)

df_c <- merge(df_c,hist_color, by ="Hist", all.x = TRUE)
df_c$ID <- paste0("<span style=\"color: ", df_c$color, "\">", df_c$ID, "</span>")

df_c$ModID <-factor(df_c$ModID, levels = c("116655-072-R", "144983-106-R", "165739-295-R", "167148-078-R", 
                                           "175126-011-R", "217283-344-R", "245127-232-R", "251568-266-R", 
                                           "255893-291-R", "261386-189-R", "278458-045-R", "282377-053-R", 
                                           "296347-364-R", "313798-341-R", "328469-098-R", "341922-053-R", 
                                           "346799-050-R", "413561-133-T", "415267-285-R", "427551-204-R", 
                                           "435261-313-R", "451658-271-R", "485368-065-R3", "485368-065-R4", 
                                           "496974-208-R", "519858-162-T", "521955-158-R3", "521955-158-R4", 
                                           "521955-158-R5", "521955-158-R6", "521955-158-R7", "522493-142-R", 
                                           "524269-139-T", "544552-058-R", "639262-121-R", "692585-246-R", 
                                           "746718-042-R", "762968-020-R", "764851-200-R", "772611-094-R", 
                                           "777334-354-R3", "813916-060-R", "825966-067-R", "931267-113-T", 
                                           "947758-054-R", "992656-260-R", "997726-040-R"))
df_c <- df_c[order(df_c$ModID), ]


df_c$Var <-factor(df_c$Var, levels = c("p.E545K","p.G12D","p.G12V","p.H1047R",
                                       "p.R175H","p.R248Q", "p.R248W","p.R273H","p.V600E"))
df_c <- df_c[order(df_c$Var), ]
df_c$Gene <-factor(df_c$Gene, levels = c("BRAF", "KRAS", "PIK3CA", "TP53"))
df_c <- df_c[order(df_c$Gene), ]

df_c$Hist <-factor(df_c$Hist, levels = c("Other","PAAD","NSCLC","MEL","HNSC","COADREAD","BLCA"))
df_c <- df_c[order(df_c$Hist), ]
df_c$ID <- factor(df_c$ID,levels = unique(df_c$ID))


k <- ggplot(df_c, aes(x=ID, y=VAF)) +
  geom_vline(xintercept = unique(df_c$ID), color="lightgray", linetype="solid") +
  coord_flip() +
  geom_point(data = df_c %>% filter(Type != "Originator"), aes(color=Type), size=3, alpha=0.8) +
  geom_point(data = df_c %>% filter(Type == "Originator"), aes(color=Type), size=3, alpha=0.8) +
  scale_color_manual(values=c("PDX"="black", "Originator"="red","PDC"="blue", "PDOrg"="green")) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.5, 0.75, 1.00)) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y = ggtext::element_markdown(size = 12,margin=margin(t=5, r=10, b=5, l=5)),
        legend.position="none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

k

