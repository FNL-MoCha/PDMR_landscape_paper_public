# PDMR Supp Figure 9 Script

library(ggplot2)
library(dplyr)
library(ggtext)
library(patchwork)


supp_fig8 <- read.table("PDMR_Supp_Fig9_table1.txt",sep = "\t",header = TRUE)
df <- data.frame(supp_fig8)

hist_color = c("BLCA"="#F39B7F","COADREAD"="#E64B35","HNSC"="#00A087","MEL"="#8491B4",
               "NSCLC"="#3C5488","PAAD"="#8e0152","SARCOMA"="#7E6148","Other"="#B09C85")
hist_color <- data.frame(Hist = names(hist_color), color = hist_color)

df <- merge(df,hist_color, by ="Hist", all.x = TRUE)
df$ID <- paste0("<span style=\"color: ", df$color, "\">", df$modID, "</span>")


df$modID <-factor(df$modID, levels = c("261386-189-R", "281475-159-R", "324938-238-R", "343268-274-T", "639262-121-R", 
                                       "665939-344-R", "746718-042-R", "772611-094-R", "116655-072-R", "147771-066-R1", 
                                       "167148-078-R", "175126-011-R", "191243-178-R", "255893-291-R", "259778-044-R", 
                                       "273589-319-R", "282377-053-R", "296347-364-R", "328469-098-R", "341922-053-R", 
                                       "413561-133-T", "431354-103-R", "435261-313-R", "451658-271-R", "483684-250-R", 
                                       "519858-162-T", "522493-142-R", "524269-139-T", "736525-025-R", "762968-020-R", 
                                       "764851-200-R", "784581-064-R", "784911-089-T", "825966-067-R", "869496-051-R", 
                                       "884544-143-R", "931267-113-T", "944381-210-T", "947725-317-R", "947758-054-R", 
                                       "974816-308-R", "217283-344-R", "245127-232-R", "266295-075-R", "323392-091-R", 
                                       "328373-195-R", "427551-204-R", "582836-169-R", "668155-338-R", "735871-273-R",
                                       "874868-142-R", "891969-043-R", "245324-029-R", "251568-266-R", "299254-011-R", 
                                       "633993-097-R", "863532-251-R", "992656-260-R", "436779-168-T", "692585-246-R", 
                                       "941728-121-R", "952719-076-R", "997726-040-R"))
df <- df[order(df$modID), ]

df$Hist <-factor(df$Hist, levels = c("Other","SARCOMA","PAAD","NSCLC","MEL","HNSC","COADREAD","BLCA"))
df <- df[order(df$Hist), ]
df$ID <- factor(df$ID,levels = unique(df$ID))


p <- ggplot(df, aes(x=ID, y=LOH)) +
  geom_vline(xintercept = unique(df$ID), color="lightgray", linetype="solid") +
  coord_flip() +
  geom_point(data = df %>% filter(sam_type != "Originator"), aes(color=sam_type), size=3, alpha=0.8) +
  geom_point(data = df %>% filter(sam_type == "Originator"), aes(color=sam_type), size=3, alpha=0.8) +
  scale_color_manual(values=c("PDX"="black", "Originator"="red","PDC"="blue", "PDOrg"="green")) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y = ggtext::element_markdown(size = 12,margin=margin(t=5, r=10, b=5, l=5)),
        legend.position="none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p


supp_fig8_o <- read.table("PDMR_Supp_Fig9_table2.txt",sep = "\t",header = TRUE)
df_o <- data.frame(supp_fig8_o)

df_o <- merge(df_o,hist_color, by ="Hist", all.x = TRUE)
df_o$ID <- paste0("<span style=\"color: ", df_o$color, "\">", df_o$modID, "</span>")


df_o$modID <-factor(df_o$modID, levels = c("158883-120-T", "194179-226-R", "197587-005-T", "292632-238-R", "377594-074-R", 
                                           "449892-231-R", "467112-017-R", "471365-239-R", "497265-261-R", "584427-346-R", 
                                           "594176-295-R", "691877-024-R", "713595-302-R", "761936-265-R", "767577-098-T", 
                                           "821394-179-R", "823721-103-R", "918122-036-R", "949626-316-R", "988467-305-R", 
                                           "165739-295-R", "485368-065-R3", "485368-065-R4", "496974-208-R", "521955-158-R3", 
                                           "521955-158-R4", "521955-158-R5", "521955-158-R6", "521955-158-R7", "777334-354-R3", 
                                           "117519-064-T", "144983-106-R", "146199-324-R", "174316-266-R", "193523-008-R", 
                                           "193832-021-R", "262622-085-R", "275155-148-R", "278458-045-R", "313798-341-R", 
                                           "326966-086-R", "331888-031-R", "346799-050-R", "388244-064-R", "415267-285-R", 
                                           "415371-026-R", "417952-058-R", "419622-098-T", "459288-142-R", "464991-137-R", 
                                           "466732-252-T", "544552-058-R", "562452-108-R", "562715-036-R", "597326-320-R", 
                                           "615535-235-R", "627122-101-R", "629538-142-R", "636577-100-R", "648538-100-R", 
                                           "682317-045-R", "748385-122-R", "813916-060-R", "824345-141-R", "862989-056-R", 
                                           "877282-077-R", "885512-296-R", "933738-175-T", "989133-093-R", "995276-142-R"))
df_o <- df_o[order(df_o$modID), ]

df_o$Hist <-factor(df_o$Hist, levels = c("Other","SARCOMA","PAAD","NSCLC","MEL","HNSC","COADREAD","BLCA"))
df_o <- df_o[order(df_o$Hist), ]
df_o$ID <- factor(df_o$ID,levels = unique(df_o$ID))


k <- ggplot(df_o, aes(x=ID, y=LOH)) +
  geom_vline(xintercept = unique(df_o$ID), color="lightgray", linetype="solid") +
  coord_flip() +
  geom_point(data = df_o %>% filter(sam_type != "Originator"), aes(color=sam_type), size=3, alpha=0.8) +
  geom_point(data = df_o %>% filter(sam_type == "Originator"), aes(color=sam_type), size=3, alpha=0.8) +
  scale_color_manual(values=c("PDX"="black", "Originator"="red","PDC"="blue", "PDOrg"="green")) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60,80)) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y = ggtext::element_markdown(size = 12,margin=margin(t=5, r=10, b=5, l=5)),
        legend.position="none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
k



combined_plot <- p + k + plot_layout(ncol = 2)

print(combined_plot)

