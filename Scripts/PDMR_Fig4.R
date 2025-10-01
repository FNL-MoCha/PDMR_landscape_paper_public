input <- read.table("/Volumes/MoCha/processedDATA/MATH.txt", header = T,sep="\t")
input2 <-input[grep("ORIGINATOR", input$Sample),]

require(ggplot2)
mysep <- function(x, ...) 
  format(x, big.mark = ' ', trim = TRUE, scientific = FALSE, ...) 

ggplot() +
  geom_boxplot(data = input,aes(x=Model, y=MATH)) +
  geom_point(data = input2,aes(x=Model, y=MATH,color = 'red'))+
  theme_classic()+
  scale_y_continuous(label= mysep,limits = c(0,150))+
  scale_x_discrete(limits= unique(input$Model))+
  theme(axis.text.x = element_text(colour = "black",angle = 90),
    #legend.text=element_text(size=10),
    legend.title = element_blank()
    #axis.title = element_text(colour = "black", size = rel(1.5))
  )

##############################################################

require(ggplot2)
mysep <- function(x, ...) 
  format(x, big.mark = ' ', trim = TRUE, scientific = FALSE, ...) 

library(ggh4x)
#> Loading required package: ggplot2


input <- read.table("ppa_pdx.txt", header = T,sep="\t")
input2<-read.table("ppa2_pdx.txt", header = T,sep = "\t")
#input2 <-input[grep("ORIGINATOR", input$PDX),]


input2$PPA_refOriginator <- input2$Common/(input2$Common+input2$Originator.Only)

input2$Model <- factor(input2$Model, levels = input2$Model[order(input2$PPA_refOriginator)])
input$Model <- factor(input$Model, levels=input2$Model)

pdf("Fig4A_PPA_Originator.pdf",width = 10,height = 5)
p<-ggplot() +
  geom_point(data = input2,aes(x=Model, y=PPA_refOriginator,color = 'red'))+
  geom_boxplot(data = input,aes(x=Model, y=PPA_refOriginator)) +
  theme_classic()+
  scale_y_continuous(label= mysep,limits = c(0,1))+
  # scale_x_discrete(limits= unique(input$Model))+
  theme(axis.text.x = element_text(colour = "black",angle = 90, size = 4),
        #legend.text=element_text(size=10),
        legend.title = element_blank()
        #axis.title = element_text(colour = "black", size = rel(1.5))
  )
print(p)
dev.off()




##############################################
# add histology information

require(ggplot2)
mysep <- function(x, ...) 
  format(x, big.mark = ' ', trim = TRUE, scientific = FALSE, ...) 

library(ggh4x)
#> Loading required package: ggplot2


input <- read.table("ppa_pdx.txt", header = T,sep="\t")
input2<-read.table("ppa2_pdx.txt", header = T,sep = "\t")
#input2 <-input[grep("ORIGINATOR", input$PDX),]

#t <- input$PDX.PDC.Organoid.Sample
#nt <- gsub("~","-", t)
#st2 <- read.csv("Supp_tables_v27_st2.csv", check.names = F, stringsAsFactors = F)


meta.data <- read.csv("Supp_tables_pdx_metadata.csv", check.names = F, stringsAsFactors = F)

t <- input$Model
nt <- gsub("~","-", t)
input$Model <- nt

t <- input2$Model
nt <- gsub("~","-", t)
input2$Model <- nt


input$oncotree <- meta.data$ONCOTREE_CODE_in_FIGURES[match(input$Model,meta.data$Tumor_Sample_Barcode)]
input2$oncotree <- meta.data$ONCOTREE_CODE_in_FIGURES[match(input2$Model,meta.data$Tumor_Sample_Barcode)]

input$oncotree <- factor(input$oncotree, levels = c("COADREAD","HNSC","NSCLC","BLCA","MEL","PAAD","SARCOMA","Other"))
input2$oncotree <- factor(input2$oncotree, levels = c("COADREAD","HNSC","NSCLC","BLCA","MEL","PAAD","SARCOMA","Other"))

require(ggplot2)
mysep <- function(x, ...) 
  format(x, big.mark = ' ', trim = TRUE, scientific = FALSE, ...) 

input2$PPA_refOriginator <- input2$Common/(input2$Common+input2$Originator.Only)

input2$Model <- factor(input2$Model, levels = input2$Model[order(input2$PPA_refOriginator)])
input$Model <- factor(input$Model, levels=input2$Model)

p<-ggplot() +
  geom_point(data = input2,aes(x=Model, y=PPA_refOriginator,color = 'red'),show.legend = FALSE)+
  geom_boxplot(data = input,aes(x=Model, y=PPA_refOriginator), outlier.size = 0.5) +
  theme_classic()+
  scale_y_continuous(label= mysep)+
  # ,limits = c(0,1)
  # scale_x_discrete(limits= unique(input$Model))+
  theme(# axis.text.x = element_text(colour = "black",angle = 90, size = 4),
        #legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank()
        #axis.title = element_text(colour = "black", size = rel(1.5))
  )



oncotree.color <- c("#E64B35FF","#00A087","#3C5488","#F39B7F","#8491B4","#8e0152","#7E6148","#B09C85")
# Only colour strips in x-direction
strip <- strip_themed(text_x=element_text(color = "white"),background_x = elem_list_rect(fill = oncotree.color))

# p <- p+facet_wrap(~oncotree,ncol=2, scales = "free")
p2 <- p+ facet_wrap2(~oncotree, ncol=2, scales="free", strip = strip)
print(p2)


pdf("Figure_PPA_Originator_PDX.pdf",width = 7.5,height = 5)
print(p2)
dev.off()


###########################
### plot for ppa_pdc and ppa_pdorg
ppa.pdc <- read.table("ppa_pdc.txt", header = T,sep="\t")

ppa.pdc$PPA_refOriginator <- ppa.pdc$Common/(ppa.pdc$Common+ppa.pdc$Originator.Only)
ppa.pdc$source <- "PDC"

ppa.pdorg <- read.table("ppa_pdorg.txt", header = T,sep="\t")
ppa.pdorg$PPA_refOriginator <- ppa.pdorg$Common/(ppa.pdorg$Common+ppa.pdorg$Originator.Only)
ppa.pdorg$source <- "PDOrg"

temp <- rbind(ppa.pdc, ppa.pdorg)
temp$Model <- factor(temp$Model, levels = unique(temp$Model[order(temp$PPA_refOriginator)]))

pdf("Figure_PPA_Originator_PDCOrg.pdf",width = 7,height = 5)

p<-ggplot() +
  geom_point(data = temp,aes(x=Model, y=PPA_refOriginator, shape=source,color=source))+
  scale_shape_manual(values=c(16, 2))+
  scale_color_manual(values=c('#E64B35', 'black'))+ 
  
  # geom_point(data = ppa.pdorg,aes(x=Model, y=PPA_refOriginator,color = 'black')) +
  # geom_boxplot(data = input,aes(x=Model, y=PPA_refOriginator)) +
  theme_classic()+
  scale_y_continuous(label= mysep,limits = c(0.5,1))+
  # scale_x_discrete(limits= unique(input$Model))+
  theme(#axis.text.x = element_text(colour = "black",angle = 90, size = 6),
    #legend.text=element_text(size=10),
    axis.text.x = element_blank(),axis.ticks.x=element_blank(),
    legend.title = element_blank()
    #axis.title = element_text(colour = "black", size = rel(1.5))
  )
print(p)

dev.off()

