
require(ConsensusClusterPlus)
require(gplots)
require(RColorBrewer)
require(mixtools)
require(ggplot2)
require(reshape2)
library('ggsci')


# figure 3B
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

####################################################
# combined figure for all passages, pdc and pdorg
op <- read.csv("CNV_stability_paper_noPforfigure.csv", check.names = F, stringsAsFactors = F)
op.ori <- op

temp <- op$PDX.Passage.or.PDC.or.PDOrg.or.Patient.Orginator.specimen

idx.pdc <- grep("PDOrg", op$PDX.Passage.or.PDC.or.PDOrg.or.Patient.Orginator.specimen)
idx.pdorg <- grep("PDC", op$PDX.Passage.or.PDC.or.PDOrg.or.Patient.Orginator.specimen)

temp[idx.pdc] <- '-100'
temp[idx.pdorg] <- '-200'

passage.dif <- as.integer(temp) + 1
passage.dif[passage.dif >= 4] <- 4
passage.dif[idx.pdc] <- 5
passage.dif[idx.pdorg] <- 6

maf.combind2 <- data.frame(passage.diff = passage.dif, value = as.numeric(op$CNVfraction.after)*100)
maf.combind2$passage.diff <- as.factor(maf.combind2$passage.diff)

p <- ggplot(maf.combind2, aes(x=passage.diff, y=value)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.2),col="black") +
  # geom_boxplot(notch=F, outlier.colour= NULL, outlier.shape= NA, outlier.size=4,alpha=0.8) +
  # scale_fill_brewer(palette="Greens") +
  # scale_fill_npg() +
  labs(x="Passage# compared to Originator",y="Change in PDX CNA fraction (%)") +
  scale_x_discrete(labels=c("P0", "P1", "P2",">=P3", "PDC", "PDOrg")) +
  theme_classic() + 
  theme(legend.position = "none") +
  # stat_summary(fun.data=data_summary,col='red')
  stat_summary(fun=median, geom="point", size=2, color="red")

pdf("Figure_CNV_stability_combined.pdf", width = 8,height = 4.5)
print(p)
dev.off()

