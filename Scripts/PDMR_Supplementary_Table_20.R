#R Scripts to calculate response
#These require exports to match file and field names used in the scripts.  Similar calculation metrics can be found using the web-based R/Shiny Tumor Volume Suite developed for PDXNet at https://tumor-volume.jax.org (Meric-Bernstam et al., 2024)

#EFSx4
### Calculate event free survival where an event equals tumor quadrupling from staging (EFSx4).  Right censored. Authors: E. Polley and L. Rubenstein

library(survival) # survival a base package, so should always be available
if(!require(ggplot2)) stop("Need to install ggplot2 package. Open R and run command `install.packages(\"ggplot2\")'")
library(ggplot2)
if(!require(RColorBrewer)) stop("Need to install RColorBrewer package. Open R and run command `install.packages(\"RColorBrewer\")'")
#library(xtable)# only required to export LaTeX tables

## Usage: Rscript PDX_analysis_script.R "Data2/" "zabj2-2" "BL0382"

#args <- commandArgs(trailingOnly = TRUE)
#if(length(args) != 3) stop("script requires 3 arguments: data directory, experiment number, and model name")
#DataDir <- args[1]  # e.g. Data/
#ExpNumber <- args[2]  # e.g. zdgj2-1
#ModelName <- args[3]  # e.g. 12346-123-R

#install.packages("ggplot2")

# CSV files used: therapy.csv, tumor.csv, groups.csv, animal.csv, and wts.csv
# Set working directory, this need to be changed depending on location of working directory
# Enter ExpNumber and ModelName below before running script. Example data ExpID = ZBTJ2-2, ModelID = 156681-154-R
DataDir <- "C:\\Users\\username\\Desktop\\R_Scripts\\Data\\datafiles\\"
ExpNumber <-"EXPID"
ModelName <- "MODELID"

#print(pwd())
setwd(DataDir) # better than doing nothing
print(getwd()) # fixed
print(list.files(path = DataDir, pattern = ExpNumber))

# 5 csv files used to process data, exported from StudyLog.  If naming or header structure of files changes, then may need to alter the read.csv() commands below
# Tumor measurements
tumor_wts <- read.csv(file = file.path(DataDir, paste0(ExpNumber, " tumor.csv")), header = TRUE, stringsAsFactors = FALSE)

# therapy
therapy <- read.csv(file = file.path(DataDir, paste0(ExpNumber, " therapy.csv")), header = TRUE, stringsAsFactors = FALSE)

# groups
therapy_groups <- read.csv(file = file.path(DataDir, paste0(ExpNumber, " groups.csv")), header = TRUE, stringsAsFactors = FALSE)

# animals
animals <- read.csv(file = file.path(DataDir, paste0(ExpNumber, " animal.csv")), header = TRUE, stringsAsFactors = FALSE)

# animal wts
animal_wts <- read.csv(file = file.path(DataDir, paste0(ExpNumber, " an wts.csv")), header = TRUE, stringsAsFactors = FALSE)


# create a group label and merge, this merges combination treatments
GROUP_CODE <- data.frame(GROUP_NBR = unique(therapy$GROUP_NBR), NSC = NA)
for(ii in seq_along(GROUP_CODE$GROUP_NBR)) {
  GROUP_CODE$NSC[ii] <- paste0(therapy$NSC[which(therapy$GROUP_NBR == GROUP_CODE$GROUP_NBR[ii])], collapse = " and ")
}

# Agent_Names can be used to decode NSC numbers, StudyLog tracks therapy by NSC number, below will return data with drug name
# NSC numbers for control material predefined to start with “99999”
Agent_Names <- data.frame(NSC = c("999999", "241240", "125973", "19893", "616348", "718781", "750927", "761431"), AgentName = c("Control", "Carboplatin", "Paclitaxel", "5-FU", "Irinotecan", "Erlotinib", "Gemcitabine", "Vemurafenib"), stringsAsFactors = FALSE)
GROUP_CODE <- merge(GROUP_CODE, Agent_Names, all.x = TRUE)

GROUP_CODE$AgentName[is.na(GROUP_CODE$AgentName)] <- GROUP_CODE$NSC[is.na(GROUP_CODE$AgentName)]

# compute time to event which is tumor volume quadrupling (2 doublings)
# function for a single mouse
time_to_quad <- function(time, volume) {
  t0 <- time[1]
  vol1 <- volume[1]
  Quad <- 4*vol1
  Event <- ifelse(max(volume, na.rm = TRUE) > Quad, 1, 0)
  if(Event == 1) {
    # find times that bracket event
    time_upper <- min(time[volume > Quad], na.rm = TRUE)
    time_lower <- time[which(time == time_upper) - 1] # step back one time point, this will give first crossing
    volume_lower <- volume[time == time_lower]
    volume_upper <- volume[time == time_upper]
    # interpolate with exponential growth
    time_event <- time_lower + (time_upper - time_lower)*log(Quad/volume_lower)/log(volume_upper/volume_lower)
  } else if(Event == 0) {
    time_event <- max(time) # last observation
  }
  out <- data.frame(time = time_event, event = Event)
  return(out)
}


# create dataset for event time
Event_df <- unique(tumor_wts[, c("GROUP_NBR", "ANIMAL_NBR")])
Event_df$time <- Event_df$event <- NA 

# Loop over all mice and compute time to tumor volume quadrupling:
# time is the quadrupling time if event = 1 (observed quadrupling), or time is last observation day when event = 0 (right-censored time to quadrupling)
for(ii in seq(nrow(Event_df))) {
  tmp <- time_to_quad(time = tumor_wts$OBS_DAY[tumor_wts$GROUP_NBR == Event_df$GROUP_NBR[ii] & tumor_wts$ANIMAL_NBR == Event_df$ANIMAL_NBR[ii]], volume = tumor_wts$TUMOR_WT[tumor_wts$GROUP_NBR == Event_df$GROUP_NBR[ii] & tumor_wts$ANIMAL_NBR == Event_df$ANIMAL_NBR[ii]])
  Event_df[ii, "time"] <- tmp$time
  Event_df[ii, "event"] <- tmp$event
}

# update events based on death code. Toxicity related deaths are censored out before export of data from StudyLog.
Event_df <- merge(Event_df, animals[, c("GROUP_NBR", "ANIMAL_NBR", "DEATH_CODE", "ANIMAL_COMMENT", "DEATH_DAY")])
print(table(Event_df$DEATH_CODE))

Event_df$event[Event_df$DEATH_CODE %in% c("TRS", "TRD", "DRS", "DRD") & Event_df$event == 0] <- 1 # make event if mouse sacrified for tumor or drug related

Event_df <- merge(Event_df, GROUP_CODE) # add agent names

write.csv(Event_df, paste0("EventTable_", ExpNumber, ".csv"))
# This will write out a file where each row is a mouse and report the time to event (doubling or quadrupling) with the event indicator ("event" column) equal to 1 if it is an observed event, or 0 if it is right censored.

# plot histograms
g <- ggplot(Event_df, aes(x = time)) + geom_histogram(binwidth = 5) + facet_wrap(~AgentName, ncol = 2) + theme_bw() + xlab("Event Time")
ggsave(paste0("Hist_event_times_", ExpNumber, "_", ModelName, ".pdf"), width = 8, height = 11)

# output a table with the event times
write.csv(Event_df, file = paste0("TimeToEvents_", ExpNumber, ".csv"), row.names = FALSE)

# now fit Kaplan Meier
survFit <- survfit(Surv(time, event)~AgentName, data = Event_df, conf.type = "log") # try log-log?
survFit # gives estimates of the median time to event

# write table with median times
tempTable <- summary(survFit)$table[, c("n.start", "events", "median", "0.95LCL", "0.95UCL")]
rownames(tempTable) <- sub("AgentName=", "", rownames(tempTable))
write.csv(tempTable, paste0("Median_TimeToEvents_", ExpNumber, ".csv"), row.names = TRUE)
# if you prefer LaTeX table, use code below
# print(xtable(tempTable, digits = c(0, 0, 0, 1, 1, 1), caption = "Median event free survival times"), type = "latex", file = paste0("Median_EFS_", ExpNumber, "_", ModelName, ".tex"), booktabs = TRUE)

# group colors
# display.brewer.all()
mypalette <- brewer.pal(length(unique(Event_df$AgentName)), "Paired") # can pick other color schemes if prefer, http://www.datavis.ca/sasmac/brewerpal.html "paired" gives 12 groups

# generates PDF with the event free survival kaplan meier with legend
pdf(paste0("EventFree_Kaplan_Meier_", ExpNumber, "_", ModelName, ".pdf"), width = 14, height = 10)
plot(survFit, col = mypalette, xlab = "Days", ylab = "Event-Free Survival", lwd = 2)
legend(1, .3, sub("AgentName=", "", names(survFit$strata)), col = mypalette, lwd = 2)
dev.off()

# get the Log Rank Test results for all pairwise comparisons, need to put together in a table
OUT <- matrix(NA, nrow = 0, ncol = 3)
colnames(OUT) <- c("Compare", "Chi-square", "p-value")

for(ii in 1:(length(unique(Event_df$GROUP_NBR)) - 1)) {
  for(jj in (ii+1):length(unique(Event_df$GROUP_NBR))) {
    foo <- survdiff(Surv(time, event)~AgentName, data = Event_df[Event_df$GROUP_NBR %in% c(ii, jj), ])
    tempOUT <- matrix(NA, ncol = 3)
    tempOUT[, 1] <- paste0(names(foo$n), collapse = " vs. ")
    tempOUT[, 2] <- round(foo$chisq, 2)
    tempOUT[, 3] <- signif(1 - pchisq(foo$chisq, df = 1), 3)
    OUT <- rbind(OUT, tempOUT)
  }
}

OUT.df <- as.data.frame(OUT, stringsAsFactors = FALSE)
OUT.df[, 2] <- as.numeric(OUT.df[, 2])
OUT.df[, 3] <- as.numeric(OUT.df[, 3])

# remove the "AgentName="
OUT.df$Compare <- gsub("AgentName=", "", OUT.df$Compare)

write.csv(OUT.df, paste0("LogRankTest_", ExpNumber, ".csv"), row.names = FALSE)
# LaTeX code:
# print(xtable(OUT.df, caption = "Log Rank Test for event free survival", align = c("l", "l", "l", "l"), digits = c(1, 1, 2, 4)), type = "latex", file = paste0("LogRank_", ExpNumber, "_", ModelName, ".tex"), booktabs = TRUE)


## Ratio of Medians
OUT <- matrix(NA, nrow = 0, ncol = 6)
OUT <- as.data.frame(OUT)


for(ii in 1:(length(unique(Event_df$GROUP_NBR)) - 1)) {
  for(jj in (ii+1):length(unique(Event_df$GROUP_NBR))) {
    cat(ii, jj, "\n")
    foo <- survfit(Surv(time, event)~AgentName, data = Event_df[Event_df$GROUP_NBR %in% c(ii, jj), ], error = "greenwood") # greenwood variance estimate
    tempOUT <- matrix(NA, ncol = 6)
    tempOUT[, 1] <- paste0(rownames(summary(foo)$table), collapse = " vs. ")
    tempOUT[, 2] <- summary(foo)$table[, "median"][[1]]
    tempOUT[, 3] <- summary(foo)$table[, "median"][[2]]
    tempOUT[, 4] <- summary(foo)$table[, "median"][[1]]/summary(foo)$table[, "median"][[2]]
    
    # CI
    # gamma = median1/median2
    gammaRange <- seq(from = 0.125, to = 8, length.out = 500)
    TimePoints <- sort(summary(foo)$time)
    outGamma <- data.frame(gamma = gammaRange, MinW = NA, LessThanChi2 = NA)
    for(gg in seq_along(gammaRange)){
      
      W_gg_tt <- rep(NA, length(TimePoints))
      for(tt in seq_along(TimePoints)){
        # If ALL times are beyond the range of dead time in any group, an error msg showed up:
        #   'Error in array(xx, dim = dd) : negative length vectors are not allowed'
        foo_gg_tt <- try(summary(foo, times = c(TimePoints[tt], gammaRange[gg]*TimePoints[tt])), silent = TRUE) 
        if (inherits(foo_gg_tt, "try-error")) next()
        if(any(foo_gg_tt$std.err == 0) | is.na(any(foo_gg_tt$std.err == 0))) next()
        if(length(foo_gg_tt$surv) < 4) next()
        W_gg_tt[tt] <- (foo_gg_tt$surv[1] - 0.5)^2/foo_gg_tt$std.err[1] + (foo_gg_tt$surv[4] - 0.5)^2/foo_gg_tt$std.err[4]
      }
      minW <- min(W_gg_tt, na.rm = TRUE)
      outGamma[gg, 2] <- minW
      outGamma[gg, 3] <- minW < qchisq(0.95, df = 1) # alpha value 0.05 (95% CI)
    }
    tempOUT[, 5] <- min(outGamma[outGamma[, 3], "gamma"])
    tempOUT[, 6] <- max(outGamma[outGamma[, 3], "gamma"])
    OUT <- rbind(OUT, tempOUT)
  }
}
colnames(OUT) <- c("Compare", "Median1", "Median2", "RelMedian", "lowerCI", "upperCI")

write.csv(OUT, paste0("Median_Ratio_", ExpNumber, ".csv"), row.names = FALSE)

#Regression Days
## Calculate total sequential days of regression; >50% animals alive; >1 consecutive time point.  
## PR >30% from staging, CR <60mm3 in NSG host mice (based on caliper thickness of skin in non-tumored animals)
## Authors: M. Konate

library(ggplot2)
options(warn = -1)

# CSV files used: therapy.csv, tumor.csv
#Set working directory, this need to be changed depending on location of working directory
DataDir <- "C:\\Users\\username\\Desktop\\R_Scripts\\Data\\datafiles\\"
setwd(DataDir)

# Enter ExpNumber below before running script. Example data studyName = ZBTJ2-2
#Set study name
studyName <- "EXPID"

#Read therapy data
txData <- read.csv(paste(studyName, "therapy.csv"), header = TRUE, stringsAsFactors = FALSE)

#Determine the number of unique treatment groups
nbGroups <- length(unique(txData$GROUP_NBR))

#Initialize output table
results <- data.frame(matrix(NA, nrow = nbGroups, ncol = 11))
colnames(results) <- c("Study", "GROUP_NBR", "Total number of animals", "Agent", "Median Staging Weight", "Median Regress >30%?", "# Days regressed >30% ", "Number of animals", "Median Regress to < 60 mm3?", "# Days regressed < 60 mm3", "Number of animals")
results[,1] <- studyName
results[,2] <- unique(txData$GROUP_NBR)

#For each group, assign NSC number(s)
for (ii in 1:nbGroups) {
  nn <- which(txData$GROUP_NBR == ii)
  if (length(nn) == 1) { results[ii,4] <- txData$NSC[nn] }
  else if (length(nn) == 2) { results[ii,4] <- paste0(txData$NSC[nn[1]], ", ", txData$NSC[nn[2]]) }
  else if (length(nn) == 3) { results[ii,4] <- paste0(txData$NSC[nn[1]], ", ", txData$NSC[nn[2]], ", ", txData$NSC[nn[3]]) }
}

#Read tumor volume data
tumorVol <- read.csv(paste(studyName, "tumor.csv"), header = TRUE, stringsAsFactors = FALSE)
tumorVol <- data.frame(tumorVol)
tumorVol <- tumorVol[(is.na(tumorVol$VALIDITY_CODE)),]

#Build table listing median tumor weight at each time point; also list number of mice at a given time point, and % of original number
summTable <- data.frame(matrix(NA, nrow = 0, ncol = 8))
colnames(summTable) <- c("STUDY", "GROUP_NBR", "NSC", "OBS_DAY", "NBR_MICE", "PCT_ORIGINAL_MICE", "MEDIAN_TUMOR_WT", "PCT_STAGING_WT")

#Row counter
nt <- 1

#Determine median tumor weight at staging for each treatment group
listGroups <- unique(tumorVol$GROUP_NBR)
for (jj in seq(listGroups)) {
  groupRows <- which(tumorVol$GROUP_NBR == listGroups[jj])
  subMat <- tumorVol[groupRows,]
  listDays <- unique(subMat$OBS_DAY)
  stagingDay <- min(listDays)
  stagingRows <- which(subMat$OBS_DAY == stagingDay)
  medianWtStaging <- median(subMat$TUMOR_WT[stagingRows])
  cutoff30 <- medianWtStaging * 0.7
  nbAnimalsOri <- length(unique(subMat$ANIMAL_NBR))#determine initial number of mice per group
  results[listGroups[jj], 3] <- nbAnimalsOri
  results[listGroups[jj],5] <- medianWtStaging
  
  #Write in table listing median tumor weight at each time point
  for (tt in 1:length(listDays)) {
    currTime <- listDays[tt]
    ttRows <- which(subMat$OBS_DAY == currTime)
    ttMedian <- median(subMat$TUMOR_WT[ttRows])
    ttNbAnimals <- length(ttRows) #determine how many animals were evaluated at the current time point
    summTable[nt,1] <- studyName
    summTable[nt,2] <- listGroups[jj] #group number
    summTable[nt,3] <- results$Agent[results$GROUP_NBR == listGroups[jj]]#NSC number
    summTable[nt,4] <- currTime #observation day
    summTable[nt,5] <- ttNbAnimals # number of mice evaluated at given time point
    summTable[nt,6] <- round(ttNbAnimals/nbAnimalsOri*100, 0) # percentage of original number of mice evaluated at given time point
    summTable[nt,7] <- ttMedian
    summTable[nt,8] <-round(ttMedian/medianWtStaging*100, 2)
    nt <- nt + 1
  }
  
  #For each treatment group, filter out the time points at which < 60% of mice are alive and only assess change in tumor size at these time points
  filtRows <- which(summTable$PCT_ORIGINAL_MICE >= 60)
  filtSummTable <- summTable[filtRows,]
  
  #for each group, determine whether tumor weight decreases by 30% or more during the experiment
  for (gp in seq(listGroups)) {
    gpRows <- which(filtSummTable$GROUP_NBR == listGroups[gp])
    subMat2 <- filtSummTable[gpRows,]
    stgWt <- subMat2$MEDIAN_TUMOR_WT[1]
    cutOff <- 0.7*stgWt #tumor decreased by at least 30%
    test <- which(subMat2$MEDIAN_TUMOR_WT <= cutOff)
    if (length(test) == 0) {
      results[listGroups[gp],6] <- "No"
      results[listGroups[gp],7] <- 0
    }
    else {
      results[listGroups[gp],6] <- "Yes"
      first <- min(subMat2$OBS_DAY[test])
      last <- max(subMat2$OBS_DAY[test])
      results[listGroups[gp],7] <- last - first
    }
    #Does median tumor volume regress to < 60 mm3?
    if (min(subMat2$MEDIAN_TUMOR_WT >= 60)) {
      results[listGroups[gp],9] <- "No"
      results[listGroups[gp],10] <- 0
    }
    else {
      results[listGroups[gp],9] <- "Yes"
      test2 <- which(subMat2$MEDIAN_TUMOR_WT < 60)
      one <- min(subMat2$OBS_DAY[test2])
      two <- max(subMat2$OBS_DAY[test2])
      results[listGroups[gp],10] <- two - one
    }
  }
  
  #Determine how many mice had tumor decrease > 30% and/or tumor volume < 60 mm3
  nbMiceEvent30 <- 0
  nbMiceEvent60 <- 0
  #Establish the valid time points for the given treatment group (i.e., time points at which at least 60% of mice in the group are alive)
  validTimes <- filtSummTable$OBS_DAY[filtSummTable$GROUP_NBR == listGroups[jj]]
  filtSubMat <- subMat[subMat$OBS_DAY %in% validTimes,]
  
  for (mm in 1:length(unique(filtSubMat$ANIMAL_NBR))) {
    currRows <- which(filtSubMat$ANIMAL_NBR == mm)
    subMat3 <- filtSubMat[currRows,]
    animalStagingWt <- subMat3$TUMOR_WT[1] #staging weight for a given mouse
    
    eventRows30 <- which(subMat3$TUMOR_WT < animalStagingWt*0.7)
    if (length(eventRows30 > 0)) { nbMiceEvent30 <- nbMiceEvent30 + 1 }
    else { nbMiceEvent30 == 0 }
    
    eventRows60 <- which(subMat3$TUMOR_WT < 60)
    if (length(eventRows60) > 0) { nbMiceEvent60 <- nbMiceEvent60 + 1}
    else { nbMiceEvent60 == 0 }
  }
  results[listGroups[jj],8] <- nbMiceEvent30
  results[listGroups[jj],11] <- nbMiceEvent60
}

#Write median tumor weights for each group and time point to a file
write.csv(summTable, file = paste(studyName, "median wts.csv"))

#Print plots of % tumor staging weight vs. ovbservation day
pdf(paste0(studyName, "_pct_tumor_wt", ".pdf"), width = 14, height = 10)
Tx_GROUP <- as.factor(summTable$NSC)
low <- 0
high <- round(max(summTable$PCT_STAGING_WT)/100,0) +1
ggplot(summTable, aes(x = OBS_DAY, y = PCT_STAGING_WT, color = Tx_GROUP)) + geom_point() + geom_line() + ggtitle(studyName) + xlab("Observation day") + ylab("% tumor staging weight") + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20), legend.text = element_text(size = 16), legend.title = element_text(size = 16), title = element_text(size = 24)) + scale_y_continuous(breaks = seq(low, high*100, 100))
dev.off()

#Write summary output to a file
write.csv(results, file = paste(studyName, "response summary.csv"))


# aAUC
#############################################################################
### aAUC analysis for SoC experiments - 22-Apr-2020.
### aAUC analysis of PDX experiments. Author L. Rubinstein and J. Subramanian
### Names of the experiments for which aAUC analysis is to be done 
### is stored in a vector (further details given as comments in the code below).
### v3 - Limit follow-up time based on inputs by Rubenstein and Evrard
### follow-up only till the time 50% or more animals alive in control group and 
### 50% or more alive in treated group, but follow up time in treated group not
### to exceed the follow-up time in the corresponding control group
### If there are groups with group nbr > 8, exclude those groups and do not
### merge their tumor data with groups 1-8 (12-Jun-2020)

### v4 - Modify code to analyze only one model at a time. Identify control groups
### by checking for presence of "99999" in the NSC name. Incorporate 
### model names in result filenames. Tabulate only for agents present in model 
### (17-Jul-2020)

### v5 - Add check to prevent t-test spitting error when there are only 0 or 1 obs 
### in a group. 
### Store a master list of primary.agents and agent.names but tabulate only for 
### agents present in model (22-Jul-2020).
#############################################################################

library(openxlsx)
library(multcomp)

## Define working directory: This is the path to the main working directory where the tumor growth 
## and therapy data are stored and results of analysis will be saved (in my laptop the data is 
## stored in subfolder named "PDM_datafiles_28Nov2019" in the case of SoC data and within "RareTumor" 
## in the case of Rare tumor data" within workdir and results are saved in 
## a subfolder named "Results"). Prior to running this code create a subfolder named "Results"
## and a subfolder named "datafiles" within workdir.

## CSV files used: therapy.csv, tumor.csv
## Set working directory, this need to be changed depending on location of working directory
workdir <- "C:\\Users\\username\\Desktop\\R_Scripts\\Data\\" 


## Store the 7 digit experiment ID, SoC.models, for which the analysis needs to be done 
## Just modify the following line by specifying the model name of interest.
## Specify only one model - 21SEPT2020

# Enter SoC.models experiment ID below before running script. Example data ExpID = ZBTJ2-2
SoC.models <- c("EXPID")

## The complete set of agents tested tumor models, expand as needed. 
## Note if the agent was given in two different dose levels for same group of animals. 
## These instances are tabulated with a "(2)" against agent names to differentiate them from the
## majority of the cases where only one dose level was given (21-Jul-2020).

primary.agents <- c("125973", "19893", "19893+19893", "241240", "241240+241240", "616348", "616348+616348", "718781", "718781+718781", "750927", "761431", "761431+761431", 
                    "0"
                    
)
agent.names <- c("Paclitaxel", "5-FU", "5-FU (2)", "Carboplatin", "Carboplatin (2)", "Irinotecan", "Irinotecan (2)", "Erlotinib", "Erlotinib (2)", "Gemcitabine", "Vemurafenib", "Vemurafenib (2)", 
                 "NoName"
                 
)

## Read therapy and tumor growth info for model of interest
## Make sure that the tumor growth and therapy files for the 
## model in Soc.Models is present in a subfolder within workdir. 
## Only one therapy and one tumor file corresponding to the model is read (17Jul2020)

therapy.files <- paste(workdir,"\\datafiles\\",SoC.models," therapy.csv",sep="")
tumor.files <- paste(workdir,"\\datafiles\\",SoC.models," tumor.csv",sep="")

## Initialize lists to store aAUC and p-value results
aAUC_GM <- list() ## to store geometric mean aAUC 
ratio.aAUC_GM <- list() ## to store GM ratio aAUC 
pval.aAUC_GM <- list() ## to store p-value for GM ratio aAUC 

pdf(paste(workdir,"\\Results\\",SoC.models," Plot.pdf",sep=""), width=11, height=8)  ## for boxplot of results 

## Do the analysis for the model of interest
therapy.info <- read.csv(therapy.files, header=TRUE)  ## read therapy groups info
therapy.gp <- aggregate(NSC~GROUP_NBR,therapy.info,paste,collapse="+")

tumorGrowth <- read.csv(tumor.files, header=TRUE, na.strings=c(""))
tumorGrowth <- tumorGrowth[is.na(tumorGrowth$VALIDITY_CODE) | tumorGrowth$VALIDITY_CODE != -1, ] ## removing obs with VALIDITY CODE = -1, 10-Jan-2019
tumorGrowth <- subset(tumorGrowth, !is.na(TUMOR_WT))

## based on requirements of study, this script expects 1 control and 7 txt groups, 8 total, therefore excluding empty groups or groups where additional studies were performed in therapy files were required for some of the models.  
## List experiment ID in appropriate line if other than 8 total groups present. Experiment IDs below are for example purposes only
if (SoC.models %in% c("ZANJ2-3","ZARJ2-2")) therapy.gp <- therapy.gp[1:6,] ## there are only 6 gps in these, therefore 2 empty rows in the therapy.info file
if (SoC.models %in% c("ZDKJ2-1","ZCHJ2-1")) therapy.gp <- therapy.gp[1:7,] ## there are only 7 gps in these, therefore 1 empty row in the therapy.info file
if (SoC.models %in% c("ZCCJ2-1")) {
  therapy.gp <- therapy.gp[1:8,] ## Removing extra Group 9, test group for training purposes
  tumorGrowth <- subset(tumorGrowth, GROUP_NBR < 9)
}
if (SoC.models %in% c("ZBHJ2-2", "ZBIJ2-1")) {
  therapy.gp <- therapy.gp[1:8,] ## Removing Group 9 and 10, test group for training purposes
  tumorGrowth <- subset(tumorGrowth, GROUP_NBR < 9)
}
if (SoC.models %in% c("ZCAJ2-1", "ZCIJ2-1")) {
  therapy.gp <- therapy.gp[1:8,] ## Removing Group 9, 10 and 11, test group for training purposes  
  tumorGrowth <- subset(tumorGrowth, GROUP_NBR < 9)
}
if (SoC.models %in% c("ZCDJ2-1")) {
  therapy.gp <- therapy.gp[-c(4,7,8),] ## groups 4,7,8 here have empty rows in the therapy.info file due to censored study data
  therapy.gp$GROUP_NBR <- 1:5 ## renumber sequentially
  tumorGrowth$GROUP_NBR <- ifelse((tumorGrowth$GROUP_NBR > 4), tumorGrowth$GROUP_NBR-1, tumorGrowth$GROUP_NBR) ## renumber in tumor growth data also
}

tumorGrowth$GROUP_NBR <- factor(tumorGrowth$GROUP_NBR)
tumorGrowth$ANIMAL_NBR <- factor(tumorGrowth$ANIMAL_NBR)
tumorGrowth$ANIMAL_ID_NEW <- factor(paste(tumorGrowth$GROUP_NBR,tumorGrowth$ANIMAL_NBR,sep="_"))

## order by group id, within each group by animal id and within each animal by obs day
tumorGrowth <- tumorGrowth[order(tumorGrowth$GROUP_NBR, tumorGrowth$ANIMAL_NBR, tumorGrowth$OBS_DAY),] 
n.gps <- length(unique(tumorGrowth$GROUP_NBR)) ## get the number of experimental groups

## preparing control group index for comparison with the correct control. NSC numbers for control material predefined to start with “99999”
controls <- as.numeric(grepl("99999", therapy.gp$NSC)) 
for(i in 1:n.gps)
  if(controls[i]) controls[i] <- i else controls[i] <- controls[i-1]

## 10-Jun-2020: Cut-off time imposed for assessing response, based on suggestions from Rubenstein and Evrard
follow_up_time <- vector(mode="integer",length=n.gps) ## to store max follow-up time for each group
for(group in 1:n.gps){
  data.group <- subset(tumorGrowth, GROUP_NBR==group)
  n.animals <- length(unique(data.group$ANIMAL_NBR))
  fups <- vector(mode="integer",length=n.animals) ## follow up time for each animal
  for(animal in 1:n.animals) {
    if(length(subset(data.group, ANIMAL_NBR==animal)$OBS_DAY) > 1){ ## if post baseline follow up exists
      fups[animal] <- max(subset(data.group, ANIMAL_NBR==animal)$OBS_DAY) ## find the follow up duration
    } else {
      n.animals <- n.animals-1 ## Animals with only baseline data are excluded
    }
  }
  fups <- fups[which(fups > 0)]
  ## at least 50% animals should be alive in a group for analysis to be carried out
  ## 0.001 added to get the NEXT higher integer if n/2 is integer
  follow_up_index <- ifelse(n.animals > 1, ceiling((n.animals/2)+0.001), NA) 
  follow_up_time[group] <- fups[order(fups)[follow_up_index]] ## time restricted to the point at which the % animals alive becomes less than 50
}
follow_up_time <- ifelse(follow_up_time > follow_up_time[controls], follow_up_time[controls], follow_up_time) ## no treated group is followed up beyond corresponding controls
tumorGrowth.limitedFUP <- do.call(rbind, lapply(1:n.gps, function(group) 
  subset(tumorGrowth, GROUP_NBR==group & OBS_DAY <= follow_up_time[group])))## groups with only 0 or 1 animal get removed here

## Normalize tumor weight and calculate aAUC from normalized weights
tumorGrowth.norm <- do.call(rbind, by(tumorGrowth.limitedFUP, tumorGrowth.limitedFUP[,c("ANIMAL_NBR", "GROUP_NBR")], ## corrected to ANIMAL+GROUP NBR on 5Apr2017
                                      function(subset) within(subset, 
                                                              { TUMOR_WT.norm <- (TUMOR_WT/TUMOR_WT[1])
                                                              OBS_DAY.norm <- (OBS_DAY - OBS_DAY[1]) 
                                                              
                                                              ## aAUC
                                                              ll <- nrow(subset)
                                                              if(ll > 1) {## at least one non-baseline obs required for AUC calculation
                                                                partial_auc <- (OBS_DAY.norm[2:ll] - OBS_DAY.norm[1:(ll-1)])*
                                                                  ((TUMOR_WT.norm[1:(ll-1)] + TUMOR_WT.norm[2:ll]) / 2) 
                                                                partial_auc[ll] <- 0
                                                                aAUC <- sum(partial_auc)/OBS_DAY.norm[ll]
                                                              } else { ## if only the baseline value is available, set aAUC=NA, this animal is not counted for t tests and CIs
                                                                partial_auc <- 0
                                                                aAUC <- NA
                                                              }
                                                              })))
rownames(tumorGrowth.norm) <- NULL
tumorGrowth.norm <- tumorGrowth.norm[,c("EXP_NBR","GROUP_NBR","ANIMAL_NBR","ANIMAL_ID_NEW","OBS_DAY",
                                        "TUMOR_LEN","TUMOR_WID","OBS_DAY.norm","TUMOR_WT.norm",
                                        "partial_auc","aAUC")]

## AUC dataset
AUCs <- tumorGrowth.norm[,c("EXP_NBR","GROUP_NBR","ANIMAL_NBR","ANIMAL_ID_NEW","aAUC")]  
AUCs <- AUCs[!duplicated(AUCs$ANIMAL_ID_NEW),]

## calculate GM of aAUC for every group
## aAUC of groups with < 50% animals post baseline are set to NA (22-Jul-2020) 
## (these groups are already removed from AUCs dataset)
aAUC_GM[[SoC.models]] <- sapply(1:n.gps,function(x) 
  ifelse(x %in% unique(AUCs$GROUP_NBR), exp(mean(log(subset(AUCs, GROUP_NBR==x)$aAUC), na.rm=TRUE)), NA))
names(aAUC_GM[[SoC.models]]) <- therapy.gp$NSC
## Ratio of GM to corr controls
ratio.aAUC_GM[[SoC.models]] <- aAUC_GM[[SoC.models]]/aAUC_GM[[SoC.models]][controls]
## Welch's test p-value
pval.aAUC_GM[[SoC.models]] <- sapply(1:n.gps, function(i){
  if(i %in% unique(AUCs$GROUP_NBR)){
    j <- controls[i]
    ifelse(i != j, t.test(log(aAUC) ~ GROUP_NBR, 
                          data=AUCs, subset=(GROUP_NBR %in% (levels(GROUP_NBR))[c(i,j)]))$p.value, NA) 
  } else NA
})
names(pval.aAUC_GM[[SoC.models]]) <- therapy.gp$NSC

## Plot
par(las=3, mar=c(8,4,4,2))
col <- ifelse(substr(therapy.gp$NSC,1,5) == "99999", "white","grey")
plot(aAUC~GROUP_NBR, data=AUCs, xlab="", ylab="aAUC", main=SoC.models, lwd=2, col=col, range=0,
     names=therapy.gp$NSC, cex.main=1.1, cex.lab=1.2, cex.axis=1.2)
dev.off() ## close plot file

## Tabulating GM ratio and its p-value for all agents in the experiment
## Modified so that tabulation will be only for the agents present in current model, instead of all agents in primary.agents - 17Jul2020
curr.primary.agents <- primary.agents[which(primary.agents %in% therapy.gp[-unique(controls),"NSC"])]
curr.agent.names <- agent.names[which(primary.agents %in% therapy.gp[-unique(controls),"NSC"])]

corder <- order(c(rep(1:length(curr.primary.agents), 2))) ## each agent has 1 ratio and 1 pval so total 2 quantities
ratio.aAUC_GM <- do.call(rbind, lapply(ratio.aAUC_GM, function(x) x[curr.primary.agents]))
pval.aAUC_GM <- do.call(rbind, lapply(pval.aAUC_GM, function(x) x[curr.primary.agents]))
colnames(ratio.aAUC_GM) <- paste(curr.primary.agents, curr.agent.names, "Ratio aAUC_GM", sep="\n")
colnames(pval.aAUC_GM) <- paste(curr.primary.agents, curr.agent.names, "pval aAUC_GM", sep="\n")

full.table <- data.frame(ratio.aAUC_GM, pval.aAUC_GM, check.names=FALSE, stringsAsFactors=FALSE)[,corder]

## Format and write to excel file
wb <- createWorkbook() # creates workbook
addWorksheet(wb, "Sheet1") # adds sheet
options("openxlsx.numFmt" = "0.0000") # 4 decimal cases formating
writeData(wb, sheet=1, keepNA=TRUE, na.string="NA", rowNames=TRUE, full.table) # writing content
cs1 <- createStyle(valign="center")
cs2 <- createStyle(textDecoration="bold", wrapText=TRUE)
cs3 <- createStyle(border="bottom", borderStyle="thick")
cs4 <- createStyle(border="right", borderStyle="thick")

addStyle(wb, sheet=1, style=cs1, rows=1, cols=1:nrow(full.table))
addStyle(wb, sheet=1, style=cs2, rows=1, cols=2:(ncol(full.table)+1), stack=TRUE)
setColWidths(wb, sheet=1, cols=1, widths="auto")
setColWidths(wb, sheet=1, cols=2:(ncol(full.table)+1), widths=16)
addStyle(wb, sheet=1, style=cs3, rows=1, cols=1:(ncol(full.table)+1), stack=TRUE)
addStyle(wb, sheet=1, style=cs4, rows=1:(nrow(full.table)+1), cols=(0:length(curr.primary.agents)*2+1), 
         gridExpand=TRUE, stack=TRUE)


saveWorkbook(wb, paste(workdir, "\\Results\\",SoC.models," tabulation.xlsx",sep=""), overwrite=TRUE) # save workbook

