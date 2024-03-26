library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyverse)
library(plyr)
library(reshape2)

# Function definition
colfunc <- colorRampPalette(c("#216285", "#7cb3d2"))
fscore <- function (dfConfMat) {
  dfConfMat[is.na(dfConfMat),"count"] = 0
  return((2*dfConfMat['TP','count'])/(2*dfConfMat['TP','count']+dfConfMat['FP','count']+dfConfMat['FN','count']))
}
mccscore <- function (dfConfMat) {
  dfConfMat[is.na(dfConfMat),"count"] = 0
  return(((dfConfMat['TP','count']*dfConfMat['TN','count'])
         - (dfConfMat['FP','count']*dfConfMat['FN','count']))
         / (sqrt((dfConfMat['TP','count'] + dfConfMat['FP','count'])
         *(dfConfMat['TP','count'] + dfConfMat['FN','count'])
         *(dfConfMat['TN','count'] + dfConfMat['FP','count'])
         *(dfConfMat['TN','count'] + dfConfMat['FN','count']))))
}

# Reading the data
results <- read.csv("all_predictions.tsv", sep="\t")
times <- read.csv("realTimes.tsv", sep="\t")

## Format the datasets ID according to the benchmarking mode
if ("proportion.mode" %in% list.files()){
  ids <- levels(factor(results[,"ID"], levels=str_sort(levels(factor(results[,"ID"])), numeric=T), ordered=T))
  for (i in 1:length(ids)) {
    results[results$ID == ids[i], "ID"] = paste("DS", i, sep="")
    times[times$X == ids[i], "X"] = paste("DS", i, sep="")
  }
  results[,"ID"] <- factor(results$ID, levels = str_sort(unique(results$ID), numeric=T))
  }

# In the case of length mode, the mean of each range is used to draw the points
if ("length.mode" %in% list.files()){
  ids = levels(as.factor(results[,"ID"]))
  for (i in 1:length(ids)) {
    results[results$ID == ids[i], "ID"] = mean(as.numeric(str_replace(strsplit(results[results$ID == ids[i], "ID"], "to")[[1]], "pb", "")))
    times[times$X == ids[i], "X"] = strsplit(times[times$X == ids[i], "X"], "to")[[1]][1]}
  results[,"ID"] <- factor(results$ID)
}


## Calculate the metrics of the Confusion Matrix
scores = data.frame(ID = "", Tool="", F1score="", MCCscore="")

for (tool in colnames(results)[c(-1, -length(colnames(results)), -(length(colnames(results))-1), -(length(colnames(results))-2))]) {
  dfTruthCheck <- cbind(as.character(results$ID), results[,1], paste(results[,tool] == results[,"division"], results[,tool] == "eukaryote"))
  colnames(dfTruthCheck) <- c("ID", "Contig ID", "Result")
  dfTruthCheck[dfTruthCheck[,3]=="TRUE FALSE",3] = "TN"
  dfTruthCheck[dfTruthCheck[,3]=="FALSE FALSE",3] = "FN"
  dfTruthCheck[dfTruthCheck[,3]=="TRUE TRUE",3] = "TP"
  dfTruthCheck[dfTruthCheck[,3]=="FALSE TRUE",3] = "FP"
  for (id in unique(results$ID)) {
    TP <- as.integer(table(dfTruthCheck[dfTruthCheck[,1]==id,3])["TP"])
    TN <- as.integer(table(dfTruthCheck[dfTruthCheck[,1]==id,3])["TN"])
    FP <- as.integer(table(dfTruthCheck[dfTruthCheck[,1]==id,3])["FP"])
    FN <- as.integer(table(dfTruthCheck[dfTruthCheck[,1]==id,3])["FN"])
    dfConfMat <- data.frame(c(TP, TN, FP, FN), row.names = c("TP", "TN", "FP", "FN"))
    colnames(dfConfMat) <- "count"
    scores <- rbind(scores, c(id, tool, fscore(dfConfMat), mccscore(dfConfMat)))}}

scores[, "F1score"] <- as.numeric(scores[, "F1score"]) 
scores[, "MCCscore"] <- as.numeric(scores[, "MCCscore"]) 
scores[, "Tool"] <- factor(scores[, "Tool"], levels= colnames(results)[c(-1, -length(colnames(results)), -(length(colnames(results))-1), -(length(colnames(results))-2))]) 
scores <- scores[-1,]
results[,"division"] <- as.factor(paste(results$division, results$subdivision, sep=" "))

# Plot the data according to the benchmarking mode
if ("proportion.mode" %in% list.files()){
  scores$ID <- factor(scores$ID)
  levels(scores$ID) <- paste("DS", 1:length(unique(scores$ID)), sep = "")
  
  # Format the division percentage
  percent <- as.data.frame(table(results$division, results$ID), margin=2)
  percent$Var2 <- factor(percent$Var2, levels=str_sort(unique(results$ID), numeric=T))
  for (i in 1:length(levels(results$ID))) {
    percent[percent$Var2==levels(percent$Var2)[i],"n"] <- sum(percent[percent$Var2==levels(percent$Var2)[i],"Freq"])
  }
  colnames(percent) <- c("Division", "ID", "Count", "CountID")
  levels(percent$ID) <- paste("DS", 1:length(unique(percent$ID)), sep = "")
  
  divisionPlot <- ggplot(percent, aes(x=ID, y=Count, fill=Division)) + 
    geom_bar(stat="identity") + 
    geom_text(aes(label= paste(round((Count/CountID)*100, 2), "%")), position=position_stack(vjust=0.5), size=4, col="white") +
    xlab("") +
    ylab("Number of contigs") +
    scale_fill_manual("Taxonomy", values=c(colfunc(length(unique(results$division))-1), "#8283a6")) +
    theme(legend.position="bottom")
  
  # Dodge barplot with the F1 score 
  fscorePlot <- ggplot(scores, aes(x=ID, y=F1score, fill=Tool)) + 
    geom_bar(position="dodge2", stat="identity") + 
    xlab("Datasets") +
    ylab("F1-Score") +
    geom_text(aes(label = round(F1score, 3)), position = position_dodge(width = .9), hjust=1.5, size=3, angle=90) +
    ggtitle("F1-Score per tool on several datasets (positive reference = eukaryote)") +
    scale_fill_manual("Tools", values=c("#eab595","#d87f81", "#ae6378", "#79616f", "#7e9680")) +
    ylim(0, 1) + 
    theme(legend.position="top")
  
  # Dodge barplot with the MCC score
  mccscorePlot <- ggplot(scores, aes(x=ID, y=MCCscore, fill=Tool)) + 
    geom_bar(position="dodge2", stat="identity") + 
    xlab("Datasets") +
    ylab("MCC score") +
    geom_text(aes(label = round(MCCscore, 3)), position = position_dodge(width = .9), hjust=1.5, size=3, angle=90) +
    ggtitle("MCC score per tool on several datasets (positive reference = eukaryote)") +
    scale_fill_manual("Tools", values=c("#eab595","#d87f81", "#ae6378", "#79616f", "#7e9680")) +
    ylim(0, 1) + 
    theme(legend.position="top")
  
  ### Times ####
  times[,"ID"] <- factor(times$X, levels=str_sort(unique(results$ID), numeric=T))
  times <- times[,-1]
  times <- melt(times, value.name = "Time", variable.name = "Tool")
  
  # Stacked barplot with the division
  timeplot <- ggplot(times, aes(x=ID, y=Time, fill=Tool)) + 
    geom_bar(position="dodge2", stat="identity") + 
    geom_text(aes(label = round(Time, 2)), position = position_dodge(width = .9), hjust=1.5, size=3, angle=90) +
    xlab("Datasets") +
    ylab("Elapsed time (s)") +
    ggtitle("Time per tool on several datasets") +
    scale_fill_manual("Tools", values=c("#eab595","#d87f81", "#ae6378", "#79616f", "#7e9680")) +
    theme(legend.position="top")
  
  print(ggarrange(fscorePlot,divisionPlot, nrow=2, heights = c(2, 1)))
  print(ggarrange(mccscorePlot,divisionPlot, nrow=2, heights = c(2, 1)))
  print(ggarrange(timeplot,divisionPlot, nrow=2, heights = c(2, 1)))}

if ("length.mode" %in% list.files()){
  results[,"division"] <- as.factor(paste(results$division, results$subdivision, sep=" "))
  scores[,"ID"] <- as.numeric(scores[,"ID"])
  times[,"ID"] <- as.factor(times$X)
  times = times[,-1]
  times = melt(times, value.name = "Time", variable.name = "Tool")
  times$ID = as.numeric(as.character(times$ID))
  
  fscorePlot <- ggplot(scores, aes(x=ID, y=F1score, group=Tool)) +
    geom_line(aes(color=Tool), size=1.1)+
    scale_color_manual("Tools", values=c("#eab595","#d87f81", "#ae6378", "#79616f", "#7e9680")) +
    xlab("Length") +
    ylab("F1-Score") +
    ylim(0.6, 1) +
    ggtitle("F1-Score per tool per length")
  
  mccscorePlot <- ggplot(scores, aes(x=ID, y=MCCscore, group=Tool)) +
    geom_line(aes(color=Tool), size=1.1)+
    scale_color_manual("Tools", values=c("#eab595","#d87f81", "#ae6378", "#79616f", "#7e9680")) +
    xlab("Length") +
    ylab("MCC Score") +
    ylim(0.6, 1) +
    ggtitle("MCC score per tool per length")
    
  timePlot <- ggplot(times, aes(x=ID, y=Time, group=Tool)) + 
      geom_line(aes(color=Tool), size=1.1)+
      scale_color_manual("Tools", values=c("#eab595","#d87f81", "#ae6378", "#79616f", "#7e9680")) +
      xlab("Length") +
      ylab("Elapsed time (s)") +
    ggtitle("Time per tool per length")
    
  print(fscorePlot)
  print(mccscorePlot)
  print(timePlot)
  }