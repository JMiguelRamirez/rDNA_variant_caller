#  Libraries 
library(ggplot2)

#Set working directory to souce file location
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

variables <- c("", "_rep", "_non_rep", "_af")
# variables <- c("", "_af")
# variables <- c("")
for(variable in variables){
  merged_scores <- read.delim(paste0("Plots/F1", variable, ".tab"), sep = " ")
  colnames(merged_scores) <- c("Precision", "Recall", "F1_score", "Simulation", "ploidy")
  
  merged_scores$ploidy <- factor(merged_scores$ploidy, levels=c("Ploidy_2", "Ploidy_5", "Ploidy_10", "Ploidy_20", "Ploidy_30", "Ploidy_40", 
                                                                "Joint_2", "Joint_5", "Joint_10", "Joint_20",
                                                                "Mutect", "Lofreq"))
  #Plotting:
  palette <- c(  "Ploidy_2"="orange", "Ploidy_5"="orangered", "Ploidy_10"="#FF6666" , "Ploidy_20"="#9e0142",
                 "Ploidy_30"="#B00B69", "Ploidy_40"="red3", 
                 "Joint_2"="#C0EEFF", "Joint_5"="#339999",
                 "Joint_10"= "royalblue", "Joint_20"="#3333FF",
                 "Mutect"="orchid2","Lofreq"="purple2") 
  
  merged_scores$F1_score <- as.numeric(merged_scores$F1_score)
  g <- ggplot(merged_scores) + geom_boxplot(aes(ploidy, F1_score, fill=ploidy), outlier.shape = NA, fatten = 1.5) + theme_bw() + xlab(NULL) +
    guides(fill="none") + scale_fill_manual(values=palette) + ylim(c(0.6,1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=11),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14, vjust = 3)) +
    scale_x_discrete(labels=c('HC - 2', 'HC - 5', 'HC - 10', 'HC - 20',
                              'HC - 30', 'HC - 40', 
                              'JG - 2', 'JG - 5', 'JG - 10', 'JG - 20',
                              'Mutect2', 'Lofreq')) +
    ylab("F1 score")
  g
  round(mean(merged_scores[merged_scores$ploidy=="Ploidy_20","F1_score"]),2)
  
  ggsave(paste0("Plots/F1", variable, ".png"), g, width = 3, height = 2.6, device = "png")
  ggsave(paste0("Plots/F1", variable, ".pdf"), g, width = 3, height = 2.6, device = "pdf")
  
  #fatten reduces the thickness of the median line in the boxplot, outlier.size reduces the size of the point outside the boxplot
  merged_scores$Precision <- as.numeric(merged_scores$Precision)
  g2 <- ggplot(merged_scores) + geom_boxplot(aes(ploidy, Precision, fill=ploidy), fatten = 1.5, outlier.shape = NA) + theme_bw() + xlab(NULL) +
    guides(fill="none") + scale_fill_manual(values=palette) + ylim(c(0.6,1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=11),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14, vjust = 3)) +
    scale_x_discrete(labels=c('HC - 2', 'HC - 5', 'HC - 10', 'HC - 20',
                              'HC - 30', 'HC - 40', 
                              'JG - 2', 'JG - 5', 'JG - 10', 'JG - 20',
                              'Mutect2', 'Lofreq')) +
    ylab("Precision")
  g2
  ggsave(paste0("Plots/Precision", variable, ".png"), g2, width = 3, height = 2.6, device = "png")
  ggsave(paste0("Plots/Precision", variable, ".pdf"), g2, width = 3, height = 2.6, device = "pdf")
  
  merged_scores$Recall <- as.numeric(merged_scores$Recall)
  g3 <- ggplot(merged_scores) + geom_boxplot(aes(ploidy, Recall, fill=ploidy), fatten = 1.5, outlier.shape = NA) + theme_bw() + xlab(NULL) +
    guides(fill="none") + scale_fill_manual(values=palette) + ylim(c(0.6,1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=11),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14, vjust = 3)) +
    scale_x_discrete(labels=c('HC - 2', 'HC - 5', 'HC - 10', 'HC - 20',
                              'HC - 30', 'HC - 40', 
                              'JG - 2', 'JG - 5', 'JG - 10', 'JG - 20',
                              'Mutect2', 'Lofreq')) +
    ylab("Recall")
  g3
  ggsave(paste0("Plots/Recall", variable, ".png"), g3, width = 3, height = 2.6, device = "png")
  ggsave(paste0("Plots/Recall", variable, ".pdf"), g3, width = 3, height = 2.6, device = "pdf")
  
  # round(mean(merged_scores[merged_scores$ploidy==20,4]),2) #Mean F1
  # round(mean(merged_scores[merged_scores$ploidy==20,3]),2) #Mean  recall
  # round(mean(merged_scores[merged_scores$ploidy==20,2]),2) #Mean precision
  
}



#Plot HC 20 for repetitive regions only:
hc_20 <- read.delim(paste0("Plots/F1_rep.tab"), sep = " ")
colnames(hc_20) <- c("Precision", "Recall", "F1_score", "Simulation", "ploidy")

hc_20 <- hc_20[hc_20$ploidy=="Ploidy_20",]
hc_20 <- hc_20[,c(1:4)] #Without the ploidy column

#Mean F1 score of repetitive regions:
round(mean(hc_20$F1_score), 2)
round(mean(hc_20$Precision), 2)
round(mean(hc_20$Recall), 2)

library(tidyr) #from wide to long format
long_data <- pivot_longer(
  hc_20,
  cols = c("Precision", "Recall", "F1_score"),
  names_to = "Statistic",   
  values_to = "Value"  
)

palette <- c(  "Precision"="#e5729b", "Recall"="#a9a3c5", "F1_score"="orchid3") 

long_data$Statistic <- factor(long_data$Statistic, levels=c("Precision", "Recall", "F1_score"))

g4 <- ggplot(long_data) + geom_boxplot(aes(Statistic, Value, fill=Statistic), fatten = 1.5, outlier.size=0.5) + theme_bw() + xlab(NULL) +
  guides(fill="none") + scale_fill_manual(values=palette) + ylim(c(0.6, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, vjust = 3)) +
  scale_x_discrete(labels=c('Precision', 'Recall', 'F1 score')) +
  ylab("Statistic")
g4
ggsave(paste0("Plots/HC_rep.png"), g4, width = 2, height = 2.6, device = "png")
ggsave(paste0("Plots/HC_rep.pdf"), g4, width = 2, height = 2.6, device = "pdf")
