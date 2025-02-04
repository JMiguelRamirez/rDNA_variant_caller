#Code to compute nucleotide diversity, as most methods do not work with polyploidy
#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez; Adapted by Winona Oliveros
# @E-mail: jose.ramirez1@bsc.es
# @software version: R=4.2.2

rm(list=ls())
#Load libraries
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)

# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Reading GTEx variants
data <- read.delim("1000G/merged.txt")

#Computing allele frequency per donor:
computing_allele_frequency <- function(data){
  
  #iterating over positions to get a new data frame with the allele frequencies
  new_colnames <- gsub(".GT", ".AF", colnames(data)[grep(".GT", colnames(data), fixed=T)], fixed = T)
  output <- matrix(data=NA, nrow=nrow(data), ncol=length(new_colnames)+3)
  colnames(output) <- c("POS", "REF", "ALT", new_colnames)
  output <- as.data.frame(output)
  output[,c(1:3)] <- data[,c(2:4)]
  
  for(pos in unique(data$POS)){
    print(pos)
    subset <- data[data$POS==pos,]
    ref <- subset$REF #All references are essentially the same, even if they are written differently
    alt <- subset$ALT 
    
    for(sample in seq(6, ncol(data), by = 2)){ #Computing AFs per donor
      if(sum(is.na(subset[,sample]))==nrow(subset)){ #This means that the positions is not called in the donor
        next
      } else{ #The positions has been called for at least one alternative allele
        donor <- gsub(".AD", ".AF", colnames(data)[sample], fixed = T)
        alt_allele_vector <- subset[,sample] #vector with counts of all alternatives
        alt_allele_vector[is.na(alt_allele_vector)] <- "0,0"
        alt_allele_vector <- do.call(rbind, strsplit(alt_allele_vector, ","))
        ref_counts <- as.numeric(unique(alt_allele_vector[alt_allele_vector[,1]!="0",1])) #All reference have the same number of counts
        alt_counts <- as.numeric(alt_allele_vector[,2])
        ref_freq <- ref_counts/sum(ref_counts, alt_counts)
        alt_freq <- alt_counts/sum(ref_counts, alt_counts)
        
        output[output$POS==pos,donor] <- paste0(ref_freq, ",", alt_freq)
      }
    }
  }
  return(output)
}

af <- computing_allele_frequency(data)


# Computing nucleotide diversity score:
# N is the Number of comparisons. The value to normalize for if we used the original formula: Number of pairwise differences / Number of combinations
# n_samples <- ncol(af) -3
# N <- (n_samples) * (n_samples-1) / 2

pis <- data.frame()
for(pos in unique(af$POS)){ #We get one value per variant and then we average over certain regions
  print(pos)
  subset <- af[af$POS==pos,]
  frequencies <- subset[,-c(1:3)]
  frequencies[is.na(frequencies)] <- "1,0"
  frequencies <- as.data.frame(t(frequencies))
  freq <- reduce(colnames(frequencies), function(df, col) { #Split a column of type REF,ALT into two columns: REF and ALT 
    separate(df, col, into = paste0(col, c("_REF", "_ALT")), sep = ",", remove = TRUE)
  }, .init = frequencies)
  
  freq <- freq[,c(1,grep("ALT",colnames(freq)))] #keep only one reference allele
  freq <- freq %>% #Transformin into numeric
    mutate(across(everything(), as.numeric))
  
  mean_af <- colMeans(freq) #Computing mean allele frequency per allele
  
  #Traditionally, for nucleotide diversity, the number of pairwise difference is computed, but to get the equivalent using frequencies we compute:
  #2*(all pairwise combinatios of allele products)
  #pi <- 2*sum(combn(mean_af, 2, prod))
  #which is equivalent to:
  pi <- 1 - sum(mean_af^2)
  # pi <- pairwise_differences/N #The original formula would be number of pairwise differences divided by number of pairwise combinations
  pis <- rbind(pis, c(pos, pi))
  
}


# Compare average nucleotide diversity in regions, as it is more reliable than per variant
assign_region <- function(position){
  # <= only in IGS at start but all other should be > only but then variants are 0 bed based
  #so variant at the end of region be annotated to be in wrong region so we should use >=
  if(position <= 9338){
    "IGS"
  }else if(position <= 12995){
    "5'ETS"
  }else if(position <= 14864){
    "18S"
  }else if(position <= 15934){
    "ITS1"
  }else if(position <= 16091){
    "5.8S"
  }else if(position <= 17258){
    "ITS2"
  }else if(position <= 22309){
    "28S"
  }else if(position <= 22670){
    "3'ETS"
  }else{
    "IGS"
  }
}

colnames(pis) <- c("Pos", "Pi")

data <- pis %>%
  rowwise() %>%
  mutate(Region = assign_region(Pos))


### Reading bed file containing the regions of the chrR that fall into each region
coordinates_rDNA_regions <- read.table("../../Pablo/SharedWork/GTEx/Ploidy20/coordinates_regions_rDNA", header = TRUE)

coordinates_rDNA_regions$Start <- coordinates_rDNA_regions$Start +1

for(row in 1:nrow(coordinates_rDNA_regions)){
  for(pos in seq(coordinates_rDNA_regions$Start[row], coordinates_rDNA_regions$End[row])){
    if(!pos %in% data$Pos){
      data <- rbind(data, c(pos, 0, NA)) #For the region I would use the other function as here some intervals overlap different regions
    }
  }
}

data <- data %>%
  rowwise() %>%
  mutate(Region = assign_region(Pos))

data$Region <- factor(data$Region, levels = unique(data$Region))
colors_regions <- c("5'ETS" = "#5A0002","18S" = "#A40606","ITS1" = "#437C90","5.8S" = "#D98324","ITS2" = "#255957","28S" = "#60D394", "3'ETS" = "#DAD6D6")

average <- data %>%
  group_by(Region) %>%
  summarize(Mean = mean(Pi, na.rm=TRUE))

#Creating and saving the plot
plot <- ggplot(average[average$Region!='IGS',]) + geom_col(aes(Region, Mean, fill=Region)) + ylab("Mean nucleotide diversity") + 
  scale_fill_manual(values=colors_regions) + theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 13), # Increase x-axis title size
    axis.title.y = element_text(size = 13), # Increase y-axis title size
    plot.title = element_text(size = 16),   # Increase plot title size
    legend.position = "none"               # Place legend on the right
  ) 
pdf("barplot_nucleotide_diversity_per_region.pdf", width = 4, height = 3)

print(plot)

dev.off()


