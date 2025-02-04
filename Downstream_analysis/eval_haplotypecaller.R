#Code to assing FP, TP or FN to the variant calling

## libraries 
library(data.table) #More efficient data handling
library(dplyr)

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Read vcf file from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
# input_file <- "Simulations_10/var/jobs/GATK_all/Sample98.20_filtering_norm_called.vcf.gz"
# input_file <- "Simulations_3/var/jobs/GATK_all/Sample0.20_filtering_norm_called.vcf.gz"
vcf <- fread(input_file)

simulation <- sub("/.*", "", input_file)
sample <- sub(".*/Sample(\\d+)\\.\\d+_filtering_norm_called\\.vcf\\.gz", "Sample\\1", input_file) #\\d+ matches any number
colnames(vcf)[10] <- "VALUES"

#Get allele frequency
vcf <- vcf %>%
  mutate(
    # Split the column into components
    parts = strsplit(VALUES, "[:,]"),
    value1 = sapply(parts, function(x) as.numeric(x[3])),
    value2 = sapply(parts, function(x) as.numeric(x[4])),
    AF = value1 / value2
  )  %>%
  select(-c("#CHROM", ID, QUAL, FILTER, INFO, FORMAT, VALUES, parts, value1, value2))


#Get the true variants for the sample and simulation:
dir <- paste0(simulation, "/var/jobs/Catalogue/")
### Get true_variants ######
true_variants <- read.delim(paste0(simulation, "/var/jobs/Catalogue/true_variants.tbl"), sep=" ")


#Annotate variants as TP, FP or FN:
# Create variant IDs 
vcf <- vcf %>%
  mutate(variant = paste(POS, REF, ALT, sep = "_"))

true_variants <- true_variants %>%
  filter(Sample == sample) %>%
  mutate(variant = paste(POS, REF, ALT, sep = "_"))


# Annotate FPs and TPs
result_1 <- vcf %>%
  mutate(
    Status = case_when(
      variant %in% true_variants$variant ~ "TP",  # True Positive
      TRUE ~ "FP"                                 # False Positive
    )
  )
#Add whether the variants in result_1 are in what we call repetitive or non-repetitive
bed <- fread("data/pre-rRNA_47S.included_intervals.bed")[,c(1:4)]
colnames(bed) <- c("POS", "START", "END", "REGION")

# Check if positions in the vcf fall within any of the non-repetitive regions in the BED file
result_1 <- result_1 %>%
  mutate(
    Callable = ifelse(
      mapply(function(pos) {
        # Check if the position is within the non-repetitive regions (defined in BED)
        any(pos >= bed$START & pos <= bed$END)
      }, result_1$POS),
      "callable",
      "blacklisted"
    )
  )


# Identify False Negatives
result_2 <- true_variants %>%
  filter(!variant %in% vcf$variant) %>%
  mutate(Status = "FN")
result_2$AF <- 0
result_2 <- result_2[,c("POS", "REF", "ALT", "AF", "variant", "Status","Class")]
colnames(result_2)[7] <- "Callable"

# Combine results
final_result <- rbind(result_1, result_2)


# tp <- sum(final_result$Status=="TP") #192
# fp <- sum(final_result$Status=="FP") #16
# fn <- sum(final_result$Status=="FN") #25

#Fixing INDELs that might be differentially annotated:
fps <- final_result[final_result$Status=="FP"]

for(row in 1:nrow(fps)){ #row=3
  ref <- as.character(fps[row,2])
  alt <- as.character(fps[row,3])
  #Check insertions: 
  #GATK and Lofreq might get "confused" with some insertions if the nucleotides after the insertion might match the first nucleotides that are inserted. 
  #This will only affects the simulated data, with real data we will call a variant that is real, the only thing that changes is in which position we call it.
  #I will fix this as they are named as FP but they are TPs.
  #In insertions we can just ignore the first nucleotide of the insertion (as we manually put it ther because it is was usually gatk does), and then
  #depending on the number of positions shifted, we append the nucleotides shifted at the end, and then we add the nucleotide that is now the first, as the end as well
  for(shifted in 1:5){ #I added as maximum a shift of 5, even though the maximum shift that happened in an indel positions is 3, but we iterate until 5 just in case
    candidate_alt <- paste0(paste0(strsplit(alt, "")[[1]][-c(1:shifted)], collapse = ""), paste0(strsplit(alt, "")[[1]][c(1:shifted+1)], collapse = "")) #Note that the way they would be added are the inverse, so not 1 and 2, but 2 and 1
    candidate_pos <- as.numeric(fps[row,1]) + shifted
    if(sum(true_variants$POS==candidate_pos & true_variants$ALT==candidate_alt)>0){
      fps[row,6] <- "TP" #The position is actually a TP
      final_result <- final_result[!(final_result$POS==candidate_pos & final_result$ALT==candidate_alt),] #remove the FN
      
    }
  }
  #Now, we can check the deletions:
  if(nchar(ref)<2){next} #In this case, for sure it is not a deletion, so we can skip it to avoid generating NAs 
  #I will create a fancy iteration such as in insertions, because I know we need it
  
  #For deletions we do something similar as with insertions but with the reference allele. We remove the first one, and we will always add at the end the same as the new first. Then we can have different shiftings
  #Examples: 10396 -> 10398, 14675 -> 14676, 15367 -> 15369, 20186 -> 20187 # test <- fps[c(1,6,8,15), ] # test_true <- true_variants[c(7,96,102,181),]
  for(shifted in 1:5){ #I added as maximum a shift of 5, even though the maximum shift that happened in an indel positions is 3, but we iterate until 5 just in case
    candidate_ref <- paste0(paste0(strsplit(ref, "")[[1]][-c(1:shifted)], collapse = ""), paste0(strsplit(ref, "")[[1]][c(1:shifted+1)], collapse = "")) #An example of deletions shifted once: instead of GTCCC > G we get TCCCT > T. This would happen if int he reference after the insertion we had a T, so this deletion reporting would be preferred by GATK
    candidate_alt <- strsplit(candidate_ref, "")[[1]][1]
    candidate_pos <- as.numeric(fps[row,1]) + shifted
    if(sum(true_variants$POS==candidate_pos & true_variants$REF==candidate_ref & true_variants$ALT==candidate_alt)>0){
      fps[row,6] <- "TP" #The position is actually a TP
      final_result <- final_result[!(final_result$POS==candidate_pos & final_result$REF==candidate_ref & final_result$ALT==candidate_alt),] #remove the FN
    }
  }
}

final_result <- final_result[final_result$Status!="FP",] #Removing the previous FP annotation
final_result <- rbind(final_result, fps) #Adding the new annotation of FPs

# tp <- sum(final_result$Status=="TP") #204
# fp <- sum(final_result$Status=="FP") #4
# fn <- sum(final_result$Status=="FN") #13

output_file <- sub("filtering_norm_called.vcf.gz", "evaluation.tab", input_file)
write.table(final_result, output_file, row.names = F, quote = F)

print(paste0(output_file, " created"))
