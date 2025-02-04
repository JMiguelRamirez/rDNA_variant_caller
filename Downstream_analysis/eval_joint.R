#Code to assing FP, TP or FN to the variant calling

## libraries 
library(data.table) #More efficient data handling
library(dplyr)

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Read vcf file from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
# input_file <- "Simulations_10/var/jobs/Joint_called_all/Simulated.20_filtering_norm.vcf.gz"
# input_file <- "Simulations_3/var/jobs/Joint_called_all/Simulated.20_filtering_norm.vcf.gz"
# input_file <- "Simulations_13/var/jobs/Joint_called_all/Simulated.20_filtering_norm.vcf.gz"
# input_file <- "Simulations_29/var/jobs/Joint_called_all/Simulated.5_filtering_norm.vcf.gz"
vcf <- fread(input_file)
vcf <- vcf[vcf$FILTER=="PASS",]

simulation <- sub("/.*", "", input_file)

#Get allele frequency
vcf_t <- vcf %>%
  mutate(
    across(
      .cols = 10:ncol(vcf), # Select columns from the 10th column onward
      ~ {
        # Split the column into components
        parts <- strsplit(as.character(.), "[:,]")
        value1 <- sapply(parts, function(x) as.numeric(x[3]))
        value2 <- sapply(parts, function(x) as.numeric(x[4]))
        value1 / value2
      },
      .names = "AF_{col}" # Create new columns for AF values
    )
  ) %>%
  select(-c("#CHROM", ID, QUAL, INFO, FORMAT))

columns <- c("POS", "REF", "ALT", grep("AF", colnames(vcf_t), value = T))
vcf_t <- vcf_t[, ..columns]

#If 0 read depth, we get NaN, I will replace with 0s:
vcf_t[is.na(vcf_t)] <- 0

#Creating one with an allele frequency filter
vcf_af <- vcf_t
vcf_af[vcf_af<0.05] <- 0

#Get the true variants for the sample and simulation:
dir <- paste0(simulation, "/var/jobs/Catalogue/")
### Get true_variants ######
true_variants <- read.delim(paste0(simulation, "/var/jobs/Catalogue/true_variants.tbl"), sep=" ")


#Annotate variants as TP, FP or FN:
# Create variant IDs 
vcf_t <- vcf_t %>%
  mutate(variant = paste(POS, REF, ALT, sep = "_"))
vcf_af <- vcf_af %>%
  mutate(variant = paste(POS, REF, ALT, sep = "_"))

true_variants <- true_variants %>%
  mutate(variant = paste(POS, REF, ALT, sep = "_"))


#Annotate existing variants
classify_variant <- function(variants, afs, sample) {
  true_vars <- true_variants$variant[true_variants$Sample == sample]
  output <- c()
  for(row in 1:length(variants)){
    variant <- variants[row]
    af <- afs[row]
    if (af == 0 && !(variant %in% true_vars)) {
      output <- c(output, "TN")
    } else if (af > 0 && variant %in% true_vars) {
      output <- c(output, "TP")
    } else if (af > 0 && !(variant %in% true_vars)) {
      output <- c(output, "FP")
    } else if (af == 0 && variant %in% true_vars) {
      output <- c(output, "FN")
    } else {
      NA  # Handle unexpected cases
    }
  }
  
  return(output)
}

# Loop through each sample and classify
sample_names <- grep("AF_", colnames(vcf_t), value = T)
for (sample in sample_names) {
  af_col <- sample
  class_col <- sub("AF_", "Class_", af_col)
  
  vcf_t[[class_col]] <- mapply(classify_variant, vcf_t$variant, vcf_t[[af_col]], sample = sub("AF_", "", sample))
  vcf_af[[class_col]] <- mapply(classify_variant, vcf_af$variant, vcf_af[[af_col]], sample = sub("AF_", "", sample))
}

columns <- c("POS", "REF", "ALT", "variant", grep("Class", colnames(vcf_t), value = T))
classification <- vcf_t[, ..columns]
classification_af <- vcf_af[, ..columns]

#Add FNs that are not even in vcf_t
true_variants_fn <- true_variants[,c("POS", "REF", "ALT", "variant")]
true_variants_fn <- true_variants_fn[!duplicated(true_variants_fn),] #Remoce duplicated rows

fns <- true_variants_fn %>% #Get variants that are FN in all individuals
  filter(!variant %in% classification$variant) 

fns <- cbind(fns, matrix("FN", nrow = nrow(fns), ncol = length(5:ncol(classification)))) #Add "FN" for the size of "classification"

evaluations <- rbind(classification, fns, use.names=F)  #Add the FNs into classidication
evaluations_af <- rbind(classification_af, fns, use.names=F)  #Add the FNs into classidication


#Add whether the variants in are in what we call repetitive or non-repetitive
bed <- fread("data/pre-rRNA_47S.included_intervals.bed")[,c(1:4)]
colnames(bed) <- c("POS", "START", "END", "REGION")

# Check if positions in the vcf fall within any of the non-repetitive regions in the BED file
evaluations <- evaluations %>%
  mutate(
    Callable = ifelse(
      mapply(function(pos) {
        # Check if the position is within the non-repetitive regions (defined in BED)
        any(pos >= bed$START & pos <= bed$END)
      }, evaluations$POS),
      "callable",
      "blacklisted"
    )
  )

evaluations_af <- evaluations_af %>%
  mutate(
    Callable = ifelse(
      mapply(function(pos) {
        # Check if the position is within the non-repetitive regions (defined in BED)
        any(pos >= bed$START & pos <= bed$END)
      }, evaluations_af$POS),
      "callable",
      "blacklisted"
    )
  )

# sum(colSums(evaluations=="TP") ) #25216
# sum(colSums(evaluations=="FP") ) #9165
# sum(colSums(evaluations=="FN") ) #6355


#Fixing INDELs that might be differentially annotated:
evaluations <- as.data.frame(evaluations)
for(variant_index in 1:nrow(evaluations)){
  cols <- c(5:(ncol(evaluations)-1))
  if(variant_index>nrow(evaluations)){ #This means that we are finished, it is just that some rows with FNs were removed :)
    break
  }
  if(sum(evaluations[variant_index, cols]=="FP")==length(cols)){
    #A position that is always FP, could it be that the way we call indels is different to the software?
    ref <- as.character(evaluations[variant_index,2])
    alt <- as.character(evaluations[variant_index,3])
    pos <- as.numeric(evaluations[variant_index,1])
    
    #Check insertions: 
    #GATK and Lofreq might get "confused" with some insertions if the nucleotides after the insertion might match the first nucleotides that are inserted. 
    #This will only affects the simulated data, with real data we will call a variant that is real, the only thing that changes is in which position we call it.
    #I will fix this as they are named as FP but they are TPs.
    #In insertions we can just ignore the first nucleotide of the insertion (as we manually put it ther because it is was usually gatk does), and then
    #depending on the number of positions shifted, we append the nucleotides shifted at the end, and then we add the nucleotide that is now the first, as the end as well
    for(shifted in 1:5){ #I added as maximum a shift of 5, even though the maximum shift that happened in an indel positions is 3, but we iterate until 5 just in case
      candidate_alt <- paste0(paste0(strsplit(alt, "")[[1]][-c(1:shifted)], collapse = ""), paste0(strsplit(alt, "")[[1]][c(1:shifted+1)], collapse = "")) #Note that the way they would be added are the inverse, so not 1 and 2, but 2 and 1
      candidate_pos <- pos + shifted
      if(sum(true_variants$POS==candidate_pos & true_variants$ALT==candidate_alt)>0){
        
        true_variant <- paste0(candidate_pos, "_", ref, "_", candidate_alt)
        truth <- true_variants[true_variants$variant==true_variant,]
        samples_with_variant <- sample_names %in% truth$Sample #Is the variant present in the different samples?
        samples_with_variant[samples_with_variant==T] <- "TP"
        samples_with_variant[samples_with_variant==F] <- "FP"
        evaluations[variant_index, cols] <- samples_with_variant #Correcting some FPs for TPs
        evaluations <- evaluations[!(evaluations$POS==candidate_pos & evaluations$ALT==candidate_alt),] #remove the FN
        
      }
    }
    #Now, we can check the deletions:
    if(nchar(ref)<2){next} #In this case, for sure it is not a deletion, so we can skip it to avoid generating NAs 
    #I will create a fancy iteration such as in insertions, because I know we need it
    
    #For deletions we do something similar as with insertions but with the reference allele. We remove the first one, and we will always add at the end the same as the new first. Then we can have different shiftings
    for(shifted in 1:5){ #I added as maximum a shift of 5, even though the maximum shift that happened in an indel positions is 3, but we iterate until 5 just in case
      candidate_ref <- paste0(paste0(strsplit(ref, "")[[1]][-c(1:shifted)], collapse = ""), paste0(strsplit(ref, "")[[1]][c(1:shifted+1)], collapse = "")) #An example of deletions shifted once: instead of GTCCC > G we get TCCCT > T. This would happen if int he reference after the insertion we had a T, so this deletion reporting would be preferred by GATK
      candidate_alt <- strsplit(candidate_ref, "")[[1]][1]
      candidate_pos <- pos + shifted
      if(sum(true_variants$POS==candidate_pos & true_variants$REF==candidate_ref & true_variants$ALT==candidate_alt)>0){
        #This positions is a FP for some samples:
        true_variant <- paste0(candidate_pos, "_", candidate_ref, "_", candidate_alt)
        truth <- true_variants[true_variants$variant==true_variant,]
        samples_with_variant <- sample_names %in% truth$Sample #Is the variant present in the different samples?
        samples_with_variant[samples_with_variant==T] <- "TP"
        samples_with_variant[samples_with_variant==F] <- "FP"
        evaluations[variant_index, cols] <- samples_with_variant #Correcting some FPs for TPs
        evaluations <- evaluations[!(evaluations$POS==candidate_pos & evaluations$REF==candidate_ref & evaluations$ALT==candidate_alt),] #remove the FN
      }
    }
  }
}

#The same for the table with af filtered:
evaluations_af <- as.data.frame(evaluations_af)
for(variant_index in 1:nrow(evaluations_af)){
  cols <- c(5:(ncol(evaluations_af)-1))
  if(variant_index>nrow(evaluations_af)){ #This means that we are finished, it is just that some rows with FNs were removed :)
    break
  }
  if(sum(evaluations_af[variant_index, cols]=="FP")==length(cols)){
    #A position that is always FP, could it be that the way we call indels is different to the software?
    ref <- as.character(evaluations_af[variant_index,2])
    alt <- as.character(evaluations_af[variant_index,3])
    pos <- as.numeric(evaluations_af[variant_index,1])
    
    #Check insertions: 
    #GATK and Lofreq might get "confused" with some insertions if the nucleotides after the insertion might match the first nucleotides that are inserted. 
    #This will only affects the simulated data, with real data we will call a variant that is real, the only thing that changes is in which position we call it.
    #I will fix this as they are named as FP but they are TPs.
    #In insertions we can just ignore the first nucleotide of the insertion (as we manually put it ther because it is was usually gatk does), and then
    #depending on the number of positions shifted, we append the nucleotides shifted at the end, and then we add the nucleotide that is now the first, as the end as well
    for(shifted in 1:5){ #I added as maximum a shift of 5, even though the maximum shift that happened in an indel positions is 3, but we iterate until 5 just in case
      candidate_alt <- paste0(paste0(strsplit(alt, "")[[1]][-c(1:shifted)], collapse = ""), paste0(strsplit(alt, "")[[1]][c(1:shifted+1)], collapse = "")) #Note that the way they would be added are the inverse, so not 1 and 2, but 2 and 1
      candidate_pos <- pos + shifted
      if(sum(true_variants$POS==candidate_pos & true_variants$ALT==candidate_alt)>0){
        
        true_variant <- paste0(candidate_pos, "_", ref, "_", candidate_alt)
        truth <- true_variants[true_variants$variant==true_variant,]
        samples_with_variant <- sample_names %in% truth$Sample #Is the variant present in the different samples?
        samples_with_variant[samples_with_variant==T] <- "TP"
        samples_with_variant[samples_with_variant==F] <- "FP"
        evaluations_af[variant_index, cols] <- samples_with_variant #Correcting some FPs for TPs
        evaluations_af <- evaluations_af[!(evaluations_af$POS==candidate_pos & evaluations_af$ALT==candidate_alt),] #remove the FN
        
      }
    }
    #Now, we can check the deletions:
    if(nchar(ref)<2){next} #In this case, for sure it is not a deletion, so we can skip it to avoid generating NAs 
    #I will create a fancy iteration such as in insertions, because I know we need it
    
    #For deletions we do something similar as with insertions but with the reference allele. We remove the first one, and we will always add at the end the same as the new first. Then we can have different shiftings
    for(shifted in 1:5){ #I added as maximum a shift of 5, even though the maximum shift that happened in an indel positions is 3, but we iterate until 5 just in case
      candidate_ref <- paste0(paste0(strsplit(ref, "")[[1]][-c(1:shifted)], collapse = ""), paste0(strsplit(ref, "")[[1]][c(1:shifted+1)], collapse = "")) #An example of deletions shifted once: instead of GTCCC > G we get TCCCT > T. This would happen if int he reference after the insertion we had a T, so this deletion reporting would be preferred by GATK
      candidate_alt <- strsplit(candidate_ref, "")[[1]][1]
      candidate_pos <- pos + shifted
      if(sum(true_variants$POS==candidate_pos & true_variants$REF==candidate_ref & true_variants$ALT==candidate_alt)>0){
        #This positions is a FP for some samples:
        true_variant <- paste0(candidate_pos, "_", candidate_ref, "_", candidate_alt)
        truth <- true_variants[true_variants$variant==true_variant,]
        samples_with_variant <- sample_names %in% truth$Sample #Is the variant present in the different samples?
        samples_with_variant[samples_with_variant==T] <- "TP"
        samples_with_variant[samples_with_variant==F] <- "FP"
        evaluations_af[variant_index, cols] <- samples_with_variant #Correcting some FPs for TPs
        evaluations_af <- evaluations_af[!(evaluations_af$POS==candidate_pos & evaluations_af$REF==candidate_ref & evaluations_af$ALT==candidate_alt),] #remove the FN
      }
    }
  }
}


# sum(colSums(evaluations=="TP") ) #25984
# sum(colSums(evaluations=="FP") ) #8397
# sum(colSums(evaluations=="FN") ) #5055

output_file <- sub("_filtering_norm.vcf.gz", "_evaluation.tab", input_file)
write.table(evaluations, output_file, row.names = F, quote = F)
output_file <- sub("_filtering_norm.vcf.gz", "_evaluation_af.tab", input_file)
write.table(evaluations_af, output_file, row.names = F, quote = F)

print(paste0(output_file, " created"))
