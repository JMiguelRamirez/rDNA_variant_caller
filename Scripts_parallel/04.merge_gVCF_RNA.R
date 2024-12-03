#Code to correlate allele frequencies between the true and the obtained with our pipeline

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

args = commandArgs(trailingOnly=TRUE)

list <- args[1] #This is a .txt file with the name of all samples that we want to merge
# list <- "WGS_ids.tab"
output <- args[2] #This is the output directory
# output <- "Results_calling"
list <- read.table(list)


suppressWarnings(library(data.table))     # For efficient data handling
suppressPackageStartupMessages(library(dplyr)) #To split columns
suppressWarnings(library(tidyr)) #to split columns
suppressWarnings(library(parallel)) #To run in parallel

process_sample <- function(sample){
  tryCatch({ #I am running this to notify if some sample failed but still generating the output for the rest
    #Reading vcf file
    path <- paste0(output, "/", sample, ".vcf.gz")
    if(file.size(path)==0){ #If the size file is 0, we return NULL
      return(NULL)
    }
    # vcf <- fread(path,select = c(1, 2, 4, 5, 10))
    vcf <- read.table(path)[,c(1,2,4,5,10)]
    #Procesing vcf file
    colnames(vcf) <- c("CHR", "POS", "REF", "ALT", sample)
    suppressWarnings(vcf <- vcf %>% #The warnings here refer to ././././, which DE is set as NA, this is what we want, so no problem. This happens very rarely
                       separate(sample, into=c("GT", "AD", "DP", "GQ"), sep=":"))
    suppressWarnings(vcf <- vcf %>% #The warnings are expected (and normal), it means that there are no reads for ALT_D, so we put the NA
                       separate(AD, into=c("REF_D", "ALT_D"), sep=","))
    
    vcf <- vcf[,c(2:8)]
    colnames(vcf)[4:7] <- paste0(sample, ".", colnames(vcf)[4:7])
    return(vcf)
  }, error = function(e) {
    # If an error occurs, print the sample ID that failed
    cat("Failed to read data for", grep(sample, list$V1), "\n")
  })

}
samples <- list$V1
print("about to start")
# result_list <- lapply(samples, process_sample)
result_list <- mclapply(samples, process_sample, mc.cores = 112)
# Remove NULLs in case any files didn't exist
print("over here")
result_list <- result_list[!sapply(result_list, is.null)]
print("about to reduce")
# Combine all data frames, including the first sample
final_vcf <- Reduce(function(x, y) merge(x, y, by = c("POS", "REF", "ALT"), all = TRUE), result_list)

#Save final vcf (even though it is no longer a vcf format, so I will save it as txt)
print(paste0(output, "/merged_with_DP.txt"))
write.table(final_vcf, paste0(output, "/merged_with_DP.txt"), row.names = F, quote = F, sep="\t")
