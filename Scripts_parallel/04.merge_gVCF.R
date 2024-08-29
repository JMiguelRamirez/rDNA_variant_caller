#Code to correlate allele frequencies between the true and the obtained with our pipeline

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

args = commandArgs(trailingOnly=TRUE)

list <- args[1] #This is a .txt file with the name of all samples that we want to merge
# list <- read.table("WGS_ids.tab")
# list <- read.table("WGS_GTEx_ids.tab")
output <- args[2] #This is the output directory
# output <- "calling_ploidy_20"
# output <- "calling_GTEx"
list <- read.table(list)
# list <- read.table("calling_ploidy_20/SRR1997411.20_filtering_norm_called.vcf.gz")

#Initialize the final vcf with the information of the first sample
first_sample <- list$V1[1]
final_vcf <- read.table(paste0(output, "/", first_sample, ".20_filtering_norm_called.vcf.gz"))[,c(1,2,4,5,10)]
colnames(final_vcf) <- c("CHR", "POS", "REF", "ALT", first_sample)
final_vcf[[paste0(first_sample, ".GT")]] <- sapply(final_vcf[,5], function(entry) strsplit(entry, ":")[[1]][1])
final_vcf[[paste0(first_sample, ".AD")]] <- sapply(final_vcf[,5], function(entry) strsplit(entry, ":")[[1]][2])
final_vcf <- final_vcf[,-5]

#Iterate over the last samples to add their info into the final vcf
for(sample in list$V1[-1]){
  print(grep(sample, list$V1))
  vcf <-  read.table(paste0(output, "/", sample, ".20_filtering_norm_called.vcf.gz"))[,c(1,2,4,5,10)]
  colnames(vcf) <- c("CHR", "POS", "REF", "ALT", sample)
  vcf[[paste0(sample, ".GT")]] <- sapply(vcf[,5], function(entry) strsplit(entry, ":")[[1]][1])
  vcf[[paste0(sample, ".AD")]] <- sapply(vcf[,5], function(entry) strsplit(entry, ":")[[1]][2])
  vcf <- vcf[,-5]
  
  #Merge the two vcf files
  final_vcf <- merge(final_vcf, vcf, by=c("CHR", "POS", "REF", "ALT"), all = T) 
  #The parameter all=T is for adding variants that might not be in all samples. At the end, we will add 0s for the rest of the samples
}
#Some samples might have NAs in GT and AD columns. I will not replace those with 0s:
# final_vcf[,grep(".GT", colnames(final_vcf), fixed=T)][is.na(final_vcf[,grep(".GT", colnames(final_vcf), fixed=T)])] <- "0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0"
# final_vcf[,grep(".AD", colnames(final_vcf), fixed=T)][is.na(final_vcf[,grep(".AD", colnames(final_vcf), fixed=T)])] <- "0,0"


#Save final vcf (even though it is no longer a vcf format, so I will save it as txt)

write.table(final_vcf, paste0(output, "/merged.txt"), row.names = F, quote = F, sep="\t")
