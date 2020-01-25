# Decontaminate and normalize amplicon data

# Things you need to change:
  # working directory - to match with your own computer
  # study name -  must be the same as used in the two previous scripts
  # name of negative controls in the DECONTAMINATION section
  # minimum reads (=sampling depths) in the SUBSAMPLING section
      # same minimum read number is also used in FRACTION OF TOTAL READS section

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns
#         Negative samples must be included and have "neg" in all columns

# set your working directory and studyName
setwd("~/ /") # CHANGE HERE TO MATCH OTHER DATA / SETUP
study <- "TSPfinal" # CHANGE HERE TO MATCH OTHER DATA / SETUP

# import raw data
data.asv <- read.csv2(sprintf("%s_DADA2_ASVtable.csv", study), row.names=1)
data.tax <- read.csv2(sprintf("%s_DADA2_TAXtable.csv", study), row.names=1)
data.meta <- read.csv2(sprintf("%s_DADA2_SAMPLEtable.csv", study), row.names=1)

# checking that everything matches up
sum(rownames(data.meta) == rownames(data.asv)) == length(rownames(data.meta))
sum(rownames(data.tax) == colnames(data.asv)) == length(rownames(data.tax))


################### DECONTAMINATION  ###########################################
neg_control_samples <- rownames(data.meta[which(data.meta$nest == "neg"), ]) # a list of sample IDs for negative control samples

# how many more times can a ASV be present in a real sample than in neg ctrl, and still be suspect.
# e.g. 0.2: has to be 5x more in real sample to be OK
# adjust this based on you data, and on how stringent you want to be
multiplication_factor <- 0.1

check_negative_controls <- function(data.asv, neg_control_samples, multiplication_factor){
  sample_names_no_negatives <- row.names(data.asv)[row.names(data.asv) %in% neg_control_samples == FALSE]
  suspect_ASVs = c()
  for(i in 1:dim(data.asv)[2]){
    sample_read_counts <- data.asv[sample_names_no_negatives,i]
    negative_read_counts <- data.asv[neg_control_samples,i]
    if(max(negative_read_counts) > multiplication_factor*max(sample_read_counts)){
      suspect_ASVs[length(suspect_ASVs)+1] <- colnames(data.asv)[i]
    }
  }
  return(suspect_ASVs)
}

suspect_ASVs <- check_negative_controls(data.asv, neg_control_samples, multiplication_factor)
suspect_ASVs

# consider the suspect ASVs and remove from your data if you want to
data.asv.contaminants <-  data.asv[, colnames(data.asv) %in% suspect_ASVs]
data.tax.contaminants <- data.tax[rownames(data.tax) %in% suspect_ASVs, ]
data.asv <- data.asv[, !colnames(data.asv) %in% suspect_ASVs]
data.tax <- data.tax[!rownames(data.tax) %in% suspect_ASVs, ]

# remove neg.ctrl sample from data
data.meta <- data.meta[!rownames(data.meta) %in% neg_control_samples, ]
data.asv <- data.asv[!rownames(data.asv) %in% neg_control_samples, ]

################### SUBSAMPLING ################################################

# method A - subsampling
# this method randomly samples x reads from each sampleID with x+ reads

library(vegan)

# set a minimun number of reads e.g. 1st Quartile of the asv data summed per sample
rowSums(data.asv)
#min_reads <- summary(rowSums(data.asv))[2]
min_reads <- 3000 # CHANGE HERE TO MATCH OTHER DATA / SETUP
# remove samples-ids below min_reads
data.asv.norm <- data.asv[rowSums(data.asv)>min_reads, ]
# randomly sample the same number of reads from all the remaining sample_ids
set.seed(42)
#min_reads <- 100
data.asv.norm <- as.data.frame(rrarefy(data.asv.norm, min_reads))
# remove asvs which are no longer represented in the data
data.asv.norm <- data.asv.norm[, colSums(data.asv.norm)>0]

# make table of sample IDs that were lost in subsampling
lost.sampleIDs <- data.meta[!(data.meta$sampleID %in% rownames(data.asv.norm)), ]
lost.sampleIDs$reads <- rowSums(data.asv[rownames(lost.sampleIDs), ])

# remove samples from meta data and tax data
data.meta.norm <- data.meta[(rownames(data.meta) %in% rownames(data.asv.norm)), ]
data.tax.norm <- data.tax[rownames(data.tax) %in% colnames(data.asv.norm), ]

# make tax-asv overview table
data.asv.norm.t <- as.data.frame(t(data.asv.norm))
data.tax.asv.norm <- cbind(data.tax.norm, data.asv.norm.t)
sample_nest_names <- paste(data.meta.norm$colony, data.meta.norm$sampleID)
colnames(data.tax.asv.norm)[8:ncol(data.tax.asv.norm)] <- sample_nest_names

# export
write.csv(data.asv.norm.t, file=sprintf("%s_DADA2_Subsampl%s.ASVtable.csv", study, min_reads))
write.csv(data.tax.norm, file=sprintf("%s_DADA2_Subsampl%s.TAXtable.csv", study, min_reads))
write.csv(data.meta.norm, file=sprintf("%s_DADA2_Subsampl%s.SAMPLEtable.csv", study, min_reads))
write.csv2(data.tax.asv.norm, file=sprintf("%s_DADA2_Subsampl%s.TAXandASVtable.csv",study, min_reads))


################### FRACTION OF TOTAL READS #####################################

# set a minimun number of reads e.g. 1st Quartile of the asv data summed per sample
# here we just use same miminum as above
# transform asv table
data.asv.t <- as.data.frame(t(data.asv))
# discard samples with less than your minimum of reads reads
data.asv.t.minreads <- data.asv.t[, which(colSums(data.asv.t) >= min_reads)]
data.meta.norm <- data.meta[which(rownames(data.meta) %in% colnames(data.asv.t.minreads)), ]

# fractionify
data.asv.t.norm <- apply(data.asv.t.minreads, 2, function(x) {x / sum(x)})

# just checking that everything sums to 1
colSums(data.asv.t.norm)

# make tax-asv overview table
rownames(data.tax) == rownames(data.asv.t.norm)
data.tax.asv.norm <- cbind(data.tax, round(data.asv.t.norm,4))
sample_nest_names <- paste(data.meta.norm$nest, data.meta.norm$sampleID)
colnames(data.tax.asv.norm)[8:ncol(data.tax.asv.norm)] <- sample_nest_names

# export
write.csv(data.asv.t.norm, sprintf("%s_DADA2_Fraction.ASVtable.csv", study))
write.csv(data.tax, sprintf("%s_DADA2_Fraction.TAXtable.csv", study))
write.csv(data.meta.norm, sprintf("%s_DADA2_Fraction.SAMPLEtable.csv", study))
write.csv(data.tax.asv.norm, sprintf("%s_DADA2_Fraction.TAXandASVtable.csv", study))

