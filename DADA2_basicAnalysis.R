# DADA2 - filtering, trimming, denoising, and merging 
# Based on DADA2 guide at http://benjjneb.github.io/dada2/dada-installation.html
# input: Trimmed fastq files

# load package
library(dada2); packageVersion("dada2")

# this will help identify output files later
study <- "run4" # CHANGE HERE TO MATCH OTHER DATA / SETUP

# get data - unzipped, trimmed data
path <- "./data"  # CHANGE HERE TO MATCH OTHER DATA / SETUP

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.trimmed.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.trimmed.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# visualize quality profiles
#plotQualityProfile(fnFs[3:4])
#plotQualityProfile(fnRs[3:4])

# filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=4) #  # CHANGE HERE TO MATCH OTHER DATA / SETUP

# learn error rates - NOTE: this takes a long time. especially if you have large samples
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# visualize error models
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# infer sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=4)
dadaRs <- dada(derepRs, err=errR, multithread=4)

# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# contruct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths (of ASVs, not reads!)
lengths <- table(nchar(getSequences(seqtab)))
write.csv(lengths, sprintf("%s_lenghtDistribuiton.csv", study))

# # remove chimeras - done in different script after merging runs
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=4, verbose=TRUE)
# dim(seqtab.nochim)
# # percentage chimeras
# sum(seqtab.nochim)/sum(seqtab)

# check how many reads made it through the different steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names
# export so you are sure you can check it later
write.csv2(track, sprintf("%s_SanityCheck.csv", study))

# CLASSIFICATION - done in different script after merging runs
# downloaded silva db following tutorial links (https://zenodo.org/record/824551#.Wh6SVktryL4)

# classify to species level
#taxa <- assignTaxonomy(seqtab.nochim, "~/2Big4Dropbox/DBs/SILVA_132_SSURef_Nr99_tax_silva_trunc_dada2.fa.gz", multithread=TRUE)  # CHANGE HERE TO MATCH OTHER DATA / SETUP
#taxa <- addSpecies(taxa, "~/2Big4Dropbox/DBs/silva_species_assignment_v132.fa.gz")  # CHANGE HERE TO MATCH OTHER DATA / SETUP
# taxa <- assignTaxonomy(seqtab.nochim, "~/DBs/SILVA_132_SSURef_Nr99_tax_silva_trunc_dada2.fa.gz", multithread=4)  # CHANGE HERE TO MATCH OTHER DATA / SETUP
# taxa <- addSpecies(taxa, "~/DBs/silva_species_assignment_v132.fa.gz")  # CHANGE HERE TO MATCH OTHER DATA / SETUP

# evaluate accuracy
# can only be done if you have mock community

# save R workspace for look at error models + general de-bugging
save.image(file=sprintf("%s_workspace.RData", study))