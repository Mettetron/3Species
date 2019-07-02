# Extract data from DADA2 run into useful tables. Add meta data (=info about samples)
# Filter ASVs to only include Bacteria and have minumum length of 400 bp
# Deal with some classification issues such as replacing "unclassified" with something more useful

# input:
#   Workspace data from DADA2_BasicAnalysis.R or DADA2_CombineRuns_ChimTax.R
#   a csv of metaData, with samples in rows and attributes in columns. 

# set your working directory and studyName
setwd("~/ /") # CHANGE HERE TO MATCH OTHER DATA / SETUP
study <- "TSPfinal" # CHANGE HERE TO MATCH OTHER DATA / SETUP

# install phyloseq - you only have to do this once. 
# look at: https://joey711.github.io/phyloseq/install.htmlif you run in to problems
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

# load phyloseq
library(phyloseq); packageVersion("phyloseq")

# load workspace data from DADA2_BasicAnalysis.R or DADA2_CombineRuns_ChimTax.R
load(sprintf("%s_workspace.RData", study))

# reset working directory and study name
setwd("~/ /") # CHANGE HERE TO MATCH OTHER DATA / SETUP
  # CHANGE HERE TO MATCH OTHER DATA / SETUP
study <- "TSPfinal"  # CHANGE HERE TO MATCH OTHER DATA / SETUP

# import sample data
# note: change to read.csv() if needed
samdf <- read.csv2(sprintf("%s_metaData.csv", study))  # CHANGE HERE TO MATCH OTHER DATA / SETUP
rownames(samdf) <- samdf[, 1]
samdf[, 1] <- NULL

# make phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  
               sample_data(samdf), 
               tax_table(taxa))

# These are different things you can do with phyloseq, see e.g.: https://joey711.github.io/phyloseq/tutorials-index.html
  # plot_richness(ps, x="state", measures=c("Shannon", "Simpson"), color="dead_or_alive")
  # ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
  # plot_ordination(ps, ord.nmds.bray, color="species", title="Bray NMDS")
  # top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
  # ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  # ps.top20 <- prune_taxa(top20, ps.top20)
  # plot_bar(ps.top20, x="species", fill="Family") + facet_wrap(~When, scales="free_x")

# get dataframes from phyloseq
data.asv <- as.data.frame(otu_table(ps))  # this one takes a while for big data
data.tax <- as.data.frame(tax_table(ps))
data.meta <- as.data.frame(sample_data(ps))

### remove ASVs with short sequences and ASVs not classified to Bacteria
nchar(as.character(rownames(data.tax)))
data.tax.longseq <- subset(data.tax, nchar(as.character(rownames(data.tax))) >= 400)
data.tax.bact <- data.tax.longseq[data.tax.longseq$Kingdom == "Bacteria", ]
data.tax <- data.tax.bact

data.asv <- data.asv[, colnames(data.asv) %in% rownames(data.tax)]

# replace "NA" with "[higherTaxonomicLevel]"
data.tax <- data.tax[!is.na(data.tax$Kingdom), ] # remove rows where kingdom is NA
data.asv <- data.asv[, colnames(data.asv) %in% rownames(data.tax)]

data.tax <- as.matrix(data.tax)
# manual unUnclassify
for (i in seq(1:nrow(data.tax))) {
  for (j in seq(1:ncol(data.tax))) {
    if (!is.na(data.tax[i, j])){
      mem.tax <- data.tax[i, j]
    } else {
      data.tax[i, j] <- mem.tax
    }
  } 
}

data.tax <- as.data.frame(data.tax)

# replace "uncultured" with "[higherTaxonomicLevel]"
data.tax <- data.tax[!is.na(data.tax$Kingdom), ] # remove rows where kingdom is NA
data.asv <- data.asv[, colnames(data.asv) %in% rownames(data.tax)]

data.tax <- as.matrix(data.tax)
# manual unUncultured
for (i in seq(1:nrow(data.tax))) {
  for (j in seq(1:ncol(data.tax))) {
    if (!data.tax[i, j] == "uncultured"){
      mem.tax <- data.tax[i, j]
    } else {
      data.tax[i, j] <- mem.tax
    }
  } 
}

data.tax <- as.data.frame(data.tax)

# replace "Ac37b" with "[higherTaxonomicLevel]"
data.tax <- data.tax[!is.na(data.tax$Kingdom), ] # remove rows where kingdom is NA
data.asv <- data.asv[, colnames(data.asv) %in% rownames(data.tax)]

data.tax <- as.matrix(data.tax)
# manual unUncultured
for (i in seq(1:nrow(data.tax))) {
  for (j in seq(1:ncol(data.tax))) {
    if (!data.tax[i, j] == "Ac37b"){
      mem.tax <- data.tax[i, j]
    } else {
      data.tax[i, j] <- mem.tax
    }
  } 
}

data.tax <- as.data.frame(data.tax)
# make ASV overview df
colnames(data.asv) == rownames(data.tax)  # check that things line up
numbers <- paste("ASV_", seq(1:length(colnames(data.asv))), sep="")  # make list of ASV_1 numbers
ASVs <- cbind(numbers, colnames(data.asv), data.tax)  # make df with seqs + tax
rownames(ASVs) <- ASVs$numbers
colnames(ASVs)[2] <- "seq"

# exchange seqs for numbers in data tax and data asv
# this makes the dataframes much faster to import, and makes them more readable in excel
colnames(data.asv) <- numbers
rownames(data.tax) <- numbers


# save sampleID as readID
data.meta$readID <- rownames(data.meta)
sum(rownames(data.meta) == rownames(data.asv)) == length(rownames(data.meta))
# make sampleID == tubeID
rownames(data.meta) <- data.meta$tube_ID
rownames(data.asv) <- data.meta$tube_ID

# trim metadata
data.meta$tube_ID <- NULL
data.meta$SRA_BioSample <- NULL
data.meta$SRA_Experiment <- NULL
data.meta$SRA_Sample <- NULL

# export dataframes
write.csv2(ASVs, sprintf("%s_DADA2_seqTAXtable.csv", study))
write.csv2(data.asv, sprintf("%s_DADA2_ASVtable.csv", study))
write.csv2(data.tax, sprintf("%s_DADA2_TAXtable.csv", study))
write.csv2(data.meta, sprintf("%s_DADA2_SAMPLEtable.csv", study))

# also make tax.asv table for easy overview
data.asv.t <- as.data.frame(t(data.asv))
colnames(data.asv.t) == (rownames(data.meta))
colnames(data.asv.t) <- paste(colnames(data.asv.t), data.meta$nest)
data.tax.asv <- cbind(data.tax, data.asv.t)
data.tax.asv$seq <- ASVs$seq
write.csv2(data.tax.asv,  sprintf("%s_DADA2_TAX_and_ASVtable.csv", study))
