# this script is based on http://benjjneb.github.io/dada2/bigdata.html
# it merges data from several sequencing runs - that data has been through:
# DADA2_BasicAnalysis.R, and you should now have a "studyName_workspace.RData" for each run.
# if you have more than five runs, this script can just be expanded.
# This script also checks for chimeras and classifies your ASVs based on a given database.
# the studyName to continue with is defined in the varibale finalName. you can use this in 
# DADA2_to_phyloseq.R and DADA2_filterAndNorm.R

# input: 
#   Workspace data from DADA2_BasicAnalysis.R
#   Databases for classification


finalName <- "TSPfinal"  # CHANGE HERE TO MATCH OTHER DATA / SETUP
study1 <- "run0"  # CHANGE HERE TO MATCH OTHER DATA / SETUP
study2 <- "run1"  # CHANGE HERE TO MATCH OTHER DATA / SETUP
study3 <- "run2"  # CHANGE HERE TO MATCH OTHER DATA / SETUP
study4 <- "run3"  # CHANGE HERE TO MATCH OTHER DATA / SETUP
study5 <- "run4"  # CHANGE HERE TO MATCH OTHER DATA / SETUP

library(dada2); packageVersion("dada2")

# data1
load(sprintf("%s_workspace.RData", study1))
saveRDS(seqtab, sprintf("%s_seqtab.rds", study1))

# data2
load(sprintf("%s_workspace.RData", study2))
saveRDS(seqtab, sprintf("%s_seqtab.rds", study2))

# data3
load(sprintf("%s_workspace.RData", study3))
saveRDS(seqtab, sprintf("%s_seqtab.rds", study3))

# data4
load(sprintf("%s_workspace.RData", study4))
saveRDS(seqtab, sprintf("%s_seqtab.rds", study4))

# data5
load(sprintf("%s_workspace.RData", study5))
saveRDS(seqtab, sprintf("%s_seqtab.rds", study5))

# Merge multiple runs (if necessary)
data1 <- readRDS(sprintf("%s_seqtab.rds", study1))
data2 <- readRDS(sprintf("%s_seqtab.rds", study2))
data3 <- readRDS(sprintf("%s_seqtab.rds", study3))
data4 <- readRDS(sprintf("%s_seqtab.rds", study4))
data5 <- readRDS(sprintf("%s_seqtab.rds", study5))

st.all <- mergeSequenceTables(data1, data2, data3, data4, data5)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# check how many reads made it through the chimera finder
track <- cbind(rowSums(st.all), rowSums(seqtab.nochim))
colnames(track) <- c("tabled", "nonchim")
# export so you are sure you can check it later
write.csv2(track, "ChimeraSurvival.csv")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/DBs/SILVA_132_SSURef_Nr99_tax_silva_trunc_dada2.fa.gz", multithread=4)  # CHANGE HERE TO MATCH OTHER DATA / SETUP
taxa <- addSpecies(taxa, "~/DBs/silva_species_assignment_v132.fa.gz")  # CHANGE HERE TO MATCH OTHER DATA / SETUP

# save workspace
save.image(file= sprintf("%s_workspace.RData", finalName))

