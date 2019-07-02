# Make an NMDS ordination plot based on ASV data. 
# Here I color with specific colors by host species, and the shapes denote different host populations.
# note: Each run of the metaMDS function will give slightly different outputs, as this is an iterative method

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#         it is recommended to use subsampled data for calculation of distance metrics
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns


setwd("~/ /")
study <- "TSPfinal"

library(vegan)

# import subsampled data
data.asv <- read.csv(sprintf("%s_DADA2_Subsampl3000.ASVtable.csv", study), row.names=1)
data.tax <- read.csv(sprintf("%s_DADA2_Subsampl3000.TAXtable.csv", study), row.names=1)
data.meta <- read.csv(sprintf("%s_DADA2_Subsampl3000.SAMPLEtable.csv", study), row.names=1)

# make pop symbols
pch.df <- as.data.frame(table(factor(data.meta$population)))
pch.df$symbol <- c(c(1:nrow(pch.df)), rep(NA, (nrow(pch.df)-nrow(pch.df))))
pch.df$symbol <- c(1,2,3,4,5,6,7,8,9,10,15,12,13,14,11,16) # order modified so that shared populations have filled symbols
pop.vector <- unlist(lapply(data.meta$population, function(x) pch.df$symbol[pch.df$Var1 %in% x]))

# make species colors
col.df <- data.frame("species" = c("dumicola", "mimosarum", "sarasinorum"), 
                     "color" = c('#d95f02', '#7570b3', '#1b9e77'))
species.vector <- unlist(lapply(data.meta$species, function(x) col.df$color[col.df$species %in% x]))

### NMDS
nmds.asv.t <- metaMDS(as.data.frame(t(data.asv)), k=2, distance = "bray", binary=FALSE) # from vegan package

# PLOT and EXPORT
jpeg(file="NMDS.jpg", width = 23, height = 23, units="cm", res=600) # this prints
plot(nmds.asv.t, display = "sites", type = "n", main = "Ordination with Bray-Curtis dissimilarity")
orditorp(nmds.asv.t, display = "sites", labels= "n", col = as.character(species.vector), pch = pop.vector, cex=1.5) 
legend("topright", y=NULL, 
       c(expression(italic("S. dumicola")),expression(italic("S. mimosarum")),expression(italic("S. sarasinorum"))), 
       pch = 0, col = c('#d95f02', '#7570b3', '#1b9e77'), cex = 1.2)
legend("bottomright", y=NULL,  pch.df$Var1, pch = pch.df$symbol, cex = 0.8, ncol=3)
dev.off()

