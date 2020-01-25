# Make an NMDS ordination plot based on ASV data. And calculate ANOSIM.
# Here I color with specific colors by host species, and the shapes denote different host populations.
# note: Each run of the metaMDS function will give slightly different outputs, as this is an iterative method

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#         it is recommended to use subsampled data for calculation of distance metrics
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns


setwd("~/")
study <- "TSPfinal"

library(vegan)

# import subsampled data
data.asv <- read.csv(sprintf("%s_DADA2_Subsampl3000.ASVtable.csv", study), row.names=1)
data.tax <- read.csv(sprintf("%s_DADA2_Subsampl3000.TAXtable.csv", study), row.names=1)
data.meta <- read.csv(sprintf("%s_DADA2_Subsampl3000.SAMPLEtable.csv", study), row.names=1)

# make population symbols
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

# margin default: par(mar=c(5.1,4.1,4.1,2.1), par(mar = c(bottom, left, top, right))
par(mar=c(5.1,4.1,4.1,10))
pdf(file="NMDS.pdf", width = 9, height = 9) # this prints
plot(nmds.asv.t, display = "sites", type = "n", main = "Ordination with Bray-Curtis dissimilarity")
orditorp(nmds.asv.t, display = "sites", labels= "n", col = as.character(species.vector), pch = pop.vector, cex=1.5)
legend("topright", y=NULL,
       c(expression(italic("S. dumicola")),expression(italic("S. mimosarum")),expression(italic("S. sarasinorum"))),
       pch = 0, col = c('#d95f02', '#7570b3', '#1b9e77'), cex = 1.2)
legend("bottomright", y=NULL,  pch.df$Var1, pch = pch.df$symbol, cex = 0.8, ncol=3)
dev.off()

################ ANOSIM ######################################################## 
# (needs package: vegan)
# anosim only gives one "between" statistic for each df, so need to split up 
# data for each comparison

# get species combinations
comb <- combn(unique(data.meta$species), 2, FUN = NULL, simplify = TRUE)  # get all combinations of your variable
# make df for results
ano.species.df <- as.data.frame(matrix(NA, nrow=1, ncol=4))
colnames(ano.species.df) <- c("ind1", "ind2", "ano_R", "ano_p")

# loop over combination pairs, and run anosim on data containing only your pair
for (i in 1:ncol(comb)) {
  ind1 <- comb[, i][1]
  ind2 <- comb[, i][2]
  data.m1 <- data.meta[data.meta$species == ind1, ] 
  data.m2 <- data.meta[data.meta$species == ind2, ]
  data.m <- rbind(data.m1, data.m2)
  data.a <- as.data.frame(t(data.asv[, rownames(data.m)])) # extract and transform asv data
  ano <- with(data.m, anosim(data.a, species))
  r.vector <- c(as.character(ind1), as.character(ind2), round(ano$statistic, 3), round(ano$signif, 3))
  ano.species.df <- rbind(ano.species.df, r.vector)
}
ano.species.df <- ano.species.df[2:nrow(ano.species.df), ]


# anosim on all population pairs
for (s in c("dumicola", "mimosarum", "sarasinorum")) {
  #s <- "dumicola"
  data.meta.s <- data.meta[data.meta$species == s, ]
  comb <- combn(unique(data.meta.s$population), 2, FUN = NULL, simplify = TRUE)  # get all combinations of your variable
  ano.pop.df <- as.data.frame(matrix(NA, nrow=1, ncol=4))
  colnames(ano.pop.df) <- c("ind1", "ind2", "ano_R", "ano_p")
  for (i in 1:ncol(comb)) {
    #i <- 1
    ind1 <- as.character(comb[, i][1])
    ind2 <- as.character(comb[, i][2])
    data.m1 <- data.meta.s[data.meta.s$population == ind1, ] 
    data.m2 <- data.meta.s[data.meta.s$population == ind2, ]
    data.m <- rbind(data.m1, data.m2)
    data.a <- as.data.frame(t(data.asv[, rownames(data.m)])) # extract and transform asv data
    print(identical(rownames(data.a), rownames(data.m)))
    ano <- with(data.m, anosim(data.a, population))
    r.vector <- c(ind1, ind2, round(ano$statistic, 3), round(ano$signif, 3))
    ano.pop.df <- rbind(ano.pop.df, r.vector)
  }
  ano.pop.df <- ano.pop.df[2:nrow(ano.pop.df), ] # get rid of that first NA row
  assign(paste0("ano.", s, ".pop.df"), ano.pop.df)
}


### make into table
s.er <- c("dumicola", "mimosarum", "sarasinorum")
for (s in s.er) {
  s <- "sarasinorum"
  ano.df <- get(paste0("ano.", s, ".pop.df"))
  ano.df <- ano.df[ano.df$ind1 != "SriLanka", ]  # remove SriLanka
  ano.df <- ano.df[ano.df$ind2 != "SriLanka", ]  # remove SriLanka
  ano.df$ano_R <- as.numeric(ano.df$ano_R)  # variable to numeric
  ano.df$ano_p <- as.numeric(ano.df$ano_p)
  data.meta.s <- data.meta[data.meta$species == s, ]
  data.meta.s <- data.meta.s[data.meta.s$population != "SriLanka", ]  # remove SriLanka
  ano.groups <- unique(data.meta.s[, "population"])
  group.nr <- length(ano.groups)
  
  # make number vector for splitting up into table:
  n <- group.nr - 1
  x <- 1:n
  id <- unlist(sapply(x, function(x) x:n))
  ano.df$id <- id
  
  library(reshape2)
  ano.table <- dcast(data=ano.df, formula = id~ind1, fun.aggregate = sum, value.var = "ano_R")
  ano.table.p <- dcast(data=ano.df, formula = id~ind1, fun.aggregate = sum, value.var = "ano_p")
  
  # table cleanup
  ano.rows <- ano.groups[2:group.nr]
  rownames(ano.table) <- ano.rows
  rownames(ano.table.p) <- ano.rows
  ano.table <- ano.table[, as.character(ano.groups[1:(group.nr-1)])]
  ano.table.p <- ano.table.p[, as.character(ano.groups[1:(group.nr-1)])]
  ano.table <- ano.table[rev(rownames(ano.table)), ]
  ano.table.p <- ano.table.p[rev(rownames(ano.table.p)), ]
  colnames(ano.table.p) <- paste0(colnames(ano.table.p), "_p")
  
  
  # table export
  ano.table.full <- cbind(ano.table, ano.table.p)
  write.csv(ano.table.full, sprintf("ANOSIMtable_pop_%s_NoSriLanka.csv", s))
}

