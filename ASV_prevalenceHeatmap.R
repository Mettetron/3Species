# Make a prevalence heatmap of your most prevalent ASVs
# here it is represented as % occurence (above a certain threshold) in host populations,
# and it is plotted in panels by host species. 
# These can be modified to fit attributes of any sample

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns


setwd("~/")

study <- "TSPfinal"

################### IMPORT #####################################################
data.asv.norm <- read.csv(sprintf("%s_DADA2_Fraction.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Fraction.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Fraction.SAMPLEtable.csv", study), row.names=1)

# merge tax and asv
data.tax.asv.norm <- cbind(data.tax.norm, data.asv.norm)

#################### FIX NAMES #################################################

# change certain population names to avoid overlap between species
species.letter <- substr(data.meta.norm$species, 1, 3)
data.meta.norm$population <- paste0("S.", species.letter, " ", data.meta.norm$population)
                                   
# count number of individualss in each population and add this to population name also
pop.n <- c()
for (pop in data.meta.norm$population) {
  n <- nrow(data.meta.norm[data.meta.norm$population == pop, ])
  pop.n <- c(pop.n, n)
}

data.meta.norm$population <- paste(data.meta.norm$population, "(", pop.n, ")", sep = "")

################# FIND CORE ##################################################
# this is a way to select ASVs to plot by their prevalence, see below for plotting of specific ASVs


# # make presence/absence data
# library(vegan)
# data.asv.pa <- data.asv.norm
# threshold <- 0.0001  # for fraction data, this is 0.01%
# data.asv.pa[data.asv.pa < threshold] <- 0 
# data.asv.pa <- decostand(data.asv.pa, method="pa")
# 
# # now sum up all the rows 
# asv.sums <- rowSums(data.asv.pa)
# 
# # just get all asvs ordered by most core-like
# data.tax.asv.corenum <- data.tax.asv.norm
# data.tax.asv.corenum$fiveplus <- asv.sums
# data.tax.asv.corenum.sorted <- data.tax.asv.corenum[order(-data.tax.asv.corenum$fiveplus), ]
# 
# # set minimum number to be called core:
# # must be present in 20%
# core.num <- ncol(data.asv.norm)/5
# 
# # only keep ASVs above core num
# data.tax.asv.corenum.sorted.core <- data.tax.asv.corenum.sorted[data.tax.asv.corenum.sorted$fiveplus > core.num, ]
# 
# # get shorter name
# data.core <- data.tax.asv.corenum.sorted.core



## core based on other data
data.core <- data.tax.asv.norm[c("ASV_2", "ASV_4", "ASV_3", 
                                 "ASV_1", "ASV_15", "ASV_6", 
                                 "ASV_8", "ASV_10", "ASV_47", 
                                 "ASV_11", "ASV_14", "ASV_7", "ASV_25", "ASV_44",
                                 "ASV_43","ASV_21", "ASV_65",
                                 "ASV_68", "ASV_45", "ASV_16"), ]

################### GET PERCENTAGES ############################################

# make a list of unique population names
pop.list <- unique(data.meta.norm$population)

# separate asv and tax and put asv column in data.core.tax
data.core.asv <- data.core[, 8:ncol(data.core)] # getting rid of the fiveplus column
data.core.tax <- data.core[, 1:7]
data.core.tax$asv <- rownames(data.core.tax)

# test that asv and meta data matches
colnames(data.core.asv) == rownames(data.meta.norm)

# make p/a again
threshold <- 0.0001  # for fraction data, this is 0.01%
data.core.asv.t <- as.data.frame(t(data.core.asv))
#data.core.asv.t[data.core.asv.t < 5] <- 0 # threshold for presence = 5
data.core.asv.t[data.core.asv.t < threshold] <- 0 # threshold for presence = 0.01%
data.core.asv.t.pa <- decostand(data.core.asv.t, method="pa")

# initialize result data frame
percent.per.pop <- data.core.asv[, 1:2]

# loop over pops and take % of indv in pop or col containing the asv
for (pop in pop.list){
  pop.matrix <- data.core.asv.t.pa[data.meta.norm$population == pop, ] # change here if you're doing colony or population!
  pop.asvs <- colSums(pop.matrix)/nrow(pop.matrix)*100
  percent.per.pop[, pop] <- pop.asvs
}

# remove those first two columns
percent.per.pop[, 1] <- NULL
percent.per.pop[, 1] <- NULL

# insert genus names
rownames(percent.per.pop) <- paste(data.core.tax$Genus, data.core.tax$asv)

### Order rows manually
rownames(percent.per.pop)
new_order <- c("Diplorickettsia ASV_2", "Mycoplasma ASV_4", "Borrelia ASV_3", 
                            "Chlamydiales ASV_1", "Acaricomes ASV_15", "Rickettsiaceae ASV_6", 
                            "Brevibacterium ASV_8", "Brevibacterium ASV_10", "Weeksellaceae ASV_47", 
                            "Entomoplasma ASV_11", "Entomoplasma ASV_14", "Entomoplasma ASV_7", "Entomoplasma ASV_25", "Entomoplasma ASV_44",
                            "Staphylococcus ASV_43","Pseudomonas ASV_21", "Enterobacter ASV_65",
                            "Enterobacter ASV_68", "Leucobacter ASV_45", "Delftia ASV_16")

percent.per.pop.ordered <- percent.per.pop[new_order,]

# make matrix
percentpop.matrix <- as.matrix(percent.per.pop.ordered)
#percentpop.matrix <- as.matrix(percent.per.pop)

percentpop.matrix.t <- t(percentpop.matrix)

################### PLOT AND SAVE ##############################################

library(reshape2)
library(ggplot2)
library(tidyverse)
dat <- percentpop.matrix.t

# fix names for some cursive, some not
df <- str_split(colnames(dat), " ", simplify = T)
df <- trimws(df)  # remove trailing spaces
df <- as.data.frame(df)
colnames(df) <- c("name", "number")

# fix italics on latin names, but only if classified to genus level
exceptions <- c("Chlamydiales", "Weeksellaceae", "Rickettsiaceae")
df$name <- ifelse(df$name %in%  exceptions,
                  paste0(df$name),
                  paste0("italic(", df$name, ")"))
df$longname <- paste(df$name, df$number)
df$longname <- str_replace_all(df$longname, " ", "~~")
colnames(dat) <- df$longname


# melt data
dat2 <- melt(dat)

# preparing for facets
dat2$species <-  dat2$Var1
species.names <- strsplit(as.character(dat2$species),' ') 
species.names <- as.data.frame(do.call(rbind, species.names))
dat2$species <- species.names$V1
dat2$species <- gsub('S.dum', 'S. dumicola', dat2$species)
dat2$species <- gsub('S.mim', 'S. mimosarum', dat2$species)
dat2$species <- gsub('S.sar', 'S. sarasinorum', dat2$species)

dat2$Var1 <- species.names$V2

# to make the legend look right
colnames(dat2)[3] <- "Prevalence"
dat2$Prevalence <- dat2$Prevalence/100


library(ggthemes)
library(scales)

## PLOT ###
# this doesn't alwayswork with RStudio plot device, but is fine if exported
p <- ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
  theme_bw() + 
  geom_tile(aes(fill = Prevalence)) + 
  scale_fill_gradient(low = "white", high = "darkgreen", labels = percent, name="Prevalence\nin population") +
  theme(text=element_text(family="Helvetica"), axis.title.y = element_blank(), axis.title.x=element_blank()) +
  #xlab("Spider populations") +
  theme(axis.title = element_blank()) +
  scale_y_discrete(breaks=unique(dat2$Var2), labels=parse(text=as.character(unique(dat2$Var2)))) +  
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10)) +
  facet_grid(.~species, scales="free_x", space="free_x") +
  theme(strip.text=element_text(face = "italic", size=10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = 'black')) +
  theme(axis.text.y = element_text(size = 8, color = 'black'))

p

# export
ggsave("Prevalence_heatmap_populations.png", p, 
       scale = 1, width = 20, height = 15, units = "cm",
       dpi = 600, limitsize = TRUE)


