# make stacked barcharts to compare relative abundance of ASVs in hosts 
# at different grouping levels.

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns
#         Here I use the sample attribute "species" to group my samples, but you can change the script to use anything

setwd("~/ /")

study <- "TSPfinal"

################### IMPORT #####################################################
data.asv.norm <- read.csv(sprintf("%s_DADA2_Fraction.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Fraction.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Fraction.SAMPLEtable.csv", study), row.names=1)

################### SUM FOR SPIDER SPECIES #####################################

# transpose
data.summed <- as.data.frame(t(data.asv.norm))

sum(rownames(data.summed) == rownames(data.meta.norm)) == length(rownames(data.meta.norm))
data.summed$species <- paste0(data.meta.norm$species)

library(dplyr)
data.summed <- data.summed %>% 
  group_by(species) %>% 
  summarise_all(funs(mean), na.rm = TRUE)

data.summed <- as.data.frame(data.summed)
rownames(data.summed) <- data.summed$species
data.summed$species <- NULL


################### COLLECT OTHER ##############################################

# order by bact abundance and sum other
data.summed <- data.summed[, order(colSums(data.summed))]
data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
other <- rowSums(data.summed[, 1:(ncol(data.summed)-10)])
data.summed.other <- data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
data.summed.other <- cbind(other, data.summed.other)

# insert genus as columnnames
colnames(data.summed.other) <- c("other", 
                                 paste(as.vector(data.tax.norm[colnames(data.summed.other[, 2:11]), ]$Genus), 
                                       as.vector(rownames(data.tax.norm[colnames(data.summed.other[, 2:11]), ]))))

# # fix order of ASVs - ### NOTE!: this is only needed if you want custom order of ASVs
# my.order <- rev(c("Diplorickettsia", "Rickettsiella", "Mycoplasma", "Borrelia", "Spirochaetaceae", "Rickettsiaceae", "Chlamydiales", "Acaricomes", "Weeksellaceae", "Brevibacterium", "Spiroplasma", "Entomoplasma", "other"))
# data.summed.other.sorted <- data.summed.other[unlist(lapply(gsub("\\d+", "", my.order), function(x) grep(x, names(data.summed.other))))]
# data.summed.leftovers <- as.data.frame(data.summed.other[, !colnames(data.summed.other) %in% colnames(data.summed.other.sorted)])
# colnames(data.summed.leftovers) <- colnames(data.summed.other)[!colnames(data.summed.other) %in% colnames(data.summed.other.sorted)]
# data.summed.other.sorted$other <- NULL
# data.summed.other.final <- cbind(data.summed.other$other, data.summed.leftovers, data.summed.other.sorted)
# colnames(data.summed.other.final)[1] <- "other"
# 
# data.summed.other<- data.summed.other.final

# insert rownames as column
data.summed.other$subject <- rownames(data.summed.other)

# count number of individuals in each species
num.vec <- c()
for (species in rownames(data.summed.other)) {
  number <- nrow(data.meta.norm[data.meta.norm$species == species, ])
  num.vec <- c(num.vec, number)
}

# fix spider names
data.summed.other$subject <- paste0("italic(", "S.~", rownames(data.summed.other),")", "(", num.vec, ")")

################### PLOT AND EXPORT ############################################
# melt
library(reshape2)
chart.data <- melt(data.summed.other, id = 'subject') 
# change columnnames
colnames(chart.data) <- c('sample', 'ASV', 'count')

# # fix colors according to "color master"    ### NOTE!: this is not needed if you just want to use random colors
# color.master <- read.csv2("3SPECIES_ColorMaster.csv")
# colnames(color.master)[1] <- "genus"
# rownames(color.master) <- color.master$genus
# 
# # vector of genus names
# my.list <- strsplit(unique(as.character(chart.data$ASV)), ' ')
# chart.genus.vector <- sapply(my.list, "[[", 1)  # take first element of each vector in a list
# 
# color.vector <- c()
# hej <- table(chart.genus.vector)
# for (gen in unique(chart.genus.vector)) {
#   if (gen %in% color.master$genus) {
#     base.gen <- rep(as.character(color.master[gen, "color"]), hej[gen])
#     color.vector <- c(color.vector, base.gen)
#   } else {
#     base.gen <- rep(as.character(color.master["extra", "color"]), hej[gen])
#     color.vector <- c(color.vector, base.gen)
#   }
# }
# 
# # make extra color different for different genera
# for (farve in unique(color.vector)) {
#   this.many <- sum(color.vector == farve)
#   colfunc <- colorRampPalette(c(as.character(farve), "white"))
#   col.replace <- colfunc(this.many+2)[1:this.many]
#   color.vector[color.vector == as.character(farve)] <- rev(col.replace)
# }

library(stringr)

# fix names for some cursive, some not
df <- str_split(chart.data$ASV, " ", simplify = T)
df <- trimws(df)  # remove trailing spaces
df <- as.data.frame(df)
colnames(df) <- c("name", "number")
# format field 1 (=ASV latin name) to be cursive, unless it is "other or "subject"
exceptions <- c("other", "subject", "Chlamydiales", "Rickettsiaceae", "Weeksellaceae")  # these can be added on to if needed
exceptions2 <- c("other", "subject")

df$name <- ifelse(df$name %in%  exceptions,
                  paste0(df$name),
                  paste0("italic(", df$name, ")"))
df$longname <- paste(df$name, df$number)
df$longname <- ifelse(df$name %in%  exceptions2,
                      paste0(df$name),
                      str_replace_all(df$longname, " ", "~~"))
chart.data$newASV <- df$longname

library(RColorBrewer)
library(ggplot2)
require(scales)
library(ggthemes)

legend.labels <- c(as.character(unique(chart.data$newASV)))

# plot (note: because of cursive, this will not plot in Rstudio - export to see)
p <- ggplot(chart.data, aes(x = sample)) + 
  geom_bar(aes(weight=count, fill = ASV), position = 'fill') + 
  theme_bw(base_family="Helvetica") +
  scale_y_continuous(labels = percent) + # percentages on y-axis
  scale_fill_manual(values = (brewer.pal(11, "Spectral")), labels=parse(text=legend.labels)) +  # generic color palette
  #scale_fill_manual(values = color.vector, labels=parse(text=legend.labels)) +  # custom color palette
  theme(legend.text.align = 0) +
  theme(axis.title.x = element_blank()) +  # removing x-axis label
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = 'black')) + # changes to x-axis text
  scale_x_discrete(labels=parse(text=sort(data.summed.other$subject))) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black')) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.margin=unit(c(0.2, 0.2, 0.2, 0.2),"cm")) +  # add margin space (top, right, bottom, left). 0.2 is nice default
  ylab("Relative Abundance")  # adding y-axis label

ggsave("ASV_barchart_grouped.jpg", p,
       scale = 1, width = 17.5, height = 10, units = "cm",
       dpi = 600, limitsize = TRUE)


