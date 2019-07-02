# make stacked barchart with each sample appearing as a column

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns
#         Here I use the sample attributes "nest" and "population" for shoving groupings in the plot,
#         but that can easily be changed or deleted.

setwd("~/ /") # CHANGE HERE TO MATCH OTHER DATA / SETUP
study <- "TSPfinal"# CHANGE HERE TO MATCH OTHER DATA / SETUP

################### IMPORT #####################################################
data.asv.norm <- read.csv(sprintf("%s_DADA2_Fraction.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Fraction.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Fraction.SAMPLEtable.csv", study), row.names=1)

################### JUST DUMICOLA ##############################################

# only dumicola
species <- "dumicola"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.asv.norm <- data.asv.norm[, which(colnames(data.asv.norm) %in% rownames(data.meta.norm))]

# change sample names to tube IDs
sum(rownames(data.meta.norm) == colnames(data.asv.norm)) == length(rownames(data.meta.norm))
rownames(data.meta.norm) <- data.meta.norm$tube_ID
colnames(data.asv.norm) <- data.meta.norm$tube_ID

################### COLLECT OTHER ##############################################
data.summed <- as.data.frame(t(data.asv.norm))

# order by bact abundance and sum other
data.summed <- data.summed[, order(colSums(data.summed))]
other <- rowSums(data.summed[, 1:(ncol(data.summed)-10)])
data.summed.other <- data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
data.summed.other <- cbind(other, data.summed.other)

# insert genus as columnnames
colnames(data.summed.other) <- c("other",
                                 paste(as.vector(data.tax.norm[colnames(data.summed.other[, 2:11]), ]$Genus),
                                       as.vector(rownames(data.tax.norm[colnames(data.summed.other[, 2:11]), ]))))


# insert nest names in rownames
sum(rownames(data.summed.other) == rownames(data.meta.norm)) == length(rownames(data.summed.other))
rownames(data.summed.other) <- paste(data.meta.norm$nest, rownames(data.summed.other))

# insert rownames as column
data.summed.other$subject <- rownames(data.summed.other)
colnames(data.summed.other)

# manually fix order so colors can fit with the grouped barchert and still look good
column.order <- c("other", "Weeksellaceae ASV_47", "Weeksellaceae ASV_19", "Acaricomes ASV_15", "Rickettsia ASV_5", "Rickettsiaceae ASV_6", "Borrelia ASV_3", "Mycoplasma ASV_26", "Mycoplasma ASV_4", "Diplorickettsia ASV_9", "Diplorickettsia ASV_2", "subject")
data.summed.other <- data.summed.other[, column.order]

# melt
library(reshape2)
chart.data <- melt(data.summed.other, id = 'subject') 
# change columnnames
colnames(chart.data) <- c('sample', 'ASV', 'count')

# fix colors according to "color master"    ### NOTE!: this is not needed if you just want to use random colors
color.master <- read.csv2("3SPECIES_ColorMaster.csv")
colnames(color.master)[1] <- "genus"
rownames(color.master) <- color.master$genus

# vector of genus names
my.list <- strsplit(unique(as.character(chart.data$ASV)), ' ')
chart.genus.vector <- sapply(my.list, "[[", 1)  # take first element of each vector in a list

color.vector <- c()
hej <- table(chart.genus.vector)
for (gen in unique(chart.genus.vector)) {
  if (gen %in% color.master$genus) {
    base.gen <- rep(as.character(color.master[gen, "color"]), hej[gen])
    color.vector <- c(color.vector, base.gen)
  } else {
    base.gen <- rep(as.character(color.master["extra", "color"]), hej[gen])
    color.vector <- c(color.vector, base.gen)
  }
}

# make extra color different for different genera
for (farve in unique(color.vector)) {
  this.many <- sum(color.vector == farve)
  colfunc <- colorRampPalette(c(as.character(farve), "white"))
  col.replace <- colfunc(this.many+2)[1:this.many]
  color.vector[color.vector == as.character(farve)] <- rev(col.replace)
}

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


################### PLOT AND EXPORT ############################################

# insert nest name for plot split
nests <- strsplit(as.character(chart.data$sample),' ') 
nests <- as.data.frame(do.call(rbind, nests))
chart.data$nest <- nests$V1
# remove nest names from sample name
chart.data$subject <- nests$V2
# change columnnames
colnames(chart.data) <- c('sample', 'ASV', 'count', 'newASV', 'nest', 'subject')

library(RColorBrewer)
library(ggplot2)
require(scales)
library(ggthemes)

pops <- strsplit(as.character(chart.data$nest),'_')
pops <- as.data.frame(do.call(rbind, pops))
chart.data$pop <- pops$V1

#### Plot and then add nested facets
library('gtable')
library('grid')
library('magrittr')

# make nest names only numbers
library(stringr)  # string things
nestdf <- as.data.frame(str_split(chart.data$nest, "_", simplify = T))
chart.data$nest <- nestdf$V2

legend.labels <- c(as.character(unique(chart.data$newASV)))

myplot <- ggplot(chart.data, aes(x = subject)) + geom_bar(aes(weight=count, fill = ASV), position = 'fill') + 
  scale_y_continuous(labels = percent) + # percentages on y-axis
  #scale_fill_manual(values = (brewer.pal(11, "Spectral"))) +  # generic color palette
  scale_fill_manual(values = color.vector, labels=parse(text=legend.labels)) + # custom color palette and names
  theme_bw(base_family="Helvetica") + # white background
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust=1, vjust=0.5, color = 'black')) + # changes to x-axis text
  facet_grid(.~pop+nest, scales="free_x", space="free_x") + # double facet, to give both nest and population (modified below)
  theme(strip.text = element_text(size=6.5)) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.text.align = 0) +
  theme(panel.spacing= unit(0.1, "cm")) +
  theme(plot.margin=unit(c(0.2, 0.2, 0.2, 0.2),"cm")) +  # add margin space (top, right, bottom, left). 0.2 is nice default
  xlab("Individual spiders, sample names") +
  ylab("Relative Abundance")  # adding y-axis label

### NOTE! this plot does not work on Rstudio built in graphics device. but works fine when exported.
# ggsave("test.png", myplot, 
#        scale = 1, width = 20, height = 10, units = "cm",
#        dpi = 600, limitsize = TRUE)

##### special for nested facets, including custom fitting = changes need to be made for different data
# nested facet code inspired by: https://stackoverflow.com/questions/40316169/nested-facets-in-ggplot2-spanning-groups

z <- ggplotGrob(myplot)

# Find the location of the strips in the main plot
locations <- grep("strip-t", z$layout$name)
# Filter out the strips (trim = FALSE is important here for positions relative to the main plot)
strip <- gtable_filter(z, "strip-t", trim = FALSE)

# Gathering our positions for the main plot
top <- strip$layout$t[1]
l   <- strip$layout$l[c(1, 8, 14, 20, 24)] # ADDO goes from 1-7, KRU goes from 8-13, PAA from 14-19, PON from 20-23, WEE from 24-28
r   <- strip$layout$r[c(7, 13, 19, 23, 28)]

mat   <- matrix(vector("list", length = 6), nrow = 2)
mat[] <- list(zeroGrob())

# The separator for the facets has zero width
res <- gtable_matrix("toprow", mat, unit(c(1, 0, 1), "null"), unit(c(1, 1), "null"))

# Adding the first group LOCATION is 1, because that is first facet with ADDO
zz <- res %>%
  gtable_add_grob(z$grobs[[locations[1]]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(z, ., t = top,  l = l[1],  b = top,  r = r[1], name = c("add-strip"))

# Adding the second group (note the indices) LOCATION is 8, because that is first facet with KRU
xx <- gtable_add_grob(res, z$grobs[[locations[8]]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(zz, ., t = top,  l = l[2],  b = top,  r = r[2], name = c("add-strip"))

# Adding the third group (note the indices) LOCATION is 14, because that is first facet with PAA
yy <- gtable_add_grob(res, z$grobs[[locations[14]]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(xx, ., t = top,  l = l[3],  b = top,  r = r[3], name = c("add-strip"))

# Adding the fouth group (note the indices) LOCATION is 20, because that is first facet with PON
mm <- gtable_add_grob(res, z$grobs[[locations[20]]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(yy, ., t = top,  l = l[4],  b = top,  r = r[4], name = c("add-strip"))

# Adding the fifth group (note the indices) LOCATION is 24, because that is first facet with WEE
pp <- gtable_add_grob(res, z$grobs[[locations[24]]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(mm, ., t = top,  l = l[5],  b = top,  r = r[5], name = c("add-strip"))

grid.newpage()
jpeg(sprintf("Barchart_individuals_%s.jpg", species), width = 30, height = 12, units="cm", res=600)
print(grid.draw(pp))
dev.off()


