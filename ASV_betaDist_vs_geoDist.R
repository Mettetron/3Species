# Make an x-y plot of beta diversity between sample ASV composition vs. 
# geographical distance between sampling sites. including linear regression line,
# correlation coefficients, and significance calculations.

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns
#         to get geographical distances between samples, you need to include 
#         latitude and longitude (in decimal degrees) as attributes in the sample data.

setwd("~/ /") 
study <- "TSPfinal"

################### IMPORT #####################################################

data.asv.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.SAMPLEtable.csv", study), row.names=1)

################### SINGLE SPECIES #############################################

# species can be changed here to run script for different species.
#species <- "mimosarum"
species <- "dumicola"
#species <- "sarasinorum"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.asv.norm <- data.asv.norm[, which(colnames(data.asv.norm) %in% rownames(data.meta.norm))]

# species colors: '#7570b3' = mim, '#d95f02' = dum, '#1b9e77' = sar
species_color <- ifelse(species=="dumicola", '#d95f02', ifelse(species=="mimosarum", '#7570b3', '#1b9e77'))

# check that things line up
sum(rownames(data.meta.norm) == colnames(data.asv.norm)) == length(rownames(data.meta.norm))

##################### SUM FOR COLONIES #########################################

# transpose
data.asv.norm.t <- as.data.frame(t(data.asv.norm))

# add column with col/pop/species you want summed
spider <- 'nest'
data.asv.norm.t$spider <- data.meta.norm[, spider]
sum.me <- data.asv.norm.t

#add together data from same spider level
require(data.table)
SUM.ME <- data.table(sum.me)
im.summed <- SUM.ME[, lapply(.SD, sum ), by = spider]

# set rownames to spider group level and order by bact abundance 
data.summed <- as.data.frame(im.summed)
rownames(data.summed) <- data.summed$spider
data.summed$spider <- NULL

# not just summed, but average of sum for each colony
hej <- table(droplevels(data.asv.norm.t$spider))
data.summed <- data.summed/hej

# reset variable
data.asv.colony <- as.data.frame(data.summed)

################### GET GEO DISTANCES ##########################################
library(sp)
# take only col, lat and lon
col.coords <- data.meta.norm[, c('nest', 'lat','lon')]
# get unique cases - so one line for each colony
col.coords <- unique(col.coords)
#set colony as rowname
rownames(col.coords) <- col.coords$nest
# remove colony column
col.coords$nest <- NULL

# make the geo distance matrix. dist in km
library(sp)
col.geo.dist <- apply(col.coords, 1, function(eachPoint) spDistsN1(as.matrix(col.coords), eachPoint, longlat=TRUE))
rownames(col.geo.dist) <- colnames(col.geo.dist)

listed.geo.dist <- melt(col.geo.dist)

##################### GET B-DIVERSITY ##########################################
library(vegan)
beta.dist <- vegdist(data.asv.colony, method="bray", binary=FALSE)
# melt
library(reshape2)
beta.dist.df <- melt(as.matrix(beta.dist), varnames = c("row", "col"))

colnames(beta.dist.df) <- c("Var1", "Var2", "value")
##################### PLOT #####################################################

# make table with both variables
beta.dist.df$Var2 == listed.geo.dist$Var2 # checking that things line up
beta.dist.df$Var1 == listed.geo.dist$Var1 # checking that things line up
listed.beta.geo.dist <- beta.dist.df
listed.beta.geo.dist$geo <- listed.geo.dist$value
colnames(listed.beta.geo.dist) <- c("colony1", "colony2", "beta", "geo")

# remove self comparison
listed.beta.geo.dist <- listed.beta.geo.dist[listed.beta.geo.dist$colony1 != listed.beta.geo.dist$colony2, ]

# remove double comparison
hej <- listed.beta.geo.dist[, 1:2]
hej.sort <- t(apply(hej, 1, sort))
listed.beta.geo.dist <- listed.beta.geo.dist[!duplicated(hej.sort), ]

# calculate correlation
hej <- cor.test(listed.beta.geo.dist$geo, listed.beta.geo.dist$beta)
cortext <- paste("cor =", round(hej$estimate, 3))
sigtext <- paste("p-val =", round(hej$p.value, 3))
label <- paste(cortext, ",", sigtext)
labels = data.frame(x = 1500, y = 0.2, label = label)

library(ggplot2)
library(grid)

# ggplot version
x <- as.vector(listed.beta.geo.dist$geo)
y <- as.vector(listed.beta.geo.dist$beta)
plot.data <- data.frame(geo = x, beta = as.numeric(y))
gg <-ggplot(plot.data, aes(x=geo, y=beta)) +
  geom_point(shape = 1, color = sprintf("%s", species_color)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid lines
  ggtitle(sprintf("Beta diversity vs Geo distance, S.%s nests", species)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Geographical distance (km)") +
  ylab("Bray-Curtis dissimilarity") +
  theme(axis.text = element_text(size = 12)) + 
  geom_smooth(method = lm, se = FALSE, color = 1, lty = "dashed", size=0.5) + # add linear regression line, don't add shaded confidence region
  annotation_custom(grob=textGrob(label, gp=gpar(fontsize=10)), xmin=1400, xmax=1700, ymin=-0.2, ymax=-0.15) + # change here for different species
  scale_y_continuous(limits = c(0, 1.0))

# sar: xmin=600, xmax=800
# mim: xmin=1700, xmax=2000
# dum: xmin=1400, xmax=1700

gg
gb <- ggplot_build(gg)
gt <- ggplot_gtable(gb)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
jpeg(sprintf("%s_BdivVSgeodiv.jpg", species), width = 18, height = 9, units="cm", res=600)
print(grid.draw(gt))
dev.off()

#  cor test. null hypothesis: slope=0
cor.test(x,as.numeric(y))

### mim
# Pearson's product-moment correlation
# 
# data:  x and as.numeric(y)
# t = 1.2861, df = 251, p-value = 0.1996
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.04284155  0.20222393
# sample estimates:
# cor 
# 0.08091394 

### dum
# Pearson's product-moment correlation
# 
# data:  x and as.numeric(y)
# t = -0.56374, df = 376, p-value = 0.5733
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.12954854  0.07201865
# sample estimates:
# cor 
# -0.02906037 

### sar
# Pearson's product-moment correlation
# 
# data:  x and as.numeric(y)
# t = 5.6446, df = 526, p-value = 2.71e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.1568533 0.3178365
# sample estimates:
# cor 
# 0.2389864 
