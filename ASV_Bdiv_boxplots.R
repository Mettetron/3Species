# make boxplot(s) comparing beta diversities between host populations and nests 

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#         it is recommended to use subsampled data for calculation of distance metrics
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns

setwd("~/") 
study <- "TSPfinal"

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
################### IMPORT #####################################################

data.asv <- read.csv(sprintf("%s_DADA2_Subsampl3000.ASVtable.csv", study), row.names=1)
data.tax <- read.csv(sprintf("%s_DADA2_Subsampl3000.TAXtable.csv", study), row.names=1)
data.meta <- read.csv(sprintf("%s_DADA2_Subsampl3000.SAMPLEtable.csv", study), row.names=1)

# threshold
#data.asv[data.asv < 3] <- 0

################### DUMICOLA #####################################################
data.asv.norm <- data.asv
data.tax.norm <- data.tax
data.meta.norm <- data.meta

species <- "dumicola"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.asv.norm <- data.asv.norm[, which(colnames(data.asv.norm) %in% rownames(data.meta.norm))]

species_color <- "black"

# transpose
data.asv.norm.t <- as.data.frame(t(data.asv.norm))

# get Bray-Curtis
library(vegan)
beta.dist <- vegdist(data.asv.norm.t, method="bray", binary=FALSE)

##### make boxplots

# melt
library(reshape2)
beta.dist.df <- melt(as.matrix(beta.dist), varnames = c("row", "col"))

# remove self comparison
beta.dist.df <- beta.dist.df[beta.dist.df$row != beta.dist.df$col, ]

# remove double comparison
hej <- beta.dist.df[, 1:2]
hej.sort <- t(apply(hej, 1, sort))
beta.dist.df.unique <- beta.dist.df[!duplicated(hej.sort), ]

################# SORT BY COMPARISON TYPE ######################################

# insert population
beta.dist.df.unique$pop1 <- data.meta.norm$population[match(beta.dist.df.unique$row, rownames(data.meta.norm))]
beta.dist.df.unique$pop2 <- data.meta.norm$population[match(beta.dist.df.unique$col, rownames(data.meta.norm))]

# insert nest
beta.dist.df.unique$nest1 <- data.meta.norm$nest[match(beta.dist.df.unique$row, rownames(data.meta.norm))]
beta.dist.df.unique$nest2 <- data.meta.norm$nest[match(beta.dist.df.unique$col, rownames(data.meta.norm))]

# renaming
bddu <- beta.dist.df.unique

# different populations
bd.diff.pop <- bddu[bddu$pop1 != bddu$pop2, ]
# same population
bd.same.pop <- bddu[bddu$pop1 == bddu$pop2, ]
# same nest
bd.same.nest <- bddu[bddu$nest1 == bddu$nest2, ]

library(ggplot2)
# prep data
same_nest <- cbind(as.numeric(as.character(bd.same.nest$value)), rep("within nests", nrow(bd.same.nest)))
same_pop <- cbind(as.numeric(as.character(bd.same.pop$value)), rep("within populations", nrow(bd.same.pop)))
diff_pop <- cbind(as.numeric(as.character(bd.diff.pop$value)), rep("between populations", nrow(bd.diff.pop)))
anova.data <- as.data.frame(rbind(same_nest, same_pop, diff_pop))   
colnames(anova.data) <- c("Bdiv", "group")
anova.data$group <- as.factor(anova.data$group)
anova.data$Bdiv <- as.numeric(as.character(anova.data$Bdiv))
anova.data$group <- factor(anova.data$group,
                           levels = c("within nests","within populations", "between populations"),ordered = TRUE)


## make significance bars
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(1.03, 1.05, 1.05, 1.03))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(1.03, 1.05, 1.05, 1.03))
df3 <- data.frame(a = c(1.1, 1.1, 2.9, 2.9), b = c(1.1, 1.12, 1.12, 1.1))

# compare groups and find p value
hej <- pairwise.t.test(anova.data$Bdiv, anova.data$group,
                            p.adjust.method = "BH")

hej2 <- hej$p.value
#View(hej2)

# insert p values
nest.spop <- paste0("p = ", ifelse(hej2[1,1] < 0.01, as.character(formatC(hej2[1,1], format = "e", digits = 2)), as.character(round(hej2[1,1], 3))))
spop.dpop <- paste0("p = ", ifelse(hej2[2,2] < 0.01, as.character(formatC(hej2[2,2], format = "e", digits = 2)), as.character(round(hej2[2,2], 3))))
nest.dpop <- paste0("p = ", ifelse(hej2[2,1] < 0.01, as.character(formatC(hej2[2,1], format = "e", digits = 2)), as.character(round(hej2[2,1], 3))))
# # OR significance indication
# nest.spop <- ifelse(hej2[1,1] <= 0.01, "*", "not significant")
# spop.dpop <- ifelse(hej2[2,2] <= 0.01, "*", "not significant")
# nest.dpop <- ifelse(hej2[2,1] <= 0.01, "*", "not significant")

plot.title <- ifelse(species == "dumicola", expression(italic("S. dumicola")), ifelse(species == "mimosarum", expression(italic("S. mimosarum")), expression(italic("S. sarasinorum"))))

p.dum <- ggplot(anova.data, aes(x=group, y=Bdiv)) +
  geom_boxplot(color=species_color, outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle(plot.title) +
  theme(plot.title = element_text(size=16)) +
  scale_y_continuous(limits=c(0,1.15), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.2)) +
  ylab(label="Bray-Curtis Dissimilarity") +  # insert yaxis label
  theme(axis.title.y = element_text(size = 14, color="black")) +  # modify size of yaxis label
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 10, color="black")) + # changes to x-axis text
  scale_x_discrete(labels=addline_format(c("within nests", "within populations", "between populations")))  +
  theme(axis.text.y = element_text(color="black")) +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=3.5, stroke = 1.5, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.08, label = nest.spop, size = 2.5) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 1.08, label = spop.dpop, size = 2.5) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.15, label = nest.dpop, size = 2.5)



################### MIMOSARUM #####################################################

data.asv.norm <- data.asv
data.tax.norm <- data.tax
data.meta.norm <- data.meta

species <- "mimosarum"

data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.asv.norm <- data.asv.norm[, which(colnames(data.asv.norm) %in% rownames(data.meta.norm))]

species_color <- "black"

# transpose
data.asv.norm.t <- as.data.frame(t(data.asv.norm))

# get Bray-Curtis
library(vegan)
beta.dist <- vegdist(data.asv.norm.t, method="bray", binary=FALSE)


##### make boxplots

# melt
library(reshape2)
beta.dist.df <- melt(as.matrix(beta.dist), varnames = c("row", "col"))

# remove self comparison
beta.dist.df <- beta.dist.df[beta.dist.df$row != beta.dist.df$col, ]

# remove double comparison
hej <- beta.dist.df[, 1:2]
hej.sort <- t(apply(hej, 1, sort))
beta.dist.df.unique <- beta.dist.df[!duplicated(hej.sort), ]

################# SORT BY COMPARISON TYPE ######################################

# insert population
beta.dist.df.unique$pop1 <- data.meta.norm$population[match(beta.dist.df.unique$row, rownames(data.meta.norm))]
beta.dist.df.unique$pop2 <- data.meta.norm$population[match(beta.dist.df.unique$col, rownames(data.meta.norm))]

# insert nest
beta.dist.df.unique$nest1 <- data.meta.norm$nest[match(beta.dist.df.unique$row, rownames(data.meta.norm))]
beta.dist.df.unique$nest2 <- data.meta.norm$nest[match(beta.dist.df.unique$col, rownames(data.meta.norm))]

# renaming
bddu <- beta.dist.df.unique

# different populations
bd.diff.pop <- bddu[bddu$pop1 != bddu$pop2, ]
# same population
bd.same.pop <- bddu[bddu$pop1 == bddu$pop2, ]
# same nest
bd.same.nest <- bddu[bddu$nest1 == bddu$nest2, ]


library(ggplot2)
# prep data
same_nest <- cbind(as.numeric(as.character(bd.same.nest$value)), rep("within nests", nrow(bd.same.nest)))
same_pop <- cbind(as.numeric(as.character(bd.same.pop$value)), rep("within populations", nrow(bd.same.pop)))
diff_pop <- cbind(as.numeric(as.character(bd.diff.pop$value)), rep("between populations", nrow(bd.diff.pop)))
anova.data <- as.data.frame(rbind(same_nest, same_pop, diff_pop))   
colnames(anova.data) <- c("Bdiv", "group")
anova.data$group <- as.factor(anova.data$group)
anova.data$Bdiv <- as.numeric(as.character(anova.data$Bdiv))
anova.data$group <- factor(anova.data$group,
                           levels = c("within nests","within populations", "between populations"),ordered = TRUE)


## make significance bars
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(1.03, 1.05, 1.05, 1.03))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(1.03, 1.05, 1.05, 1.03))
df3 <- data.frame(a = c(1.1, 1.1, 2.9, 2.9), b = c(1.1, 1.12, 1.12, 1.1))

# compare groups and find p value
hej <- pairwise.t.test(anova.data$Bdiv, anova.data$group,
                            p.adjust.method = "BH")

hej2 <- hej$p.value
#View(hej2)

# insert p values
nest.spop <- paste0("p = ", ifelse(hej2[1,1] < 0.01, as.character(formatC(hej2[1,1], format = "e", digits = 2)), as.character(round(hej2[1,1], 3))))
spop.dpop <- paste0("p = ", ifelse(hej2[2,2] < 0.01, as.character(formatC(hej2[2,2], format = "e", digits = 2)), as.character(round(hej2[2,2], 3))))
nest.dpop <- paste0("p = ", ifelse(hej2[2,1] < 0.01, as.character(formatC(hej2[2,1], format = "e", digits = 2)), as.character(round(hej2[2,1], 3))))
# # OR significance indication
# nest.spop <- ifelse(hej2[1,1] <= 0.01, "*", "not significant")
# spop.dpop <- ifelse(hej2[2,2] <= 0.01, "*", "not significant")
# nest.dpop <- ifelse(hej2[2,1] <= 0.01, "*", "not significant")

plot.title <- ifelse(species == "dumicola", expression(italic("S. dumicola")), ifelse(species == "mimosarum", expression(italic("S. mimosarum")), expression(italic("S. sarasinorum"))))

p.mim <- ggplot(anova.data, aes(x=group, y=Bdiv)) +
  geom_boxplot(color=species_color, outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle(plot.title) +
  theme(plot.title = element_text(size=16)) +
  scale_y_continuous(limits=c(0,1.15), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.2)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +  
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 10, color="black")) + # changes to x-axis text
  scale_x_discrete(labels=addline_format(c("within nests", "within populations", "between populations")))  +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=3.5, stroke = 1.5, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.08, label = nest.spop, size = 2.5) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 1.08, label = spop.dpop, size = 2.5) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.15, label = nest.dpop, size = 2.5)



################### SARASINORUM #####################################################

data.asv.norm <- data.asv
data.tax.norm <- data.tax
data.meta.norm <- data.meta


species <- "sarasinorum"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.asv.norm <- data.asv.norm[, which(colnames(data.asv.norm) %in% rownames(data.meta.norm))]

species_color <- "black"

# transpose
data.asv.norm.t <- as.data.frame(t(data.asv.norm))

# get Bray-Curtis
library(vegan)
beta.dist <- vegdist(data.asv.norm.t, method="bray", binary=FALSE)


##### make boxplots

# melt
library(reshape2)
beta.dist.df <- melt(as.matrix(beta.dist), varnames = c("row", "col"))

# remove self comparison
beta.dist.df <- beta.dist.df[beta.dist.df$row != beta.dist.df$col, ]

# remove double comparison
hej <- beta.dist.df[, 1:2]
hej.sort <- t(apply(hej, 1, sort))
beta.dist.df.unique <- beta.dist.df[!duplicated(hej.sort), ]

################# SORT BY COMPARISON TYPE ######################################

# insert population
beta.dist.df.unique$pop1 <- data.meta.norm$population[match(beta.dist.df.unique$row, rownames(data.meta.norm))]
beta.dist.df.unique$pop2 <- data.meta.norm$population[match(beta.dist.df.unique$col, rownames(data.meta.norm))]

# insert nest
beta.dist.df.unique$nest1 <- data.meta.norm$nest[match(beta.dist.df.unique$row, rownames(data.meta.norm))]
beta.dist.df.unique$nest2 <- data.meta.norm$nest[match(beta.dist.df.unique$col, rownames(data.meta.norm))]

# renaming
bddu <- beta.dist.df.unique

# different populations
bd.diff.pop <- bddu[bddu$pop1 != bddu$pop2, ]
# same population
bd.same.pop <- bddu[bddu$pop1 == bddu$pop2, ]
# same nest
bd.same.nest <- bddu[bddu$nest1 == bddu$nest2, ]


library(ggplot2)
# prep data
same_nest <- cbind(as.numeric(as.character(bd.same.nest$value)), rep("within nests", nrow(bd.same.nest)))
same_pop <- cbind(as.numeric(as.character(bd.same.pop$value)), rep("within populations", nrow(bd.same.pop)))
diff_pop <- cbind(as.numeric(as.character(bd.diff.pop$value)), rep("between populations", nrow(bd.diff.pop)))
anova.data <- as.data.frame(rbind(same_nest, same_pop, diff_pop))   
colnames(anova.data) <- c("Bdiv", "group")
anova.data$group <- as.factor(anova.data$group)
anova.data$Bdiv <- as.numeric(as.character(anova.data$Bdiv))
anova.data$group <- factor(anova.data$group,
                           levels = c("within nests","within populations", "between populations"),ordered = TRUE)


## make significance bars
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(1.03, 1.05, 1.05, 1.03))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(1.03, 1.05, 1.05, 1.03))
df3 <- data.frame(a = c(1.1, 1.1, 2.9, 2.9), b = c(1.1, 1.12, 1.12, 1.1))

# compare groups and find p value
hej <- pairwise.t.test(anova.data$Bdiv, anova.data$group,
                            p.adjust.method = "BH")

hej2 <- hej$p.value
#View(hej2)

# insert p values
nest.spop <- paste0("p = ", ifelse(hej2[1,1] < 0.01, as.character(formatC(hej2[1,1], format = "e", digits = 2)), as.character(round(hej2[1,1], 3))))
spop.dpop <- paste0("p = ", ifelse(hej2[2,2] < 0.01, as.character(formatC(hej2[2,2], format = "e", digits = 2)), as.character(round(hej2[2,2], 3))))
nest.dpop <- paste0("p = ", ifelse(hej2[2,1] < 0.01, as.character(formatC(hej2[2,1], format = "e", digits = 2)), as.character(round(hej2[2,1], 3))))
# # OR significance indication
# nest.spop <- ifelse(hej2[1,1] <= 0.01, "*", "not significant")
# spop.dpop <- ifelse(hej2[2,2] <= 0.01, "*", "not significant")
# nest.dpop <- ifelse(hej2[2,1] <= 0.01, "*", "not significant")

plot.title <- ifelse(species == "dumicola", expression(italic("S. dumicola")), ifelse(species == "mimosarum", expression(italic("S. mimosarum")), expression(italic("S. sarasinorum"))))

p.sar <- ggplot(anova.data, aes(x=group, y=Bdiv)) +
  geom_boxplot(color=species_color, outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle(plot.title) +
  theme(plot.title = element_text(size=16)) +
  scale_y_continuous(limits=c(0,1.15), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.2)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +  
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 10, color="black")) + # changes to x-axis text
  scale_x_discrete(labels=addline_format(c("within nests", "within populations", "between populations")))  +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=3.5, stroke = 1.5, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.08, label = nest.spop, size = 2.5) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 1.08, label = spop.dpop, size = 2.5) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.15, label = nest.dpop, size = 2.5)


library(gridExtra)
p.all <- gtable_cbind(ggplotGrob(p.dum + theme(plot.margin = unit(c(5.1, 2, 4.1, 5.1), "pt"))),
                      ggplotGrob(p.mim + theme(plot.margin = unit(c(5.1, 2, 4.1, 1), "pt"))),
                      ggplotGrob(p.sar + theme(plot.margin = unit(c(5.1, 4.1, 4.1, 1), "pt"))))

ggsave("BdiversityBoxplot_all3_BrayCurtis_anova.png", p.all,
       scale = 1, width =25, height = 10, units = "cm",
       dpi = 600, limitsize = TRUE)

