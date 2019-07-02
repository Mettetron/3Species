# make boxplot(s) comparing beta diversities between host populations and nests 

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#         it is recommended to use subsampled data for calculation of distance metrics
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns

setwd("~/ /") 
study <- "TSPfinal"

################### IMPORT #####################################################

data.asv.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.SAMPLEtable.csv", study), row.names=1)

################### SINGLE SPECIES #############################################

# species can be changed here to run script for different species.
species <- "mimosarum"
#species <- "dumicola"
#species <- "sarasinorum"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.asv.norm <- data.asv.norm[, which(colnames(data.asv.norm) %in% rownames(data.meta.norm))]

# also change color for plot
# species colors: '#7570b3' = mim, '#d95f02' = dum, '#1b9e77' = sar
#species_color <- ifelse(species=="dumicola", '#d95f02', ifelse(species=="mimosarum", '#7570b3', '#1b9e77'))
species_color <- "black" # boring colors for publication

# check that things line up
sum(rownames(data.meta.norm) == colnames(data.asv.norm)) == length(rownames(data.meta.norm))

################### GET DIVERSITIES ############################################
# transpose
data.asv.norm.t <- as.data.frame(t(data.asv.norm))

# get Bray-Curtis
library(vegan)
beta.dist <- vegdist(data.asv.norm.t, method="bray", binary=FALSE)
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
sample.to.pop <- as.data.frame(cbind(rownames(data.meta.norm), as.character(data.meta.norm$population)))
colnames(sample.to.pop) <- c("sample", "population")
# row pops
new.frame <- as.data.frame(beta.dist.df.unique[, "row"])
new.frame[] <- lapply(new.frame, function(x) sample.to.pop$population[match(x, sample.to.pop$sample)])
colnames(new.frame) <- "pop1"
beta.dist.df.unique$pop1 <- new.frame$pop1
# col pops
new.frame <- as.data.frame(beta.dist.df.unique[, "col"])
new.frame[] <- lapply(new.frame, function(x) sample.to.pop$population[match(x, sample.to.pop$sample)])
colnames(new.frame) <- "pop2"
beta.dist.df.unique$pop2 <- new.frame$pop2

# insert nest
sample.to.nest <- as.data.frame(cbind(rownames(data.meta.norm), as.character(data.meta.norm$nest)))
colnames(sample.to.nest) <- c("sample", "nest")
# row pops
new.frame <- as.data.frame(beta.dist.df.unique[, "row"])
new.frame[] <- lapply(new.frame, function(x) sample.to.nest$nest[match(x, sample.to.nest$sample)])
colnames(new.frame) <- "nest1"
beta.dist.df.unique$nest1 <- new.frame$nest1
# row pops
new.frame <- as.data.frame(beta.dist.df.unique[, "col"])
new.frame[] <- lapply(new.frame, function(x) sample.to.nest$nest[match(x, sample.to.nest$sample)])
colnames(new.frame) <- "nest2"
beta.dist.df.unique$nest2 <- new.frame$nest2

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
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(1.08, 1.1, 1.1, 1.08))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(1.08, 1.1, 1.1, 1.08))

# set significance level (calculated below plot)
nest.spop <- "**"  # it's <2e-16
spop.dpop <- ifelse(species == "dumicola", "*", "**") # it's <2e-16 if sar or mim, 0.01 if dum
plot.title <- ifelse(species == "dumicola", expression(paste(italic("S. dumicola"), " beta diversity")), ifelse(species == "mimosarum", expression(paste(italic("S. mimosarum"), " beta diversity")), expression(paste(italic("S. sarasinorum"), " beta diversity"))))

ggplot(anova.data, aes(x=group, y=Bdiv)) +
  geom_boxplot(color=species_color, outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle(plot.title) +
  theme(plot.title = element_text(size=16)) +
  scale_y_continuous(limits=c(0,1.12), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(panel.grid.major = element_blank()) +  # remove vertical grid lines
  theme(panel.grid.minor = element_blank()) +  # remove horizontal grid lines
  ylab(label="Bray-Curtis dissimilarity") +  # insert yaxis label
  theme(axis.title.y = element_text(size = 14, color="black")) +  # modify size of yaxis label
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 14, color="black")) + # changes to x-axis text
  theme(axis.text.y = element_text(color="black")) +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=3.5, stroke = 1.5, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.1, label = nest.spop, size = 8) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 1.1, label = spop.dpop, size = 8)

ggsave(sprintf("diversityBoxplot_%s.jpg", species),
       scale = 1, width = 20, height = 10, units = "cm",
       dpi = 600, limitsize = TRUE)

################# ANOVA TO COMPARE GROUPS ######################################
# method from: http://www.sthda.com/english/wiki/one-way-anova-test-in-r
same_nest <- cbind(as.numeric(as.character(bd.same.nest$value)), rep("s_nest", nrow(bd.same.nest)))
same_pop <- cbind(as.numeric(as.character(bd.same.pop$value)), rep("s_pop", nrow(bd.same.pop)))
diff_pop <- cbind(as.numeric(as.character(bd.diff.pop$value)), rep("d_pop", nrow(bd.diff.pop)))
anova.data <- as.data.frame(rbind(same_nest, same_pop, diff_pop))   
colnames(anova.data) <- c("Bdiv", "group")
anova.data$Bdiv <- as.numeric(as.character(anova.data$Bdiv))

# averages
average.same.nest <- mean(bd.same.nest$value)
average.same.nest
average.same.pop <- mean(bd.same.pop$value)
average.same.pop
average.diff.pop <- mean(bd.diff.pop$value)
average.diff.pop

levels(anova.data$group)
anova.data$group <- ordered(anova.data$group, levels = c("d_pop", "s_nest", "s_pop" ))
library(dplyr)
group_by(anova.data, group) %>%
  summarise(
    count = n(),
    mean = mean(Bdiv, na.rm = TRUE),
    sd = sd(Bdiv, na.rm = TRUE)
  )

bdiv.anova <- aov(Bdiv ~ group, data = anova.data)
summary(bdiv.anova)
TukeyHSD(bdiv.anova)

pairwise.t.test(anova.data$Bdiv, anova.data$group,
                p.adjust.method = "BH")
#### sar

# > # averages
#   > average.same.nest <- mean(bd.same.nest$value)
# > average.same.nest
# [1] 0.3334707
# > average.same.pop <- mean(bd.same.pop$value)
# > average.same.pop
# [1] 0.7178737
# > average.diff.pop <- mean(bd.diff.pop$value)
# > average.diff.pop
# [1] 0.9168016
# > 
#   > levels(anova.data$group)
# [1] "d_pop"  "s_nest" "s_pop" 
# > anova.data$group <- ordered(anova.data$group, levels = c("d_pop", "s_nest", "s_pop" ))
# > library(dplyr)
# > group_by(anova.data, group) %>%
#   +   summarise(
#     +     count = n(),
#     +     mean = mean(Bdiv, na.rm = TRUE),
#     +     sd = sd(Bdiv, na.rm = TRUE)
#     +   )
# # A tibble: 3 x 4
# group  count  mean    sd
# <ord>  <int> <dbl> <dbl>
#   1 d_pop   4006 0.917 0.168
# 2 s_nest   131 0.333 0.274
# 3 s_pop    747 0.718 0.291
# > 
#   > bdiv.anova <- aov(Bdiv ~ group, data = anova.data)
# > summary(bdiv.anova)
# Df Sum Sq Mean Sq F value Pr(>F)    
# group          2  63.77   31.88   834.6 <2e-16 ***
#   Residuals   4881 186.47    0.04                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > TukeyHSD(bdiv.anova)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Bdiv ~ group, data = anova.data)
# 
# $group
# diff        lwr        upr p adj
# s_nest-d_pop -0.5833309 -0.6240158 -0.5426460     0
# s_pop-d_pop  -0.1989279 -0.2171900 -0.1806659     0
# s_pop-s_nest  0.3844030  0.3409987  0.4278073     0
# 
# > 
#   > pairwise.t.test(anova.data$Bdiv, anova.data$group,
#                     +                 p.adjust.method = "BH")
# 
# Pairwise comparisons using t tests with pooled SD 
# 
# data:  anova.data$Bdiv and anova.data$group 
# 
# d_pop  s_nest
# s_nest <2e-16 -     
#   s_pop  <2e-16 <2e-16
# 
# P value adjustment method: BH 
### mim

# > # averages
#   > average.same.nest <- mean(bd.same.nest$value)
# > average.same.nest
# [1] 0.3836474
# > average.same.pop <- mean(bd.same.pop$value)
# > average.same.pop
# [1] 0.764483
# > average.diff.pop <- mean(bd.diff.pop$value)
# > average.diff.pop
# [1] 0.8922544
# > 
#   > levels(anova.data$group)
# [1] "d_pop"  "s_nest" "s_pop" 
# > anova.data$group <- ordered(anova.data$group, levels = c("d_pop", "s_nest", "s_pop" ))
# > library(dplyr)
# > group_by(anova.data, group) %>%
#   +   summarise(
#     +     count = n(),
#     +     mean = mean(Bdiv, na.rm = TRUE),
#     +     sd = sd(Bdiv, na.rm = TRUE)
#     +   )
# # A tibble: 3 x 4
# group  count  mean    sd
# <ord>  <int> <dbl> <dbl>
#   1 d_pop   1427 0.892 0.206
# 2 s_nest    52 0.384 0.316
# 3 s_pop    343 0.764 0.293
# > 
#   > bdiv.anova <- aov(Bdiv ~ group, data = anova.data)
# > summary(bdiv.anova)
# Df Sum Sq Mean Sq F value Pr(>F)    
# group          2  16.34   8.170   156.5 <2e-16 ***
#   Residuals   1819  94.99   0.052                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > TukeyHSD(bdiv.anova)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Bdiv ~ group, data = anova.data)
# 
# $group
# diff        lwr         upr p adj
# s_nest-d_pop -0.5086069 -0.5842830 -0.43293086     0
# s_pop-d_pop  -0.1277714 -0.1600055 -0.09553727     0
# s_pop-s_nest  0.3808356  0.3010658  0.46060527     0
# 
# > 
#   > pairwise.t.test(anova.data$Bdiv, anova.data$group,
#                     +                 p.adjust.method = "BH")
# 
# Pairwise comparisons using t tests with pooled SD 
# 
# data:  anova.data$Bdiv and anova.data$group 
# 
# d_pop  s_nest
# s_nest <2e-16 -     
#   s_pop  <2e-16 <2e-16
# 
# P value adjustment method: BH 

### dum

# > # averages
#   > average.same.nest <- mean(bd.same.nest$value)
# > average.same.nest
# [1] 0.2087193
# > average.same.pop <- mean(bd.same.pop$value)
# > average.same.pop
# [1] 0.7819878
# > average.diff.pop <- mean(bd.diff.pop$value)
# > average.diff.pop
# [1] 0.8336191
# > 
#   > levels(anova.data$group)
# [1] "d_pop"  "s_nest" "s_pop" 
# > anova.data$group <- ordered(anova.data$group, levels = c("d_pop", "s_nest", "s_pop" ))
# > library(dplyr)
# > group_by(anova.data, group) %>%
#   +   summarise(
#     +     count = n(),
#     +     mean = mean(Bdiv, na.rm = TRUE),
#     +     sd = sd(Bdiv, na.rm = TRUE)
#     +   )
# # A tibble: 3 x 4
# group  count  mean    sd
# <ord>  <int> <dbl> <dbl>
#   1 d_pop   1325 0.834 0.319
# 2 s_nest    38 0.209 0.269
# 3 s_pop    328 0.782 0.357
# > 
#   > bdiv.anova <- aov(Bdiv ~ group, data = anova.data)
# > summary(bdiv.anova)
# Df Sum Sq Mean Sq F value Pr(>F)    
# group          2  14.73   7.367   69.39 <2e-16 ***
#   Residuals   1688 179.21   0.106                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > TukeyHSD(bdiv.anova)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Bdiv ~ group, data = anova.data)
# 
# $group
# diff         lwr          upr     p adj
# s_nest-d_pop -0.62489982 -0.75065459 -0.499145053 0.0000000
# s_pop-d_pop  -0.05163131 -0.09876895 -0.004493675 0.0277144
# s_pop-s_nest  0.57326851  0.44229360  0.704243412 0.0000000
# 
# > 
#   > pairwise.t.test(anova.data$Bdiv, anova.data$group,
#                     +                 p.adjust.method = "BH")
# 
# Pairwise comparisons using t tests with pooled SD 
# 
# data:  anova.data$Bdiv and anova.data$group 
# 
# d_pop  s_nest
# s_nest <2e-16 -     
#   s_pop  0.01   <2e-16
# 
# P value adjustment method: BH 