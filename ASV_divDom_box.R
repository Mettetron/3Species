# Make Boxplots of ASV richness, alpha diversity and dominance

# input:
#     Data frames (CSVs) possibly made with DADA2_filterAndNorm.R
#       Normalized ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns

setwd("~/") 
study <- "TSPfinal"

################### IMPORT #####################################################

data.asv.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Subsampl3000.SAMPLEtable.csv", study), row.names=1)


################### COLLECT METRICS ############################################
threshold <- 1  # this can be changed if you want an abundance threshold
library(vegan)
metrics <- c("n", "OTUs_total", "shannon_mean", "shannon_SD", "shannon_median", "simpson_mean", "simpson_SD", "dmn_mean", "dmn_SD", "dmn_median", "OTUs_per_sample", "OTUs_per_sample_SD", "OTUs_per_sample_median")
metrics.df <- as.data.frame(metrics)
asv.mat <- as.matrix(data.asv.norm)
asv.mat <- asv.mat[, colSums(asv.mat) > 0]
data.meta <- data.meta.norm
data.meta <- data.meta[colnames(asv.mat), ]
for (spec in unique(data.meta$species)) {
    data.meta.spec <- data.meta[data.meta$species == spec, ]
    asv.mat.spec <- asv.mat[, rownames(data.meta.spec)]
    ncol(asv.mat.spec)
    if (is.null(ncol(asv.mat.spec))) {
      sample.names <- rownames(data.meta.spec)
      asv.mat.spec <- asv.mat.spec[asv.mat.spec > 0]
    } else {
      asv.mat.spec <- as.matrix(asv.mat.spec[rowSums(asv.mat.spec) > 0, ])
      sample.names <- colnames(asv.mat.spec)
    }
    results.df <- data.frame("shannon"=c(1), "simpson"=c(1), 
                             "dmn"=c(1), "asv_numbers_thresholded"=c(1))
    results.df <- as.data.frame(t(results.df))
    for (n in sample.names) {
      # get thresholded asvs for one sample at a time
      if (is.null(ncol(asv.mat.spec))) {
        s <- asv.mat.spec
      } else {
        s <- asv.mat.spec[, n]
      }
      s2 <- s[s >= threshold]
      # calculate metrics
      shannon <- diversity(s2, index="shannon")
      simpson <- diversity(s2, index="simpson")
      s.ordered <- rev(sort(s2))
      dmn <- (s.ordered[1]+s.ordered[2])/3000
      asv.num <- length(s2)
      results.df[, n] <- c(shannon, simpson, dmn, asv.num)
    }
    results.df$V1 <- NULL
    sds <- apply(results.df, 1, sd, na.rm = TRUE)
    means <- apply(results.df, 1, mean, na.rm = TRUE)
    medians <- apply(results.df, 1, median, na.rm = TRUE)
    results.df$Mean <- means
    results.df$SD <- sds
    results.df$Median <- medians
    results.df.short <- results.df[, c("Mean", "SD", "Median")]
    assign(paste0(study,"_", spec, "_samples"), results.df)
    assign(paste0(study, "_", spec, "_overview"), results.df.short)
}


# get ready for boxplots
data.meta <- data.meta.norm
study.species <- paste0(study, "_", unique(data.meta$species))

################### DMN PLOT ################################################### 
m <- "dmn"
chart.data <- data.frame(m = c(1), "study" = c(1))
for (s in study.species) {
  dat <- get(paste0(s, "_samples"))
  m.vec <- dat[m, 1:(ncol(dat)-2)]
  m.dat <- rbind(m.vec, rep(s, (ncol(dat)-2)))
  m.dat <- as.data.frame(t(m.dat))
  colnames(m.dat) <- c("m", "study")
  chart.data <- rbind(chart.data, m.dat)
}
chart.data <- chart.data[2:nrow(chart.data), ]
chart.data$m <- as.numeric(as.character(chart.data$m))
study.species.list <- strsplit(as.character(chart.data$study), '_')
chart.data$species <- paste0("Stegodyphus", "_",sapply(study.species.list, "[[", 2))
chart.data$species <- factor(chart.data$species,
                             levels = c("Stegodyphus_dumicola", "Stegodyphus_mimosarum", "Stegodyphus_sarasinorum"), ordered = TRUE)

# make significance bars
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(1.08, 1.1, 1.1, 1.08))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(1.08, 1.1, 1.1, 1.08))
df3 <- data.frame(a = c(1.1, 1.1, 2.9, 2.9), b = c(1.15, 1.17, 1.17, 1.15))

# t-test
hej <- pairwise.t.test(chart.data$m, chart.data$species,
                p.adjust.method = "BH")
hej2 <- hej$p.value

dum.mim <- paste0("p = ", as.character(formatC(hej2[1,1], format = "e", digits = 2)))
mim.sar <- paste0("p = ", as.character(formatC(hej2[2,2], format = "e", digits = 2)))
dum.sar <- paste0("p = ", as.character(formatC(hej2[2,1], format = "e", digits = 2)))

ggplot(chart.data, aes(x=species, y=m)) +
  geom_boxplot(color="black", outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle("C") +
  theme(plot.title = element_text(size=16)) +
  scale_y_continuous(limits=c(0,1.2), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(panel.grid.major.x = element_blank()) +  # remove vetical grid lines
  ylab(label="DMN") +  # insert yaxis label
  theme(axis.title.y = element_text(size = 14, color="black")) +  # modify size of yaxis label
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 9.3, color="black")) + # changes to x-axis text
  theme(axis.text.y = element_text(color="black")) +
  scale_x_discrete(labels=parse(text=c("italic(S.~dumicola)", "italic(S.~mimosarum)", "italic(S.~sarasinorum)"))) + # cursive x-axis 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=2, stroke = 1, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.12, label = dum.mim, size = 3) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 1.12, label = mim.sar, size = 3) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.19, label = dum.sar, size = 3)

ggsave(sprintf("DMN_boxplot_3sp_signifBars_pval.jpg"),
       scale = 1, width = 10, height = 10, units = "cm",
       dpi = 600, limitsize = TRUE)


################### SHANNON PLOT ###############################################
m <- "shannon"
chart.data <- data.frame(m = c(1), "study" = c(1))
for (s in study.species) {
  dat <- get(paste0(s, "_Gini_samples"))
  m.vec <- dat[m, 1:(ncol(dat)-2)]
  m.dat <- rbind(m.vec, rep(s, (ncol(dat)-2)))
  m.dat <- as.data.frame(t(m.dat))
  colnames(m.dat) <- c("m", "study")
  chart.data <- rbind(chart.data, m.dat)
}
chart.data <- chart.data[2:nrow(chart.data), ]
chart.data$m <- as.numeric(as.character(chart.data$m))
study.species.list <- strsplit(as.character(chart.data$study), '_')
chart.data$species <- paste0("Stegodyphus", "_",sapply(study.species.list, "[[", 2)) 
chart.data$species <- factor(chart.data$species,
                             levels = c("Stegodyphus_dumicola", "Stegodyphus_mimosarum", "Stegodyphus_sarasinorum"), ordered = TRUE)


# make significance bars
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(4.5, 4.6, 4.6, 4.5))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(4.5, 4.6, 4.6, 4.5))
df3 <- data.frame(a = c(1.1, 1.1, 2.9, 2.9), b = c(4.8, 4.9, 4.9, 4.8))

# t-test
hej <- pairwise.t.test(chart.data$m, chart.data$species,
                       p.adjust.method = "BH")
hej2 <- hej$p.value

dum.mim <- paste0("p = ", as.character(formatC(hej2[1,1], format = "e", digits = 2)))
mim.sar <- paste0("p = ", as.character(formatC(hej2[2,2], format = "e", digits = 2)))
dum.sar <- paste0("p = ", as.character(formatC(hej2[2,1], format = "e", digits = 2)))

ggplot(chart.data, aes(x=species, y=m)) +
  geom_boxplot(color="black", outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle("B") +
  theme(plot.title = element_text(size=16)) +
  scale_y_continuous(limits=c(0,5.2), breaks=c(0, 1, 2, 3, 4, 5)) +
  theme(panel.grid.major.x = element_blank()) +  # remove vetical grid lines
  ylab(label="Shannon Index") +  # insert yaxis label
  theme(axis.title.y = element_text(size = 14, color="black")) +  # modify size of yaxis label
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 9.3, color="black")) + # changes to x-axis text
  theme(axis.text.y = element_text(color="black")) +
  scale_x_discrete(labels=parse(text=c("italic(S.~dumicola)", "italic(S.~mimosarum)", "italic(S.~sarasinorum)"))) + # cursive x-axis 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=2, stroke = 1, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 4.7, label = dum.mim, size = 3) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 4.7, label = mim.sar, size = 3) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2, y = 5, label = dum.sar, size = 3)

ggsave(sprintf("Shannon_boxplot_3sp_signifBars.jpg"),
       scale = 1, width = 10, height = 10, units = "cm",
       dpi = 600, limitsize = TRUE)


################### RICHNESS PLOT ##############################################
m <- "asv_numbers_thresholded"
chart.data <- data.frame(m = c(1), "study" = c(1))
for (s in study.species) {
  dat <- get(paste0(s, "_Gini_samples"))
  m.vec <- dat[m, 1:(ncol(dat)-2)]
  m.dat <- rbind(m.vec, rep(s, (ncol(dat)-2)))
  m.dat <- as.data.frame(t(m.dat))
  colnames(m.dat) <- c("m", "study")
  chart.data <- rbind(chart.data, m.dat)
}
chart.data <- chart.data[2:nrow(chart.data), ]
chart.data$m <- as.numeric(as.character(chart.data$m))
study.species.list <- strsplit(as.character(chart.data$study), '_')
chart.data$species <- paste0("Stegodyphus", "_",sapply(study.species.list, "[[", 2))
chart.data$species <- factor(chart.data$species,
                             levels = c("Stegodyphus_dumicola", "Stegodyphus_mimosarum", "Stegodyphus_sarasinorum"), ordered = TRUE)

# make significance bars
df1 <- data.frame(a = c(1.1, 1.1, 1.9, 1.9), b = c(270, 275, 275, 270))
df2 <- data.frame(a = c(2.1, 2.1, 2.9, 2.9), b = c(270, 275, 275, 270))
df3 <- data.frame(a = c(1.1, 1.1, 2.9, 2.9), b = c(285, 290, 290, 285))

# t-test
hej <- pairwise.t.test(chart.data$m, chart.data$species,
                       p.adjust.method = "BH")
hej2 <- hej$p.value

dum.mim <- paste0("p = ", as.character(formatC(hej2[1,1], format = "e", digits = 2)))
mim.sar <- paste0("p = ", as.character(formatC(hej2[2,2], format = "e", digits = 2)))
dum.sar <- paste0("p = ", as.character(formatC(hej2[2,1], format = "e", digits = 2)))

ggplot(chart.data, aes(x=species, y=m)) +
  geom_boxplot(color="black", outlier.shape=1) + 
  theme_bw(base_family="Helvetica") +
  ggtitle("A") +
  theme(plot.title = element_text(size=16)) +
  theme(panel.grid.major.x = element_blank()) +  # remove vetical grid lines
  ylab(label="Observed ASV numbers") +  # insert yaxis label
  theme(axis.title.y = element_text(size = 14, color="black")) +  # modify size of yaxis label
  theme(axis.title.x = element_blank()) +  # remove xaxis label
  theme(axis.text.x = element_text(size = 9.3, color="black")) + # changes to x-axis text
  theme(axis.text.y = element_text(color="black")) +
  scale_x_discrete(labels=parse(text=c("italic(S.~dumicola)", "italic(S.~mimosarum)", "italic(S.~sarasinorum)"))) + # cursive x-axis 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=4, size=2, stroke = 1, show.legend = FALSE) +  # add mean to plot as "X"
  geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 280, label = dum.mim, size = 3) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 280, label = mim.sar, size = 3) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2, y = 295, label = dum.sar, size = 3)

ggsave(sprintf("ASVnum_boxplot_3sp_signifBars.jpg"),
       scale = 1, width = 10, height = 10, units = "cm",
       dpi = 600, limitsize = TRUE)
