# Make heatmap indicating how many samples contain each ASV at certain abundances.

# input:
#     Data frames (CSVs) made with DADA2_to_Phyoloseq.R
#       ASV table with samples in rows and ASVs in columns
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns

setwd("~/ /")
study <- "TSPfinal"

library(vegan)
library(stringr)
library(ggplot2)
library(gridExtra)

################### IMPORT #####################################################
data.asv.norm.raw <- read.csv(sprintf("%s_DADA2_Fraction.ASVtable.csv", study), row.names=1)
data.tax.norm.raw <- read.csv(sprintf("%s_DADA2_Fraction.TAXtable.csv", study), row.names=1)
data.meta.norm.raw <- read.csv(sprintf("%s_DADA2_Fraction.SAMPLEtable.csv", study), row.names=1)

data.tax.norm.raw$ASV <- paste(data.tax.norm.raw$Genus, rownames(data.tax.norm.raw))

# count how many samples each ASV appear in (above 0.1% threshold)
data.asv.norm.present <- data.asv.norm.raw
data.asv.norm.present[data.asv.norm.present < 0.001] <- 0
data.asv.norm.present[data.asv.norm.present > 0] <- 1
ASV.count <- rowSums(data.asv.norm.present)

# insert bactrial genus and ASV count in rownames
sum(rownames(data.asv.norm.raw) == rownames(data.tax.norm.raw)) == length(rownames(data.asv.norm.raw))
rownames(data.asv.norm.raw) <- paste0(data.tax.norm.raw$Genus, " ", rownames(data.tax.norm.raw), " (", ASV.count, ")")

# find the most dominant ASVs across all 3 spider species
data.asv.norm.raw.t <- as.data.frame(t(data.asv.norm.raw))

breaks <- c(1, 0.8, 0.6, 0.4)
number <- 3

abundance.cut <- cut(as.vector(data.asv.norm.raw.t[, 1]), breaks, right=FALSE)
results <- as.data.frame(table(abundance.cut))
for (ASVs in rownames(data.asv.norm.raw)){
  abundance.cut <- cut(as.vector(data.asv.norm.raw.t[, ASVs]), breaks, right=FALSE)
  abundance.freq <- as.data.frame(table(abundance.cut))
  results <- cbind(results, abundance.freq$Freq)
}

rownames(results) <- results$abundance.cut
results$abundance.cut <- NULL
results$Freq <- NULL
colnames(results) <- rownames(data.asv.norm.raw)

# Get the most dominant ASVs across all spider species.
hej <- results[2, ] + results[3, ]
hej2 <- results[1, ] + results[2, ] + results[3, ]
results <- rbind(results, hej, hej2)
results <- results[,order(-results[(nrow(results)-1),], -results[nrow(results), ])]
resultsSmall <- results[1:number, 1:15]
resultsSmall <- resultsSmall[,order(-resultsSmall[nrow(resultsSmall),])]
ASVs.for.all <- colnames(resultsSmall)


# loop through spider species and save a plot for each
for (species in c("dumicola", "mimosarum", "sarasinorum")) {
  
  data.meta.norm <- data.meta.norm.raw[which(data.meta.norm.raw$species == species), ]
  data.asv.norm <- data.asv.norm.raw[, which(colnames(data.asv.norm.raw) %in% rownames(data.meta.norm))]
  data.tax.norm <- data.tax.norm.raw
  
  data.asv.norm.t <- as.data.frame(t(data.asv.norm))
  
  breaks <- c(1, 0.8, 0.6, 0.4)
  number <- 3
  
  abundance.cut <- cut(as.vector(data.asv.norm.t[, 1]), breaks, right=FALSE)
  results <- as.data.frame(table(abundance.cut))
  for (ASVs in rownames(data.asv.norm)){
    abundance.cut <- cut(as.vector(data.asv.norm.t[, ASVs]), breaks, right=FALSE)
    abundance.freq <- as.data.frame(table(abundance.cut))
    results <- cbind(results, abundance.freq$Freq)
  }
  
  rownames(results) <- results$abundance.cut
  results$abundance.cut <- NULL
  results$Freq <- NULL
  colnames(results) <- rownames(data.asv.norm)
  
  # order by most dominant ASVs for ll 3 species
  resultsSmall <- results[1:number, ASVs.for.all]
  
  # fix names for some cursive, some not
  df <- str_split(colnames(resultsSmall), " ", simplify = T)
  df <- trimws(df)  # remove trailing spaces
  df <- as.data.frame(df)
  colnames(df) <- c("name", "number", "count")
  # format field 1 (=ASV latin name) to be cursive, unless it is "other or "subject"
  exceptions <- c("other", "subject", "Weeksellaceae", "Rickettsiaceae", "Lachnospiraceae", "Spirochaetaceae", "Chlamydiales", "Orbaceae")  # these can be added on to if needed
  exceptions2 <- c("other", "subject")
  
  df$name <- ifelse(df$name %in%  exceptions,
                    paste0(df$name),
                    paste0("italic(", df$name, ")"))
  df$longname <- paste(df$name, df$number, df$count)
  df$longname <- ifelse(df$name %in%  exceptions2,
                        paste0(df$name),
                        str_replace_all(df$longname, " ", "~~"))
  colnames(resultsSmall) <- df$longname
  rownames(resultsSmall) <- c("[0.4, 0.6)", "[0.6, 0.8)", "[0.8, 1]")
  library(reshape2)
  library(ggplot2)
  dat <- as.matrix(resultsSmall)
  dat2 <- melt(dat)
  
  if (species == "mimosarum"){
    plot.title <- expression(italic("S. mimosarum"))
    p.mim <- ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
      geom_tile(aes(fill = value)) +
      geom_text(size = 3, aes(label = round(value, 1))) +  # numbers in plot
      theme_bw() +
      theme(text=element_text(family="Helvetica"), axis.title.y = element_blank()) +
      #scale_y_discrete(breaks=unique(dat2$Var2), labels=parse(text=as.character(unique(dat2$Var2)))) +
      scale_fill_gradient(low = "white", high = "steelblue", name=NULL, limits=c(0,20), breaks=c(0,5,10,15,20), guide=FALSE) +
      theme(axis.title.y = element_blank()) +
      xlab("ASV dominance in fraction per sample") +
      ggtitle(plot.title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.margin=margin(5,5,-1,5), legend.text = element_text(size=10)) +
      theme(axis.title.x = element_text(size = 12, margin=margin(10,0,0,0))) +
      theme(axis.text.x = element_text(hjust = 0.5, size = 7, color = 'black')) +
      theme(plot.margin=margin(5,0,5,0)) +
      theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
  }
  if (species == "sarasinorum") {
    plot.title <- expression(italic("S. sarasinorum"))
    p.sar <- ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
      geom_tile(aes(fill = value)) +
      geom_text(size = 3, aes(label = round(value, 1))) +  # numbers in plot
      theme_bw() +
      theme(text=element_text(family="Helvetica"), axis.title.y = element_blank()) +
      #scale_y_discrete(breaks=unique(dat2$Var2), labels=parse(text=as.character(unique(dat2$Var2)))) +
      scale_fill_gradient(low = "white", high = "steelblue", name="Number\nof samples", limits=c(0,20),breaks=c(0,5,10,15,20)) +
      theme(axis.title.y = element_blank()) +
      #xlab(sprintf("Relative abundance of ASV \n in S.%s", species)) +
      xlab("Relative abundance of ASV") +
      ggtitle(plot.title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      #theme(legend.position = "right") +
      theme(legend.margin=margin(5,5,-1,5), legend.text = element_text(size=10)) +
      #theme(axis.title.x = element_text(size = 12, margin=margin(10,0,0,0))) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_text(hjust = 0.5, size = 7, color = 'black')) +
      theme(plot.margin=margin(5,5,5,0)) +
      theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
  }
  if (species == "dumicola") {
    plot.title <- expression(italic("S. dumicola"))
    p.dum <- ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
      geom_tile(aes(fill = value)) +
      geom_text(size = 3, aes(label = round(value, 1))) +  # numbers in plot
      theme_bw() +
      theme(text=element_text(family="Helvetica")) +
      scale_y_discrete(breaks=unique(dat2$Var2), labels=parse(text=as.character(unique(dat2$Var2)))) +
      scale_fill_gradient(low = "white", high = "steelblue", name=NULL, limits=c(0,20), breaks=c(0,5,10,15,20), guide=FALSE) +
      theme(axis.title.y = element_blank()) +
      #xlab("Relative abundance of ASV") +
      ggtitle(plot.title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "top") +
      theme(legend.margin=margin(5,5,-1,5), legend.text = element_text(size=10)) +
      #theme(axis.title.x = element_text(size = 12, margin=margin(10,0,0,0))) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_text(hjust = 0.5, size = 7, color = 'black')) +
      theme(plot.margin=margin(5,0,5,5)) +
      #theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
      theme(axis.text.y = element_text(size = 10, color = 'black', face="italic"))
  }
}

# collect plots in same image
p.dum <- ggplotGrob(p.dum)
p.mim <- ggplotGrob(p.mim)
p.sar <- ggplotGrob(p.sar)
p.all <- gtable_cbind(p.dum, p.mim, p.sar)
plot(p.all) # this does not work in Rstudio, you have to export to see the final plot

# export plot
ggsave("Dominance_plot.jpg", p.all,
       scale = 1,
       dpi = 600, limitsize = TRUE)


