# Make scatterplots of ASV prevalence vs ASV abundance 

# input:
#     Data frames (CSVs) made with DADA2_filterAndNorm.R
#       ASV table with samples in rows and ASVs in columns (normalized data)
#       TAX table with ASVs in rows and Taxonomic levels in columns
#       Sample data table with samples in rows and sample attributes in columns


setwd("~/")

study <- "TSPfinal"

################### IMPORT #####################################################
data.asv.norm <- read.csv(sprintf("%s_DADA2_Fraction.ASVtable.csv", study), row.names=1)
data.tax.norm <- read.csv(sprintf("%s_DADA2_Fraction.TAXtable.csv", study), row.names=1)
data.meta.norm <- read.csv(sprintf("%s_DADA2_Fraction.SAMPLEtable.csv", study), row.names=1)

# For loops with ggplot are annoying, so just run the stuff below 3 times, once with each species
species <- "dumicola"
#species <- "mimosarum"
#species <- "sarasinorum"

#for (species in c("dumicola", "mimosarum", "sarasinorum")) {
  data.meta.norm <- data.meta.norm[data.meta.norm$species == species, ]
  data.asv.norm <- data.asv.norm[, colnames(data.asv.norm) %in% rownames(data.meta.norm)]
  n.spec <- nrow(data.meta.norm)
  
  # find prevalence
  data.asv.pa <- data.asv.norm
  base.prev.thresh <- 0.0001
  base.prev.perc <- "0.01%"
  data.asv.pa[data.asv.pa < base.prev.thresh] <- 0 
  data.asv.pa[data.asv.pa > 0] <- 1
  rowSums(data.asv.pa)
  
  # get average abundance
  data.asv.abund <- rowSums(data.asv.norm)/ncol(data.asv.norm)
  
  # collect prevalence and abundance in same df, and get genus in rownames
  data.prev.abund <- data.frame(prevalence = rowSums(data.asv.pa),
                                abundance = data.asv.abund)
  rownames(data.prev.abund) <- paste0(data.tax.norm$Genus, "\n(", rownames(data.asv.norm), ")")
  
  # how many individual spiders with certain abundance of each asv
  individual.abundance.treshold <- 0.30
  dom.num <- c()
  for (i in 1:nrow(data.asv.norm)) {
    dom.num <- c(dom.num, sum(data.asv.norm[i, ] >= individual.abundance.treshold))
  }
  
  data.prev.abund$dom <- dom.num
  
  # make color vector
  prev.thresh <- "50"
  prev.threshold <- ncol(data.asv.norm)/2  # 50% prevalence
  #abundance.threshold <- 0.001
  ind.abund.threshold <- 3
  
  # individual abundance threshold
  hej.c <- ifelse(data.prev.abund$prevalence >= prev.threshold & data.prev.abund$dom >= ind.abund.threshold, "darkorchid", ifelse(data.prev.abund$prevalence > prev.threshold, "red", ifelse(data.prev.abund$dom >= ind.abund.threshold, "blue", "black"))) 

  data.prev.abund$col <- hej.c
  
  
  # make labels for the colored ones
  hej.l <- c()
  for (i in 1:nrow(data.prev.abund)) {
    if (data.prev.abund[i, "col"] != "black") {
      hej.l <- c(hej.l, rownames(data.prev.abund)[i])
    } else {
      hej.l <- c(hej.l, NA)
    }
  }
  
  data.prev.abund$lab <- hej.l
  #data.prev.abund$lab <- rownames(data.prev.abund)
  
  # only keep asvs with prevalence > 1
  data.prev.abund <- data.prev.abund[data.prev.abund$prevalence > 1, ]
  
  # change prevalence to percent of total n
  data.prev.abund$prevalence <- (data.prev.abund$prevalence/n.spec) * 100
  data.prev.abund$abundance <- data.prev.abund$abundance*100
  
  library(tidyverse)
  library(ggrepel)
  
  data.prev.abund$dom[data.prev.abund$dom < 3] <- NA
  #max(data.prev.abund$dom, na.rm=TRUE)
  library(tidyverse)
  library(ggrepel)
  
  my.title <- ifelse(species == "mimosarum", expression(paste(italic("S. mimosarum"), "(n=60)")), ifelse(species == "dumicola", expression(paste(italic("S. dumicola"), "(n=58)")),  expression(paste(italic("S. sarasinorum"), "(n=98)"))))
  if (species == "dumicola") {
    dum.p <- ggplot(data.prev.abund, aes(abundance, prevalence)) +
      theme_bw() +
      geom_hline(yintercept=50, linetype="dashed", color = "blue") +
      geom_point(aes(fill=dom), color= "black", pch = 21, size=3) +
      scale_fill_gradient(low = "yellow", high = "red", space = "Lab", 
                          na.value = "black", guide = "colourbar", aesthetics = "fill", 
                          limits = c(3, 20), name = "Individuals\ndominated") +
      ylab("Prevalence (% of individuals)") +
      xlab("Average Relative Abundance (%)") +
      scale_y_continuous(limits = c(0, 100)) + # y-axis limits
      geom_text_repel(aes(label=lab), hjust=0, vjust=0, size=3, xlim = c(3, NA)) +
      ggtitle(my.title) +
      theme(legend.position = "none")
  } else if (species == "mimosarum") {
    mim.p <- ggplot(data.prev.abund, aes(abundance, prevalence)) +
      theme_bw() +
      geom_hline(yintercept=50, linetype="dashed", color = "blue") +
      geom_point(aes(fill=dom), color= "black", pch = 21, size=3) +
      scale_fill_gradient(low = "yellow", high = "red", space = "Lab", 
                          na.value = "black", guide = "colourbar", aesthetics = "fill", 
                          limits = c(3, 20), name = "Individuals\ndominated") +
      ylab("Prevalence (% of individuals)") +
      xlab("Average Relative Abundance (%)") +
      scale_y_continuous(limits = c(0, 100)) + # y-axis limits
      geom_text_repel(aes(label=lab), hjust=0, vjust=0, size=3, xlim = c(3, NA)) +
      ggtitle(my.title) +
      theme(legend.position = "none") +
      theme(axis.title.y = element_blank())
  } else if (species == "sarasinorum") {
    sar.p <- ggplot(data.prev.abund, aes(abundance, prevalence)) +
      theme_bw() +
      geom_hline(yintercept=50, linetype="dashed", color = "blue") +
      geom_point(aes(fill=dom), color= "black", pch = 21, size=3) +
      scale_fill_gradient(low = "yellow", high = "red", space = "Lab", 
                            na.value = "black", guide = "colourbar", aesthetics = "fill", 
                            limits = c(3, 20), name = "Individuals\ndominated") +
      ylab("Prevalence (% of individuals)") +
      xlab("Average Relative Abundance (%)") +
      scale_y_continuous(limits = c(0, 100)) + # y-axis limits
      geom_text_repel(aes(label=lab), hjust=0, vjust=0, size=3, xlim = c(2, NA)) +
      ggtitle(my.title) +
      theme(axis.title.y = element_blank())
  }
#}

  
  
  
  
  
################ collect plots and export ######################################

library(gridExtra)
p3 <- grid.arrange(cbind(ggplotGrob(dum.p + theme(plot.margin = unit(c(5.1, 2, 4.1, 5.1), "pt")) +
                                      scale_x_continuous(breaks=seq(0, 30, 5))), 
                         ggplotGrob(mim.p + theme(plot.margin = unit(c(5.1, 2, 4.1, 1), "pt"), 
                                                  axis.text.y = element_blank(), 
                                                  axis.ticks.y = element_blank())), 
                         ggplotGrob(sar.p + theme(plot.margin = unit(c(5.1, 4.1, 4.1, 1), "pt"), 
                                                  axis.text.y = element_blank(), 
                                                  axis.ticks.y = element_blank()))))


ggsave("prev_vs_abund_domColor_all3_uglylabels_new.pdf", p3,
              scale = 1, width = 40, height = 10, units = "cm",
              dpi = 300, limitsize = TRUE)
       
