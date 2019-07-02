# plot tree based on phylogeny calculated with Mr. Bayes

# input:
#   tree file in newick format (example file provided on GitHub)
#   csv with name and formatting information, tip.label column must match sequence IDs in tree
#     (example file provided on GitHub)  
#   you have to manually input posterior probabilities from e.g. Mr. Bayes output file

# directory containing tree files
setwd("~/ /")

# import tree in newick format
library(ape)
tree <- read.tree("Borrelia_3s_shortNames.tree")  # newick format tree, exported from ARB with only sequence IDs

# import csv with tree name and formatting information, 
# tip.label column must match sequence IDs in tree
df.annotation <- read.csv2("bortreenames.csv", row.names = 1)

# To add posterior probabilities from other source - takes a few iterations to get it right. plot tree first.
shown <- 1:tree$Nnode  # run first, to see which is which
shown <- rep(100,tree$Nnode)  # if no other info, number is 100
shown[c(4)] <- NA  # get rid of number at root and at ASV branches
shown[c(6,7)] <- c(63, 69)  # special placement for numbers not = 100 (based on Mr. Bayes output)

# plot
tip.linetypes <- c(1,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)  # add broken lines for ASV-insertions
library(ggtree)
tree <- groupClade(tree, .node=13)
library(ggplot2) 
ggplot(tree, aes(color=group, linetype=I(tip.linetypes))) %<+% df.annotation + # this plots the short name tree, but uses the long names from df.annotation
  geom_tree(size=0.5) +
  theme_tree() +
  geom_tiplab(aes(colour = group, label = name), parse = T, hjust=-0.02, family="Helvetica") +  # tip labels  added and parsed
  scale_x_continuous(limits = c(0, 0.8)) +  # plot dimensions
  geom_nodelab(aes(label=c(1:length(tree$tip.label), shown)), vjust=-0.4, hjust=1.25, size=2.5) + # id, size and location of posterior probablilites
  scale_color_manual(values=c("black","blue")) +  # choose colors
  geom_treescale(x=0.1, y=-0.5, offset=0.3, width=0.1)  # add small scale, and choose position
  
ggsave("BorreliaTree.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 28, height = length(tree$tip.label), units = "cm",
       dpi = 300, limitsize = TRUE)

