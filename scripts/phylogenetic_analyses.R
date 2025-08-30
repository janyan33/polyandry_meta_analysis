library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(DHARMa)
library(metafor)
library(rotl)
library(ape)
library(ggcorrplot)

# We used the function tnrs_match_names() from the `rotl` package to build a phylogenetic tree
# with taxonomic data from Open Tree Taxonomy

# Tree for fecundity analyses was created and saved on August 30th, 2025
load("data/fecund_tree.Rdata")

# Plotting fecundity tree
plot(tree, cex = 0.5, label.offset = 0.25, no.margin = TRUE)

### Checking that tree is complete
#tree$tip.label <- gsub("_", " ", tree$tip.label)

data_fecund <- read.csv("data/fecund_data.csv")

# listed in tree but not in data
setdiff(as.character(tree$tip.label),
          as.character(data_fecund$species))

# listed in data but not in tree
setdiff(as.character(data_fecund$species),
           as.character(tree$tip.label))

# Compute branch lengths of tree
phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)

is.ultrametric(phylo_branch)

phylo_cor <- vcv(phylo_branch, cor = T)

# visual exploration of the phylogenetic correlation matrix
ggcorrplot::ggcorrplot(phylo_cor, sig.level = 0.05, lab_size = 1,
                       p.mat = NULL,insig = c("pch", "blank"),
                       pch = 1, pch.col = "black", pch.cex = 1, tl.cex = 1.5) +
  theme(axis.text.x = element_text(size = 5, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_text(size = 5, margin = margin(0, -2, 0, 0)),
        panel.grid.minor = element_line(size = 3)) +
  geom_tile(fill = "white") +
  geom_tile(height = 0.8, width = 0.8) +
  scale_fill_gradient2(low = "#E69F00",mid = "white", high = "#56B4E9",
                       midpoint = 0.5, breaks = c(0, 1),
                       limit = c(0,1)) + labs(fill = "Correlation")
