library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(DHARMa)
library(metafor)
library(rotl)
library(ape)
library(ggcorrplot)

# We used the function tnrs_match_names() from the `rotl` package to build a phylogenetic tree
# with taxonomic data from Open Tree Taxonomy

######## FECUNDITY ANALYSES TREE ###########

# Tree for fecundity analyses was created and saved on August 30th, 2025
load("data/fecund_tree.Rdata")

# Plotting fecundity tree
plot(tree, cex = 0.5, label.offset = 0.25, no.margin = TRUE)

### Checking that tree is complete
data_fecund <- read.csv("data/fecund_data.csv") # load in dataset to compare

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

# Save phylogenetic correlation matrix 
save(phylo_cor, file = "data/phylo_cor_fecund.Rdata")

######## LONGEVITY ANALYSES TREE ###########
# Creating tree
taxa_long <- tnrs_match_names(c("Aacanthocnema dobsoni",
                                "Adalia bipunctata",
                                "Aedes aegypti",
                                "Agrilus coxalis auroguttatus",
                                "Agriphila aeneociliella",
                                "Allonemobius socius",
                                "Amblyseius largoensis",
                                "Amblyseius zaheri",
                                "Anaphes nitens",
                                "Anastrepha ludens",
                                "Anegleis cardoni",
                                "Anthocoris minki",
                                "Aphelinus asychis",
                                "Arctia plantaginis",
                                "Atrophaneura alcinous",
                                "Bactrocera dorsalis",
                                "Cadra cautella",
                                "Callosobruchus analis",
                                "Callosobruchus chinensis",
                                "Callosobruchus maculatus",
                                "Callosobruchus phaseoli",
                                "Callosobruchus rhodesianus",
                                "Callosobruchus subinnotatus",
                                "Chrysochus cobaltinus",
                                "Cimex lectularius",
                                "Cnaphalocrocis medinalis",
                                "Coelopa frigida",
                                "Colaphellus bowringi",
                                "Copitarsia decolora",
                                "Cotinis nitida",
                                "Danaus plexippus",
                                "Diabrotica barberi",
                                "Diaphorina citri",
                                "Drosophila arizonae",
                                "Drosophila mauritiana",
                                "Drosophila melanogaster",
                                "Drosophila mojavensis",
                                "Drosophila sechellia",
                                "Drosophila simulans",
                                "Echinothrips americanus",
                                "Ellychnia corrusca",
                                "Engytatus varians",
                                "Eublaberus posticus",
                                "Galerucella birmanica",
                                "Gnatocerus cornutus",
                                "Grapholita molesta",
                                "Gryllodes sigillatus",
                                "Gryllus bimaculatus",
                                "Gryllus lineaticeps",
                                "Gryllus vocalis",
                                "Helicoverpa armigera",
                                "Histiostoma feroniarum",
                                "Kampimodromus aberrans",
                                "Lasioderma serricorne",
                                "Lygaeus equestris",
                                "Lygocoris pabulinus",
                                "Microplitis rufiventris",
                                "Mythimna unipuncta",
                                "Nasonia vitripennis",
                                "Neoseiulus californicus",
                                "Neoseiulus cucumeris",
                                "Ophraella communa",
                                "Ostrinia nubilalis",
                                "Parasitus fimetorum",
                                "Pardosa astrigera",
                                "Photinus ignitus",
                                "Pieris napi",
                                "Plagiodera versicolora",
                                "Plodia interpunctella",
                                "Plutella xylostella",
                                "Podisus nigrispinus",
                                "Protohermes grandis",
                                "Pyrocoelia pectoralis",
                                "Riptortus clavatus",
                                "Saltella sphondylii",
                                "Sancassania berlesei",
                                "Sangalopsis microleuca",
                                "Sepsis cynipsea",
                                "Sitophilus oryzae",
                                "Spodoptera litura",
                                "Tenebrio molitor",
                                "Tetranychus urticae",
                                "Thaumatotibia leucotreta",
                                "Tribolium castaneum",
                                "Trichogramma evanescens",
                                "Trichoplusia ni",
                                "Xylocoris flavipes",
                                "Yponomeuta cagnagellus",
                                "Yponomeuta padellus"))

long_tree <- tol_induced_subtree(ott_ids =
                                 taxa_long[,"ott_id"],
                                 label_format = "name")

long_tree$tip.label <- gsub("_", " ", long_tree$tip.label)

plot(long_tree, cex = 0.5, label.offset = 0.25, no.margin = TRUE)

save(long_tree, file = "data/long_tree.Rdata")

######## Checking with longevity data to see that all taxa are there #########
data_long <- read.csv("data/long_data.csv") # load in longevity data to compare

# listed in tree but not in data
setdiff(as.character(long_tree$tip.label),
        as.character(data_long$species))

# listed in data but not in tree
setdiff(as.character(data_long$species),
        as.character(long_tree$tip.label))


#################### Extracting correlation matrix ###########################
# Compute branch lengths of tree
phylo_branch_long <- compute.brlen(long_tree, method = "Grafen", power = 1)

is.ultrametric(phylo_branch_long)

phylo_cor_long <- vcv(phylo_branch_long, cor = T)

# visual exploration of the phylogenetic correlation matrix
ggcorrplot::ggcorrplot(phylo_cor_long, sig.level = 0.05, lab_size = 1,
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

# Save phylogenetic correlation matrix 
save(phylo_cor, file = "data/phylo_cor_long.Rdata")








