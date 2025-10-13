library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(DHARMa)
library(metafor)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 16))

# For installing the orchaRd package
#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

# SCRIPT FOR ANALYZING EFFECT OF POLYANDRY ON LONGEVITY
# Load data in
data_long <- read.csv("data/long_data.csv") |>
             filter(nzchar(source))

############################ CLEANING TREATMENT DATA BY COMBINING INTO FOUR LEVELS #############################
unique(data_long$treatment)
length(unique(data_long$study)) # 109 studies

# once-few treatments
data_long$treatment[data_long$treatment %in% paste0("once_", c("thrice", "three", "four"))] <- "once_few"

## once-many treatments
many_vec <- c("five", "six", "seven", "eight", "nine",
              "ten", "twelve", "fifteen", "twenty", "twentyfive")
data_long$treatment[data_long$treatment %in% paste0("once_", many_vec)] <- "once_many"

## few-many treatments
fewmany_vec <- c("twice_many", "three_five", "twice_ten", "five_ten", "five_fifteen")
data_long$treatment[data_long$treatment %in% fewmany_vec] <- "few_many"

# Inspect number of factor levels and N for each level for treatment
data_long$treatment <- as.factor(data_long$treatment)
aggregate(data_long$treatment, by = list(data_long$treatment), FUN = length) # get count of N per treatment

# Reorder factor levels for treatment
data_long$treatment <- factor(data_long$treatment, levels = c("few_many", "once_many", "once_few", "once_twice"))

############################ CALCULATE EFFECT SIZES ##################################

data_long <- escalc(measure = "ROM", data = data_long,
                      m1i = exp_long, 
                      m2i = con_long, 
                      sd1i = exp_long_SD, 
                      sd2i = con_long_SD, 
                      n1i = exp_N, n2i = adjusted_con_N, 
                      slab = source, vtype = "AVHO")

length(data_long$yi) # 242 effect sizes

# CREATING A VARIANCE-COVARIANCE MATRIX 
# Code taken from Mentesana et al. 2025: 10.5281/zenodo.14930059
# Creating a var-covar matrix assuming a 0.5 correlation between effect sizes 
# from the same study. covariance = (0.5 * sqrt(vi.1) * sqrt(vi.2))
# Creates a matrix (called 'VCV_ESVar') with the dimensions =  
# n(effect_sizes) x n(effect_sizes)

VCV_ESVar <- matrix(0, nrow = nrow(data_long), 
                    ncol = nrow(data_long))

# Names rows and columns for each obsID
rownames(VCV_ESVar) <- data_long[, "X"]
colnames(VCV_ESVar) <- data_long[, "X"]

# Finds effect sizes that come from the same study
shared_coord <- which(data_long[, "study"] %in% 
                        data_long[duplicated(data_long[, "study"]), 
                                    "study"] == TRUE)

combinations <- do.call("rbind", tapply(shared_coord, 
                                        data_long[shared_coord, "study"], 
                                        function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations)[1]) {
  p1 <- combinations[i, 1]
  p2 <- combinations[i, 2]
  p1_p2_cov <- 0.5 * sqrt(data_long[p1, "vi"]) * 
    sqrt(data_long[p2, "vi"])
  VCV_ESVar[p1, p2] <- p1_p2_cov
  VCV_ESVar[p2, p1] <- p1_p2_cov
} 

diag(VCV_ESVar) <- data_long[, "vi"]

########################### OVERALL MODEL ############################################
data_long$experiment <- as.factor(data_long$experiment)
data_long$study <- as.factor(data_long$study)
data_long$species <- as.factor(data_long$species)

# Loading in phylogenetic data
load("data/long_tree.Rdata")
load("data/phylo_cor_long.Rdata")

# Creating a duplicate species variable for the phylogenetic analysis
data_long$species_phylo <- data_long$species

overall_model <- rma.mv(yi, VCV_ESVar, data = data_long, 
                        random = list( ~ 1|study/experiment,
                                       ~ 1|species,
                                       ~ 1|species_phylo),
                        R = list(species_phylo = phylo_cor_long),
                        method = "REML")

summary(overall_model)
forest(overall_model)

# Calculating heterogeneity
i2_ml(overall_model)

#### TREATMENT MODEL (sig)
treatment_model <- rma.mv(yi, VCV_ESVar, data = data_long, 
                          random = list( ~ 1|study/experiment,
                                         ~ 1|species,
                                         ~ 1|species_phylo),
                          R = list(species_phylo = phylo_cor_long),
                          method = "REML",
                          mods = ~ 0 + treatment)

summary(treatment_model)

(treatment_fig <- orchard_plot(treatment_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "treatment", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("#004e89","#2A8EDB","#ABDBFF","#E9F5FF")) +
             scale_color_manual(values = c("grey4", "grey4", "grey4", "grey4")))

ggsave(treatment_fig, filename = "fig_long_treatment.png", width = 6, height = 4)

#### HARASSMENT MODEL (not sig)
harass_model <- rma.mv(yi, VCV_ESVar, data = data_long, 
                       random = list( ~ 1|study/experiment,
                                      ~ 1|species,
                                      ~ 1|species_phylo),
                       R = list(species_phylo = phylo_cor_long),
                       method = "REML",
                       mods = ~ 1 + harass.)

summary(harass_model)

(harass_fig <- orchard_plot(harass_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "harass.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")))

ggsave(harass_fig, filename = "fig_long_harass.png", width = 6, height = 4)

#### NUPTIAL GIFT MODEL (marg sig)
gift_model <- rma.mv(yi, VCV_ESVar, data = data_long, 
                     random = list( ~ 1|study/experiment,
                                    ~ 1|species,
                                    ~ 1|species_phylo),
                     R = list(species_phylo = phylo_cor_long),
                     method = "REML",
                     mods = ~ 1 + nup_gift.)

summary(gift_model)

(gift_fig <- orchard_plot(gift_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "nup_gift.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
  scale_fill_manual(values = c("grey100", "grey25")) +
  scale_color_manual(values = c("grey4", "grey4")))

ggsave(gift_fig, filename = "fig_long_gift.png", width = 6, height = 4)


#### SELECTION BIAS MODEL (sig)
bias_model <- rma.mv(yi, VCV_ESVar, data = data_long, 
                     random = list( ~ 1|study/experiment,
                                    ~ 1|species,
                                  ~ 1|species_phylo),
                     R = list(species_phylo = phylo_cor_long),
                     method = "REML",
                     mods = ~ 1 + bias.)

summary(bias_model)

(bias_fig <- orchard_plot(bias_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "bias.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")))

ggsave(bias_fig, filename = "fig_long_bias.png", width = 6, height = 4)

################################ PUBLICATION BIAS #########################################
### Funnel plot visual inspection
funnel(overall_model, yaxis="seinv", xlab="Effect size (log odds ratio)", 
       back="white", col=rgb(0,153,76, max=255, alpha=125), digits = 1) # nicer funnel plot

data_long$precision <- sqrt(1/data_long$vi)

# # meta-regressions testing for publication bias
egger_model_long <- rma.mv(yi, VCV_ESVar, data = data_long, 
                             random = list( ~ 1|study/experiment,
                                            ~ 1|species,
                                            ~ 1|species_phylo),
                             R = list(species_phylo = phylo_cor_long),
                             method = "REML",
                             mods = ~ precision)
summary(egger_model_long)


#### Small study bias
# Compute effective sample size
data_long$treatment <- factor(data_long$treatment, levels = c("once_twice","once_few", "once_many", "few_many"))

data_long$inv_ESS <- (data_long$exp_N + data_long$con_N) / (data_long$exp_N *
                                                                    data_long$con_N)
data_long$sqrt_inv_ESS <- sqrt(data_long$inv_ESS)

# Small-study effects (SME) model
small_study_long_model <- rma.mv(yi, vi, data = data_long,
                                   random = list( ~ 1|study/experiment,
                                                  ~ 1|species,
                                                  ~ 1|species_phylo),
                                   R = list(species_phylo = phylo_cor_long),
                                   mods = ~ sqrt_inv_ESS*treatment)

summary(small_study_long_model)

anova(small_study_long_model, btt = 6:9)

# Plot for small-study effects model
orchaRd::bubble_plot(small_study_long_model, 
                     mod = "sqrt_inv_ESS", group = "study",
                     xlab = "square root inverse effective sample size", 
                     legend.pos = "bottom.right", 
                     by = "treatment") + geom_blank(data = data_long, aes(color = treatment)) + 
                     scale_fill_manual(values = c("#E9F5FF","#ABDBFF","#2A8EDB","#004e89"))

##### Time-lag bias or decline effects
# Extract year from source and center it
data_long$year <- as.integer(unlist(str_extract_all(data_long$source, "\\d+")))
data_long$year.c <- data_long$year - mean(data_long$year)

# Decline effects (DE) model
time_long_model <- rma.mv(yi, vi, data = data_long,
                            random = list( ~ 1|study/experiment,
                                           ~ 1|species,
                                           ~ 1|species_phylo),
                            R = list(species_phylo = phylo_cor_long),
                            mods = ~ 1 + year.c*treatment) # significant decline effects

summary(time_long_model)
anova(time_long_model, btt = 6:8)

orchaRd::bubble_plot(time_long_model, mod = "year.c", 
                     group = "study",
                     xlab = "Mean-centered year",
                     legend.pos = "bottom.left", 
                     by = "treatment") + geom_blank(data = data_long, aes(color = treatment)) + 
                     scale_fill_manual(values = c("#E9F5FF","#ABDBFF","#2A8EDB","#004e89"))

