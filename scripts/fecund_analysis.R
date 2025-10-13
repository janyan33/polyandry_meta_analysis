library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(DHARMa)
library(metafor)
library(clubSandwich)
library(orchaRd)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 16))

par(mar = c(5.1, 4.1, 4.1, 2.1)) 

########################### CLEANING TREATMENT DATA BY COMBINING INTO FOUR LEVELS ##############################
data_fecund <- read.csv("data/fecund_data.csv") |>
               filter(nzchar(source))

unique(data_fecund$treatment)
length(unique(data_fecund$study)) # 157 studies

# once-few treatments
data_fecund$treatment[data_fecund$treatment %in% 
                     paste0("once_", c("thrice", "three", "four"))] <- "once_few"

## once-many treatments
many_vec <- c("five", "six", "seven", "eight", "nine",
              "ten", "twelve", "fifteen", "twenty", "twentyfive")
data_fecund$treatment[data_fecund$treatment %in% paste0("once_", many_vec)] <- "once_many"

## few-many treatments
fewmany_vec <- c("twice_many", "three_five", "twice_ten", "five_ten", "five_fifteen")
data_fecund$treatment[data_fecund$treatment %in% fewmany_vec] <- "few_many"

# Inspect number of factor levels and N for each level for treatment
data_fecund$treatment <- as.factor(data_fecund$treatment)
aggregate(data_fecund$treatment, by = list(data_fecund$treatment), FUN = length) # get count of N per treatment

# Reorder factor levels for treatment
data_fecund$treatment <- factor(data_fecund$treatment, levels = c("few_many", "once_many", "once_few", "once_twice"))

######################################### CALCULATE EFFECT SIZES ###############################################
data_fecund$exp_fert <- as.numeric(data_fecund$exp_fert)
data_fecund$con_fert <- as.numeric(data_fecund$con_fert)

data_fecund <- escalc(measure = "ROM", data = data_fecund,
                      m1i = ifelse(is.na(exp_fert), exp_fecund, exp_fert), 
                      m2i = ifelse(is.na(con_fert), con_fecund, con_fert), 
                      sd1i = ifelse(is.na(exp_fert), exp_fecund_SD, exp_fert_SD), 
                      sd2i = ifelse(is.na(con_fert_SD), con_fecund_SD, con_fert_SD), 
                      n1i = exp_N, n2i = adjusted_con_N, 
                      slab = source, 
                      vtype = "AVHO")

length(data_fecund$yi) # 325 effect sizes

# CREATING A VARIANCE-COVARIANCE MATRIX 
# Code taken from Mentesana et al. 2025: 10.5281/zenodo.14930059
# Creating a var-covar matrix assuming a 0.5 correlation between effect sizes 
# from the same study. covariance = (0.5 * sqrt(vi.1) * sqrt(vi.2))
# Creates a matrix (called 'VCV_ESVar') with the dimensions =  
# n(effect_sizes) x n(effect_sizes)

VCV_ESVar <- matrix(0, nrow = nrow(data_fecund), 
                    ncol = nrow(data_fecund))

# Names rows and columns for each obsID
rownames(VCV_ESVar) <- data_fecund[, "X"]
colnames(VCV_ESVar) <- data_fecund[, "X"]

# Finds effect sizes that come from the same study
shared_coord <- which(data_fecund[, "study"] %in% 
                        data_fecund[duplicated(data_fecund[, "study"]), 
                                    "study"] == TRUE)

combinations <- do.call("rbind", tapply(shared_coord, 
                                        data_fecund[shared_coord, "study"], 
                                        function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations)[1]) {
  p1 <- combinations[i, 1]
  p2 <- combinations[i, 2]
  p1_p2_cov <- 0.5 * sqrt(data_fecund[p1, "vi"]) * 
    sqrt(data_fecund[p2, "vi"])
  VCV_ESVar[p1, p2] <- p1_p2_cov
  VCV_ESVar[p2, p1] <- p1_p2_cov
} 

diag(VCV_ESVar) <- data_fecund[, "vi"]

############################################# OVERALL MODEL ####################################################
data_fecund$experiment <- as.factor(data_fecund$experiment)
data_fecund$study <- as.factor(data_fecund$study)
data_fecund$species <- as.factor(data_fecund$species)

# Loading in phylogenetic data
load("data/fecund_tree.Rdata")
load("data/phylo_cor_fecund.Rdata")

# Creating a duplicate species variable for the phylogenetic analysis
# I don't entirely understand why this is necessary
data_fecund$species_phylo <- data_fecund$species

overall_model <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                        random = list( ~ 1|study/experiment,
                                       ~ 1|species,
                                       ~ 1|species_phylo),
                        R = list(species_phylo = phylo_cor),
                        method = "REML")

summary(overall_model)
forest(overall_model) # checked the extreme values, not errors

# Calculating heterogeneity
i2_ml(overall_model)

######################################### TREATMENT MODERATOR MODEL #############################################
treatment_model <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                          random = list( ~ 1|study/experiment,
                                         ~ 1|species,
                                         ~ 1|species_phylo),
                          R = list(species_phylo = phylo_cor),
                          method = "REML",
                          mods = ~ 0 + treatment)

summary(treatment_model)
r2_ml(treatment_model) # Get marginal r2

(treatment_fig <- orchard_plot(treatment_model, xlab = "Effect size (log response ratio)", 
             group = "source",  
             mod = "treatment", twig.size = 0.5, branch.size = 2, trunk.size = 0.6, angle = 45, flip = TRUE, 
             alpha = 0.4, g = T) + theme(legend.position = "top") + ylim(-2, 2) +
             scale_fill_manual(values = c("#004e89","#2A8EDB","#ABDBFF", "#E9F5FF")) +
             scale_color_manual(values = c("grey20", "grey20", "grey20", "grey20")))

ggsave(treatment_fig, filename = "fig_fecund_treatment.png", width = 6, height = 4)

#################################### CONTINOUS HOUSING MODERATOR MODEL ###########################################
harass_model <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                       random = list( ~ 1|study/experiment,
                                      ~ 1|species,
                                      ~ 1|species_phylo),
                       R = list(species_phylo = phylo_cor),
                       mods = ~ 1 + harass.)

summary(harass_model)
r2_ml(harass_model) # Get marginal r2

(harass_fig <- orchard_plot(harass_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "harass.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65))

ggsave(harass_fig, filename = "fig_fecund_harass.png", width = 6, height = 4)

#################################### PARTIAL VS. LIFETIME FITNESS MODEL ###########################################
lifetime_model <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                         random = list( ~ 1|study/experiment,
                                        ~ 1|species,
                                        ~ 1|species_phylo),
                         R = list(species_phylo = phylo_cor),
                         mods = ~ 1 + lifetime.) # change the intercept from 1 to 0 to get model results relative to zero

summary(lifetime_model)
r2_ml(lifetime_model) # Get marginal r2

(lifetime_plot <- orchard_plot(lifetime_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "lifetime.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65))

ggsave(lifetime_plot, filename = "fig_fecund_lifetime.png", width = 6, height = 4)

############################################## NUPTIAL GIFT MODEL #################################################
gift_model <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                     random = list( ~ 1|study/experiment,
                                    ~ 1|species,
                                    ~ 1|species_phylo),
                     R = list(species_phylo = phylo_cor),
                     mods = ~ 1 + nup_gift.) # change the intercept from 1 to 0 to get model results relative to zero

summary(gift_model)
r2_ml(gift_model) # Get marginal r2

(gift_fig <- orchard_plot(gift_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "nup_gift.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") +  
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65))

ggsave(gift_fig, filename = "fig_fecund_gift.png", width = 6, height = 4)
            
############################################## SELECTION BIAS MODEL ##############################################
bias_model <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                     random = list( ~ 1|study/experiment,
                                    ~ 1|species,
                                    ~ 1|species_phylo),
                     R = list(species_phylo = phylo_cor),
                     mods = ~ 1 + bias.) # change the intercept from 1 to 0 to get model results relative to zero

summary(bias_model)
r2_ml(bias_model) # Get marginal r2

(bias_fig <- orchard_plot(bias_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "bias.", twig.size = 0.5, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65))

ggsave(bias_fig, filename = "fig_fecund_bias.png", width = 6, height = 4)

################################ PUBLICATION BIAS ###############################################
### Funnel plot visual inspection
funnel(overall_model, yaxis="seinv", xlab="Effect size (log odds ratio)", 
       back="white", col=rgb(0,153,76, max=255, alpha=125), digits = 1)

### Meta-regressions testing for publication bias
data_fecund$precision <- sqrt(1/data_fecund$vi)

egger_model_fecund <- rma.mv(yi, VCV_ESVar, data = data_fecund, 
                                  random = list( ~ 1|study/experiment,
                                                 ~ 1|species,
                                                 ~ 1|species_phylo),
                                                   R = list(species_phylo = phylo_cor),
                                                   method = "REML",
                                                   mods = ~ precision)

summary(egger_model_fecund)

#### Small study bias
# Compute effective sample size
data_fecund$treatment <- factor(data_fecund$treatment, levels = c("once_twice","once_few", "once_many", "few_many"))
data_fecund$inv_ESS <- (data_fecund$exp_N + data_fecund$con_N) / (data_fecund$exp_N *
                                                                    data_fecund$con_N)
data_fecund$sqrt_inv_ESS <- sqrt(data_fecund$inv_ESS)

# Small-study effects (SME) model
small_study_fecund_model <- rma.mv(yi, vi, data = data_fecund,
                            random = list( ~ 1|study/experiment,
                                           ~ 1|species,
                                           ~ 1|species_phylo),
                           R = list(species_phylo = phylo_cor),
                           mods = ~ sqrt_inv_ESS*treatment)

summary(small_study_fecund_model)

anova(small_study_fecund_model, btt = 6:8) # testing interaction (not significant)

# Plot for small-study effects model
orchaRd::bubble_plot(small_study_fecund_model, 
                     mod = "sqrt_inv_ESS", group = "study",
                     xlab = "squareroot inverse effective sample size", 
                     legend.pos = "bottom.right", 
                     by = "treatment") + geom_blank(data = data_fecund, aes(color = treatment)) + 
                     scale_fill_manual(values = c("#E9F5FF","#ABDBFF","#2A8EDB","#004e89"))

##### Time-lag bias or decline effects
# Extract year from source and center it
data_fecund$year <- as.integer(unlist(str_extract_all(data_fecund$source, "\\d+")))
data_fecund$year.c <- data_fecund$year - mean(data_fecund$year)

# Decline effects (DE) model
time_fecund_model <- rma.mv(yi, vi, data = data_fecund,
                            random = list( ~ 1|study/experiment,
                                           ~ 1|species,
                                           ~ 1|species_phylo),
                            R = list(species_phylo = phylo_cor),
                            mods = ~  year.c*treatment) # significant decline effects

summary(time_fecund_model)
anova(time_fecund_model, btt = 6:8)

orchaRd::bubble_plot(time_fecund_model, mod = "year.c", 
                                        group = "study",
                                        xlab = "Mean-centered year",
                                        legend.pos = "bottom.left", 
                                        by = "treatment") + geom_blank(data = data_fecund, aes(color = treatment)) + 
                                        scale_fill_manual(values = c("#E9F5FF","#ABDBFF","#2A8EDB","#004e89"))
