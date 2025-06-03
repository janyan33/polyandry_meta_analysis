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
data_long <- read.csv("long_data.csv") |>
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
                      n1i = exp_N, n2i = con_N, 
                      slab = source, vtype = "AVHO")

length(data_long$yi) # 242 effect sizes

########################### OVERALL MODEL ############################################
data_long$experiment <- as.factor(data_long$experiment)
data_long$study <- as.factor(data_long$study)

overall_model <- rma.mv(yi, vi, data = data_long, random = ~ 1|study/experiment,
                        method = "REML")

summary(overall_model)
forest(overall_model)

# Calculating heterogeneity
i2_ml(overall_model)

#### TREATMENT MODEL (sig)
treatment_model <- rma.mv(yi, vi, data = data_long, random = ~ 1|study/experiment,
                          method = "REML",
                          mods = ~ 0 + treatment)

summary(treatment_model)

treatment_fig <- orchard_plot(treatment_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "treatment", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("#004e89","#2A8EDB","#ABDBFF","#E9F5FF")) +
             scale_color_manual(values = c("grey4", "grey4", "grey4", "grey4"))

#### HARASSMENT MODEL (no sig)
harass_model <- rma.mv(yi, vi, data = data_long, random = ~ 1|study/experiment,
                       method = "REML",
                       mods = ~ 1 + harass.)

summary(harass_model)

harass_fig <- orchard_plot(harass_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "harass.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4"))

#### NUPTIAL GIFT MODEL (sig)
gift_model <- rma.mv(yi, vi, data = data_long, random = ~ 1|study/experiment,
                     method = "REML",
                     mods = ~ 0 + nup_gift.)

summary(gift_model)

gift_fig <- orchard_plot(gift_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "nup_gift.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
  scale_fill_manual(values = c("grey100", "grey25")) +
  scale_color_manual(values = c("grey4", "grey4"))


#### SELECTION BIAS MODEL (sig)
bias_model <- rma.mv(yi, vi, data = data_long, random = ~ 1|study/experiment,
                     method = "REML",
                     mods = ~ 0 + bias.)

summary(bias_model)

bias_fig <- orchard_plot(bias_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "bias.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4"))

## TAXONOMIC GROUP MODEL
aggregate(data_long$order, by = list(data_long$order), FUN = length) 

# Filter for only categories with N > 5
data_long_order <- data_long %>% 
                   filter(order != "Blattodea" &
                          order != "Megaloptera" &
                          order != "Thysanoptera")


order_model <- rma.mv(yi, vi, data = data_long_order, random = ~ 1|study/experiment,
                      method = "REML",
                      mods = ~ 1 + order) 

summary(order_model)

orchard_plot(order_model, xlab = "Effect size (log response ratio)", 
             group = "study",
             mod = "order", twig.size = NA, branch.size = 2, trunk.size = 0.75, 
             angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top")





