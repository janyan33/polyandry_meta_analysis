library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(DHARMa)
library(metafor)

# For installing the orchaRd package
install.packages("pacman")
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 16))

########################### CLEANING TREATMENT DATA BY COMBINING INTO FOUR LEVELS ##############################
data_fecund <- read.csv("fecund_data.csv") |>
               filter(nzchar(source))

unique(data_fecund$treatment)
length(unique(data_fecund$study)) # 157 studies

# once-few treatments
data_fecund$treatment[data_fecund$treatment %in% paste0("once_", c("thrice", "three", "four"))] <- "once_few"

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
                      n1i = exp_N, n2i = con_N, 
                      slab = source, 
                      vtype = "AVHO")

length(data_fecund$yi) # 325 effect sizes

############################################# OVERALL MODEL ####################################################
data_fecund$experiment <- as.factor(data_fecund$experiment)
data_fecund$study <- as.factor(data_fecund$study)

overall_model <- rma.mv(yi, vi, data = data_fecund, random = ~ 1|study/experiment,
                       method = "REML")

summary(overall_model)
forest(overall_model) # checked the extreme values, not errors

# Calculating heterogeneity
i2_ml(overall_model)

######################################### TREATMENT MODERATOR MODEL #############################################
treatment_model <- rma.mv(yi, vi, data = data_fecund, random = ~ 1|study/experiment,
                        method = "REML",
                        mods = ~ 0 + treatment)

summary(treatment_model)

(treatment_fig <- orchard_plot(treatment_model, xlab = "Effect size (log response ratio)", 
             group = "source",  
             mod = "treatment", twig.size = NA, branch.size = 2, trunk.size = 0.6, angle = 45, flip = TRUE, 
             alpha = 0.4, g = T) + theme(legend.position = "top") + ylim(-1.5, 1.65) +
             scale_fill_manual(values = c("#004e89","#2A8EDB","#ABDBFF", "#E9F5FF")) +
             scale_color_manual(values = c("grey20", "grey20", "grey20", "grey20")))

ggsave(treatment_fig, filename = "fig_fecund_treatment.png", width = 6, height = 4)

#################################### CONTINOUS HOUSING MODERATOR MODEL ###########################################
harass_model <- rma.mv(yi, vi, data = data_fecund, random = ~ 1|study/experiment,
                          method = "REML",
                          mods = ~ 1 + harass.) # change the intercept from 1 to 0 to get model results relative to zero

summary(harass_model)

harass_fig <- orchard_plot(harass_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "harass.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65)

ggsave(harass_fig, filename = "fig_fecund_harass.png", width = 6, height = 4)

#################################### PARTIAL VS. LIFETIME FITNESS MODEL ###########################################
lifetime_model <- rma.mv(yi, vi, data = data_fecund, random = ~ 1|study/experiment,
                       method = "REML",
                       mods = ~ 1 + lifetime.) # change the intercept from 1 to 0 to get model results relative to zero

summary(lifetime_model)

(lifetime_plot <- orchard_plot(lifetime_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "lifetime.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65))

ggsave(lifetime_plot, filename = "fig_fecund_lifetime.png", width = 6, height = 4)

############################################## NUPTIAL GIFT MODEL #################################################
gift_model <- rma.mv(yi, vi, data = data_fecund, random = ~ 1|study/experiment,
                         method = "REML",
                         mods = ~ 1 + nup_gift.) # change the intercept from 1 to 0 to get model results relative to zero

summary(gift_model)

(gift_fig <- orchard_plot(gift_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "nup_gift.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") +  
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65))

ggsave(gift_fig, filename = "fig_fecund_gift.png", width = 6, height = 4)
            
############################################## SELECTION BIAS MODEL ##############################################
bias_model <- rma.mv(yi, vi, data = data_fecund, random = ~ 1|study/experiment,
                         method = "REML",
                         mods = ~ 1 + bias.) # change the intercept from 1 to 0 to get model results relative to zero

summary(bias_model)

bias_fig <- orchard_plot(bias_model, xlab = "Effect size (log response ratio)", 
             group = "study",  
             mod = "bias.", twig.size = NA, branch.size = 2, trunk.size = 0.75, angle = 45, flip = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + 
             scale_fill_manual(values = c("grey100", "grey25")) +
             scale_color_manual(values = c("grey4", "grey4")) + ylim(-1.65, 1.65)

ggsave(bias_fig, filename = "fig_fecund_bias.png", width = 6, height = 4)

################################################# ORDER MODEL #####################################################
# Filter for only categories with N > 5
aggregate(data_fecund$order, by = list(data_fecund$order), FUN = length) 

data_fecund_order <- data_fecund %>% 
                     filter(order != "Blattodea" & # N = 3
                     order != "Megaloptera" & # N = 2
                     order != "Thysanoptera" & # N = 1
                     order != "Decapoda" & # N = 1
                     order != "Mantodea" & # N = 2
                     order != "Neuroptera" & # N = 3
                     order != "Siphonaptera") # N = 1

order_fecund_model <- rma.mv(yi, vi, data = data_fecund_order, random = ~ 1|study/experiment,
                      method = "REML",
                      mods = ~ 1 + order) # change the intercept from 1 to 0 to get model results relative to zero

summary(order_fecund_model)

order_fig <- orchard_plot(order_fecund_model, xlab = "Effect size (log response ratio)", 
             group = "study",
             mod = "order", twig.size = NA, branch.size = 2, trunk.size = 0.75, 
             angle = 45, flip = TRUE, cb = TRUE, 
             alpha = 0.35, g = TRUE) + theme(legend.position = "top") + ylim(-1.65, 1.65)

ggsave(order_fig, filename = "fig_fecund_order.png", width = 6, height = 4)

