library(tidyverse)

options(digits = 2)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 16))


### SCRIPT FOR GENERATING A CLEAN LONGEVITY DATA SET
long_data <- read.csv("data/polyandry_meta_full_dataset.csv") %>% 
            filter(is.na(con_long) == F & is.na(exp_long) ==F) %>% # remove rows with no longevity data
            select(-con_fecund, -con_fecund_SD, -exp_fecund, -exp_fecund_SD, 
                   -con_fert, -con_fert_SD, -exp_fert, -exp_fert_SD)

write.csv(long_data, "data/long_data.csv")

### SCRIPT FOR GENERATING A CLEAN FECUNDITY DATA SET
fecund_data <- read.csv("data/polyandry_meta_full_dataset.csv") %>%
               filter(!(is.na(con_fecund) == T & is.na(con_fert) == T)) %>% # remove rows with no fecundity data
               select(-con_long, -con_long_SD, -exp_long, -exp_long_SD)

write.csv(fecund_data, "data/fecund_data.csv")


### GLOBAL EFFECT SIZE PLOT
global_results <- read.csv("data/global_model_results.csv")

global_results$fitness_metric <- factor(global_results$fitness_metric, levels = c("longevity", "fecundity"))

ggplot(data = global_results, aes(x = estimate, y = fitness_metric)) + My_Theme + 
  xlab("Effect size (Log response ratio)") + ylab("") + theme(legend.position = "none") + 
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(data = global_results, aes(x = estimate, xmin = ci.lb, xmax = ci.ub), color = "grey20", 
                  lwd = 2, fatten = 5, shape = 23)
