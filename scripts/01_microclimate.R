####Libraries####
library(lme4)
library(MuMIn)
library(patchwork)
library(tidyverse)

source("scripts/_internals.R")

####Import data####
env <- read.table("data/microclimate.txt", 
                  header = T) %>%
  # remove some columns
  dplyr::select(!c(Wpt, Alt, TimeH, nn1, nn2, TempA, CovU, CovC, RichU, RichC)) %>%
  left_join(.,
            fg_cover,
            by = c("Site" = "Site", "Pair" = "Pair")) %>% 
  # data in long format
  pivot_longer(.,
               cols = TempU:richness_grass_U,
               names_to = "variable") %>%
  # create variable for microsite (and drop from variable)
  mutate(microsite = as.factor(str_extract(variable, ".{1}$")),
         variable = str_replace(variable, ".{1}$", ""))

####Mixed effect models####

# get names of all the variables
vars <- distinct(env, variable) %>%
  pull()

# create empty data frame
model_results = data.frame(variable = vars, 
                           intercept = NA, 
                           intercpet_stderror = NA, 
                           intercpet_tval = NA,
                           micrositeU = NA, 
                           micrositeU_stderror = NA,
                           micrositeU_tval = NA,
                           R_marginal = NA,
                           R_conditional = NA)

for (i in 1:length(vars)) {
  
  dat <- env %>%
    filter(variable == vars[i])
  
  # model
  lmm <- lmer(formula = value ~ microsite + (1 | Site), 
              dat)
  summ = summary(lmm)$coefficients
  
  # extract relevant summary stats
  model_results[i, 2] <- summ[1,1]
  model_results[i, 3] <- summ[1,2]
  model_results[i, 4] <- summ[1,3]
  model_results[i, 5] <- summ[2,1]
  model_results[i, 6] <- summ[2,2]
  model_results[i, 7] <- summ[2,3]
  model_results[i, 8] <- r.squaredGLMM(lmm)[1]
  model_results[i, 9] <- r.squaredGLMM(lmm)[2]
  
}

write.csv(model_results,
          "outputs/mixed_effect.csv")

# plot pairwise data

env_plot <- 
  env %>%
  mutate(variable = case_when(variable == 'Temp' ~ "Temperature (C)",
                              variable == 'STemp' ~ "Soil temperature (C)",
                              variable == 'SMoist' ~ "Soil moisture",
                              variable == 'cover_forb_' ~ "Forb cover (%)",
                              variable == 'cover_grass_' ~ "Grass cover (%)",
                              variable == 'richness_forb_' ~ "Forb species richness",
                              variable == 'richness_grass_' ~ "Grass species richness"),
         microsite = case_when(as.character(microsite) == "U" ~ "Under",
                               .default = "Away"))

plots <- vector('list', 3)
grp_vars <- vector('list', 3)
grp_vars[[1]] <- c("Temperature (C)", "Soil temperature (C)", "Soil moisture")
grp_vars[[3]] <- c("Forb cover (%)", "Grass cover (%)")
grp_vars[[2]] <- c("Forb species richness", "Grass species richness")

for (i in 1:length(plots)) {
  
  dat <- env_plot %>% 
          filter(variable %in% grp_vars[[i]])
  
  plots[[i]] <- ggplot(dat,
       aes(x = variable,
           y = value)) +
  geom_point(data = dat %>%
               filter(microsite == "Under"),
             aes(x = 1.2,
                 y = value),
             alpha = 0.3,
             fill = 'forestgreen',
             colour = "white",
             shape = 21,
             position = position_jitternormal(sd_x = 0.05, sd_y = 0)) +
  geom_point(data = dat %>%
               filter(microsite == "Away"),
             aes(x = 0.8,
                 y = value),
             alpha = 0.3,
             fill = 'goldenrod1',
             colour = "white",
             shape = 21,
             position = position_jitternormal(sd_x = 0.05, sd_y = 0)) +
  geom_boxplot(aes(colour = microsite),
               outliers = FALSE,
               fill = NA) +
  facet_wrap(vars(variable),
             scales = "free") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = c('goldenrod1','forestgreen'),
                      name = "Microsite") +
  labs(y = "Value",
       x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.caption = element_markdown(),
    legend.position = 'bottom'
  )
  
}

plots[[1]] + 
  labs(tag = "A") + 
plots[[2]] + 
  labs(tag = "B") + 
plots[[3]] + 
  labs(tag = "C") +
  plot_layout(ncol = 1,
  guides = 'collect') +
  plot_annotation(theme = theme(
    legend.position = 'bottom'))


ggsave("figures/microclimate_boxplot.png",
       width = 7,
       height = 12)


