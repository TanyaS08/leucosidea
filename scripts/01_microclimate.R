####Libraries####
library(car)
library(ggforce)
library(ggtext)
library(lme4)
library(MuMIn)
library(patchwork)
library(tidyverse)

source("scripts/_internals.R")

set.seed(66)

####Import data####
env <- left_join(read.table("data/leucosidea_metadata.txt",
                            header = TRUE),
                 read.table("data/microclimate.txt", 
                            header = TRUE)) %>%
  # average the two nearest neighbours as a proxy for 'stand density'
  mutate(density = (nn1 + nn2)/2) %>% 
  # move density column
  relocate(density, .after = ContHei) %>% 
  # remove some columns
  dplyr::select(!c(Wpt, TimeH, nn1, nn2, nna, TempA, CovU, CovC, RichU, RichC)) %>%
  left_join(.,
            fg_cover,
            by = c("Site" = "Site", "Pair" = "Pair")) %>% 
  group_by(Site, Pair) %>% 
  mutate(rich_all_C = sum(richness_forb_C) + sum(richness_grass_C),
         rich_all_U = sum(richness_forb_U) + sum(richness_grass_U),
         cov_all_C = sum(cover_forb_C) + sum(cover_grass_C),
         cov_all_U = sum(cover_forb_U) + sum(cover_grass_U)) %>% 
  ungroup() %>% 
  # data in long format
  pivot_longer(.,
               cols = TempU:cov_all_U,
               names_to = "variable") %>%
  # create variable for microsite (and drop from variable)
  mutate(microsite = as.factor(str_extract(variable, ".{1}$")),
         variable = str_replace(variable, ".{1}$", ""))

###Richness & cover summary statistics
read.table("data/microclimate.txt", 
           header = TRUE) %>%
  select(CovU, RichU, CovC, RichC) %>%
  pivot_longer(everything()) %>%
  mutate(category = case_when(name == "RichU" ~ "richness",
                              name == "CovU" ~ "cover",
                              name == "RichC" ~ "richness",
                              name == "CovC" ~ "cover"),
         plot = case_when(name == "RichU" ~ "under",
                          name == "CovU" ~ "under",
                          name == "RichC" ~ "away",
                          name == "CovC" ~ "away"),
         taxon = "all") %>%
  rbind(fg_cover %>%
          select(-Site, -Pair)%>%
          pivot_longer(everything()) %>%
          mutate(category = case_when(name == "cover_forb_C" ~ "cover",
                                      name == "cover_grass_C" ~ "cover",
                                      name == "cover_forb_U" ~ "cover",
                                      name == "cover_grass_U" ~ "cover",
                                      .default = "richness"),
                 plot = case_when(name == "cover_forb_C" ~ "away",
                                  name == "cover_grass_C" ~ "away",
                                  name == "richness_forb_C" ~ "away",
                                  name == "richness_grass_C" ~ "away",
                                  .default = "under"),
                 taxon = case_when(name == "cover_forb_C" ~ "forb",
                                   name == "cover_forb_U" ~ "forb",
                                   name == "richness_forb_C" ~ "forb",
                                   name == "richness_forb_U" ~ "forb",
                                   .default = "grass"))) %>%
  select(-name) %>%
  group_by(category, taxon) %>%
  summarise(
    mean_under = mean(value[plot == "under"], na.rm = TRUE),
    mean_away  = mean(value[plot == "away"], na.rm = TRUE),
    sd_away    = sd(value[plot == "under"], na.rm = TRUE),
    sd_under   = sd(value[plot == "away"], na.rm = TRUE),
    t_stat     = t.test(value ~ plot)$statistic,
    p.value    = t.test(value ~ plot)$p.value,
    .groups = "drop"
  )


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
                           R_conditional = NA,
                           p_val = NA)

for (i in 1:length(vars)) {
  
  dat <- env %>%
    filter(variable == vars[i])
  
  # model
  lmm <- lmer(formula = value ~ microsite + (1 | Site), 
              dat)
  summ = summary(lmm)$coefficients
  annova = Anova(lmm)
  
  # extract relevant summary stats
  model_results[i, 2] <- summ[1,1]
  model_results[i, 3] <- summ[1,2]
  model_results[i, 4] <- summ[1,3]
  model_results[i, 5] <- summ[2,1]
  model_results[i, 6] <- summ[2,2]
  model_results[i, 7] <- summ[2,3]
  model_results[i, 8] <- r.squaredGLMM(lmm)[1]
  model_results[i, 9] <- r.squaredGLMM(lmm)[2]
  model_results[i, 10] <- annova[,3]
  
}

write.csv(model_results,
          "outputs/mixed_effect_community.csv")

# plot pairwise data

env_plot <- 
  env %>%
  mutate(variable = case_when(variable == 'Temp' ~ "Temperature (°C)",
                              variable == 'STemp' ~ "Soil temperature (°C)",
                              variable == 'SMoist' ~ "Soil moisture (VWC %)",
                              variable == 'cover_forb_' ~ "Forb cover (%)",
                              variable == 'cover_grass_' ~ "Grass cover (%)",
                              variable == 'richness_forb_' ~ "Forb species richness",
                              variable == 'richness_grass_' ~ "Grass species richness",
                              variable == 'rich_all_' ~ "All species richness",
                              variable == 'cov_all_' ~ "All species cover (%)"),
         microsite = case_when(as.character(microsite) == "U" ~ "Under",
                               .default = "Away"))

plots <- vector('list', 3)
grp_vars <- vector('list', 3)
grp_vars[[1]] <- c("Temperature (°C)", "Soil temperature (°C)", "Soil moisture (VWC %)")
grp_vars[[3]] <- c("All species cover (%)", "Forb cover (%)", "Grass cover (%)")
grp_vars[[2]] <- c("All species richness", "Forb species richness", "Grass species richness")

# stars (*) for sig diffs)
starstruck <- data.frame(variable = env_plot %>% distinct(variable),
                         p_val = model_results$p_val) %>%
  full_join(env_plot %>%
              group_by(variable) %>%
              reframe(maxy = max(value))) %>%
  filter(p_val < 0.05)

for (i in 1:length(plots)) {
  
  dat <- env_plot %>% 
    filter(variable %in% grp_vars[[i]])
  
  starstruck_temp <- starstruck %>% 
    filter(variable %in% grp_vars[[i]])
  
  plots[[i]] <- ggplot(dat,
                       aes(x = variable,
                           y = value)) +
    geom_point(data = dat %>%
                 filter(microsite == "Under"),
               aes(x = 1.2,
                   y = value),
               alpha = 0.3,
               fill = "#046A38",
               colour = "white",
               shape = 21,
               position = position_jitternormal(sd_x = 0.05, sd_y = 0)) +
    geom_point(data = dat %>%
                 filter(microsite == "Away"),
               aes(x = 0.8,
                   y = value),
               alpha = 0.3,
               fill = "#FFB81C",
               colour = "white",
               shape = 21,
               position = position_jitternormal(sd_x = 0.05, sd_y = 0)) +
    geom_boxplot(aes(colour = microsite),
                 outliers = FALSE,
                 fill = NA) +
    geom_point(data = starstruck_temp,
               aes(x = 1,
                   y = maxy),
               shape = 8,
               size = 4) +
    facet_wrap(vars(variable),
               scales = "free",
               strip.position = "left") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_colour_manual(values = c("#FFB81C","#046A38"),
                        name = "Microsite") +
    labs(y = NULL,
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

####'full' Mixed effect models####

# this is just to see if a full model is going to have an impact on R2 vals

# create empty data frame
model_results = data.frame(variable = vars, 
                           R_marginal = NA,
                           R_conditional = NA)

for (i in 1:length(vars)) {
  
  dat <- env %>%
    filter(variable == vars[i])
  
  # model
  lmm <- lmer(formula = value ~ microsite + Alt + TrHei + Circ +
                LeuCov + density + (1 | Site),
              dat)
  
  # extract R2 vals
  model_results[i, 2] <- r.squaredGLMM(lmm)[1]
  model_results[i, 3] <- r.squaredGLMM(lmm)[2]
  
}

write.csv(model_results,
          "outputs/full_mixed_effect_community.csv")
