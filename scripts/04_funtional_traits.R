####Libraries####
library(ggforce)
library(ggtext)
library(tidyverse)

####Import data####
FT = read.table("data/FT.txt", header=T) %>%
  mutate(Pair = rep(1:286, each = 2), #unique pair code
         LDMC = Dmass/WMass, 
         SLA = SA/LDMC,
         microsite = case_when(Habitat == "Under" ~ "Under",
                               TRUE ~ "Away"))

FT_long = FT %>%
  select(Site, Species, SA, Chlorophyll, Toughness, PHeight, LDMC, SLA, microsite) %>%
  mutate(Toughness = as.numeric(Toughness),
         PHeight = as.numeric(PHeight)) %>%
  pivot_longer(
    cols = -c(Site, Species, microsite),
    names_to = "trait",
    values_to = "trait_val",
    values_drop_na = TRUE
  )


####Mixed effect models####

# get names of all the variables
traits <- distinct(FT_long, trait) %>%
  pull()

spp_names <- distinct(FT_long, Species) %>%
  pull()

# create empty data frame
model_results = data.frame(species = rep(spp_names, each = length(traits)),
                           trait = rep(traits, times = length(spp_names)), 
                           intercept = NA, 
                           intercpet_stderror = NA, 
                           intercpet_tval = NA,
                           micrositeU = NA, 
                           micrositeU_stderror = NA,
                           micrositeU_tval = NA,
                           R_marginal = NA,
                           R_conditional = NA,
                           p_val = NA)

for (i in 1:length(spp_names)) {
  for (j in 1:length(traits)) {

    index <- j + length(traits)*(i-1)
    
    dat <- FT_long %>%
    filter(Species == spp_names[i]) %>%
    filter(trait == traits[j])

    # model
    lmm <- lmer(formula = trait_val ~ microsite + (1 | Site),
                dat)
    summ = summary(lmm)$coefficients
    annova = Anova(lmm)

    # extract relevant summary stats
    model_results[index, 3] <- summ[1,1]
    model_results[index, 4] <- summ[1,2]
    model_results[index, 5] <- summ[1,3]
    model_results[index, 6] <- summ[2,1]
    model_results[index, 7] <- summ[2,2]
    model_results[index, 8] <- summ[2,3]
    model_results[index, 9] <- r.squaredGLMM(lmm)[1]
    model_results[index, 10] <- r.squaredGLMM(lmm)[2]
    model_results[index, 11] <- annova[,3]
  }
}

write.csv(model_results,
          "outputs/mixed_effect_FT.csv")

# plot 

FT_long <- 
  FT_long %>%
  mutate(Species = case_when(Species == "Pseudognath" ~ "**Pseudognaphalium**",
                             Species == "Oxalis" ~ "**Oxalis**",
                             Species == "Helichrysum" ~ "**Helichrysum**",
                             Species == "Commelina" ~ "**Commelina**",
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = factor(Species, levels=c('**Commelina**', '**Helichrysum**', '**Oxalis**', '**Pseudognaphalium**', 
                                            'Miscanthus', 'Themeda', 'Tristachya')),
         trait = case_when(trait == 'Chlorophyll' ~ "Chlorophyll~(mg/m^2)",
                           trait == 'PHeight' ~ "Plant~height~(m)",
                           trait == 'SLA' ~ "SLA~(cm^2%.%g)",
                           trait == 'LDMC' ~ "LDMC",
                           trait == 'Toughness' ~ "Toughness~(N)",
                           trait == 'SA' ~ "Leaf~surface~area~(cm^2)"))

# stars (*) for sig diffs)
starstruck <-
  data.frame(Species = spp_names,
             trait = traits,
             p_val = model_results$p_val) %>%
  mutate(Species = case_when(Species == "Pseudognath" ~ "**Pseudognaphalium**",
                             Species == "Oxalis" ~ "**Oxalis**",
                             Species == "Helichrysum" ~ "**Helichrysum**",
                             Species == "Commelina" ~ "**Commelina**",
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = factor(Species, levels=c('**Commelina**', '**Helichrysum**', '**Oxalis**', '**Pseudognaphalium**', 
                                            'Miscanthus', 'Themeda', 'Tristachya')),
         trait = case_when(trait == 'Chlorophyll' ~ "Chlorophyll~(mg/m^2)",
                           trait == 'PHeight' ~ "Plant~height~(m)",
                           trait == 'SLA' ~ "SLA~(cm^2%.%g)",
                           trait == 'LDMC' ~ "LDMC",
                           trait == 'Toughness' ~ "Toughness~(N)",
                           trait == 'SA' ~ "Leaf~surface~area~(cm^2)")) %>%
  full_join(FT_long %>%
              group_by(trait) %>%
              reframe(maxy = max(trait_val))) %>%
  filter(p_val < 0.05)

ggplot(FT_long,
       aes(x = Species,
           y = trait_val)) +
  geom_boxplot(aes(colour = microsite),
               outliers = FALSE) +
  geom_point(data = FT_long %>%
               filter(microsite == "Under"),
             aes(x = as.numeric(Species) + 0.2,
                 y = trait_val),
             alpha = 0.3,
             fill = 'forestgreen',
             colour = "white",
             shape = 21,
             position = position_jitternormal(sd_x = 0.05, sd_y = 0)) +
  geom_point(data = FT_long %>%
               filter(microsite == "Away"),
             aes(x = as.numeric(as.factor(Species)) - 0.2,
                 y = trait_val),
             alpha = 0.3,
             fill = 'goldenrod1',
             colour = "white",
             shape = 21,
             position = position_jitternormal(sd_x = 0.05, sd_y = 0)) +
  geom_point(data = starstruck,
             aes(x = Species,
                 y = maxy),
             shape = 8,
             size = 4) +
  facet_wrap(vars(trait),
             scales = "free",
             labeller = label_parsed,
             ncol = 2) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_colour_manual(values = c('goldenrod1','forestgreen'),
                      name = "Microsite") +
  labs(y = "Trait value",
       caption = "Species in **bold** are forbs and non-bold species are grasses") +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(),
    plot.caption = element_markdown(),
    legend.position = 'bottom'
  )

ggsave("figures/FT_boxplot.png",
       width = 8,
       height = 13)
