####Libraries####
library(ggforce)
library(ggtext)
library(tidyverse)

####Import data####
FT = read.table("data/FT.txt",header=T) %>%
  mutate(Pair = rep(1:286, each = 2), #unique pair code
         LDMC = Dmass/WMass, 
         SLA = SA/LDMC,
         microsite = case_when(Habitat == "Under" ~ "Under",
                               TRUE ~ "Away"))

FT_long = FT %>%
  select(Species, SA, Chlorophyll, Toughness, PHeight, LDMC, SLA, microsite) %>%
  mutate(Toughness = as.numeric(Toughness),
         PHeight = as.numeric(PHeight)) %>%
  pivot_longer(
    cols = -c(Species, microsite),
    names_to = "trait",
    values_to = "trait_val",
    values_drop_na = TRUE
  ) %>%
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


####t-test####

FT_ttest = 
  FT_long %>%
  group_by(trait, Species) %>%
  group_split() %>%
  map_dfr(. %>%
            summarise(trait = trait[1],
                      Species = Species[1],
                      p_val = t.test(trait_val ~ microsite, data = .)$p.value)) %>%
  full_join(FT_long %>%
              group_by(trait) %>%
              reframe(maxy = max(trait_val))) %>%
  filter(p_val < 0.05)

# plot 

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
  geom_point(data = FT_ttest,
             aes(x = Species,
                 y = maxy),
             shape = 8,
             size = 4) +
  facet_wrap(vars(trait),
             scales = "free",
             labeller = label_parsed) +
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
       width = 13,
       height = 8)