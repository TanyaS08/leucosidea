####Libraries####
library(ggforce)
library(ggtext)
library(patchwork)
library(tidyverse)
library(vegan)

set.seed(66)

# for spp growth forms
source("scripts/_internals.R")

####Import data####

# array for species community
COMMarray <- read.table("data/spp_richness.txt", header = T, row.names = 1) %>%
  #remove unneeded cols
  select(-c(Site, Dead)) %>% 
  # get microsite from ID
  mutate(Habitat = str_extract(row.names(.), ".{1}$")) %>%
  # standardise naming
  mutate(Habitat = case_when(Habitat == "U" ~ "Under",
                             TRUE ~ "Away")) %>%
  # reorder
  select(Habitat, everything())
#remove NA
COMMarray <- na.omit(COMMarray)

# array for traits
FTarray = read.table("data/FT.txt", header = T, row.names = 1) %>%
  select(Plot, Species, Habitat, Chlorophyll, Toughness, PHeight, SLA, LDMC) %>%
  mutate(Habitat = case_when(Habitat == "Under" ~ "Under",
                             TRUE ~ "Away"),
         Species = case_when(Species == "Pseudognath" ~ "Pseudognaphalium",
                             TRUE ~ Species))
#remove NA
FTarray <- na.omit(FTarray)

####Analysis - Community####

comm_mds <- metaMDS(COMMarray[2:ncol(COMMarray)] %>% 
                      # remove problem (outlier) site
                      filter(!row.names(.) %in% "A11C"), 
                    distance = "bray")

####Plot####

comm_sp_nmds <- as.data.frame(comm_mds$species) %>% 
  mutate(species = row.names(.)) %>% 
  left_join(.,
            growth_forms)  %>% 
  right_join(.,
             spp_fidelity) %>% 
  filter(cover > 0)

comm_site_nmds <- as.data.frame(comm_mds$points) %>%
  mutate(Microsite = COMMarray$Habitat,
         Site = str_extract(row.names(.), "^.{1}"))

ggplot(comm_site_nmds,
       aes(x = MDS1,
           y = MDS2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.8,
             alpha = 0.6,
             colour = "grey50") +
  stat_density_2d(data = comm_sp_nmds,
                  geom = "polygon",
                  aes(x = MDS1,
                      y = MDS2,
                      fill = growth_form,
                      alpha = ..nlevel..),
                  contour_var = "ndensity",
                  breaks = c(0.5, 0.9)) +
  stat_density_2d(data = comm_sp_nmds,
                  geom = "polygon",
                  aes(x = MDS1,
                      y = MDS2,
                      colour = growth_form),
                  contour_var = "ndensity",
                  fill = NA,
                  breaks = c(0.5)) +
  guides(alpha = FALSE,
         fill = NULL) +
  facet_grid(cols = vars(Microsite)) +
  scale_fill_manual(values = c('#01665E', '#BF822E'),
                    name = "Growth form") +
  scale_colour_manual(values = c('#01665E', '#BF822E'),
                      name = "Growth form") +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20)) +
  xlim(-0.65,0.65) +
  ylim(-0.65,0.65) +
  labs(x = "MDS1",
       y = "MDS2",
       caption = paste0( "Stress = ", round(comm_mds$stress*100, digits = 2), "%"))

ggsave("figures/community_pca.png",
       width = 11,
       height = 7)

####PERMANOVA####

permanova_all <- adonis2(COMMarray[2:ncol(COMMarray)] ~ Site*Microsite,
                         data = comm_site_nmds, perm = 999)

write.csv(permanova_all,
          "outputs/permanova_all.csv")

####Analysis - Forbs/Grasses####

FORBarray <-
  COMMarray %>% 
  # remove problem (outlier) site
  filter(!row.names(.) %in% c("A1U", "A11C", "I13C", "A5C", 
                              "A15C", "G13C", "C5U", "G4U",
                              "B1C", "H4U")) %>%
  select(c(growth_forms %>%
             filter(growth_form == "forb") %>%
             filter(species != 'Dead') %>%
             pull(species))) %>%
  # remove plots were there are no forbs
  filter(rowSums(across(where(is.numeric))) != 0)

forb_mds <- metaMDS(FORBarray,
                    distance = "bray")

####Plot####

forb_sp_nmds <- as.data.frame(forb_mds$species) %>%
  mutate(species = row.names(.)) %>%
  left_join(.,
            growth_forms)  %>%
  left_join(.,
            spp_fidelity) %>%
  filter(cover > 0)

forb_site_nmds <- as.data.frame(forb_mds$points) %>%
  mutate(Site = str_extract(row.names(.), "^.{1}")) %>% 
  # get microsite from ID
  mutate(Microsite = str_extract(row.names(.), ".{1}$")) %>%
  # standardise naming
  mutate(Microsite = case_when(Microsite == "U" ~ "Under",
                               TRUE ~ "Away"))

forb_nmds_plot <-
  ggplot(forb_site_nmds,
         aes(x = MDS1,
             y = MDS2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.8,
             alpha = 0.6,
             colour = "grey50") +
  stat_density_2d(data = forb_sp_nmds,
                  geom = "polygon",
                  aes(x = MDS1,
                      y = MDS2,
                      fill = Microsite,
                      alpha = after_stat(nlevel)),
                  contour_var = "ndensity",
                  breaks = c(0.5, 0.9)) +
  stat_density_2d(data = forb_sp_nmds,
                  geom = "polygon",
                  aes(x = MDS1,
                      y = MDS2,
                      colour = Microsite),
                  contour_var = "ndensity",
                  fill = NA,
                  breaks = c(0.5)) +
  guides(alpha = "none",
         fill = NULL) +
  scale_fill_manual(values = c('#01665E', '#BF822E'),
                    name = "Microsite",
                    aesthetics = c("colour", "fill")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20)) +
  labs(x = "MDS1",
       y = "MDS2",
       caption = paste0( "Stress = ", round(forb_mds$stress*100, digits = 2), "%"))

# grasses

GRASSarray <-
  COMMarray %>% 
  # remove problem (outlier) site
  filter(!row.names(.) %in% c("C4C", "C7C","C11U","A10C",
                              "A11C")) %>%
  select(c(growth_forms %>%
             filter(growth_form == "grass") %>%
             pull(species))) %>%
  # remove plots were there are no grasses
  filter(rowSums(across(where(is.numeric)))!= 0)

grass_mds <- metaMDS(GRASSarray,
                     distance = "bray")

####Plot####

grass_sp_nmds <- as.data.frame(grass_mds$species) %>%
  mutate(species = row.names(.)) %>%
  left_join(.,
            growth_forms)  %>%
  left_join(.,
            spp_fidelity) %>%
  filter(cover > 0)

grass_site_nmds <- as.data.frame(grass_mds$points) %>%
  mutate(Site = str_extract(row.names(.), "^.{1}")) %>% 
  # get microsite from ID
  mutate(Microsite = str_extract(row.names(.), ".{1}$")) %>%
  # standardise naming
  mutate(Microsite = case_when(Microsite == "U" ~ "Under",
                               TRUE ~ "Away"))

grass_nmds_plot <-
  ggplot(grass_site_nmds,
         aes(x = MDS1,
             y = MDS2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.8,
             alpha = 0.6,
             colour = "grey50") +
  stat_density_2d(data = grass_sp_nmds,
                  geom = "polygon",
                  aes(x = MDS1,
                      y = MDS2,
                      fill = Microsite,
                      alpha = after_stat(nlevel)),
                  contour_var = "ndensity",
                  breaks = c(0.5, 0.9)) +
  stat_density_2d(data = grass_sp_nmds,
                  geom = "polygon",
                  aes(x = MDS1,
                      y = MDS2,
                      colour = Microsite),
                  contour_var = "ndensity",
                  fill = NA,
                  breaks = c(0.5)) +
  guides(alpha = "none",
         fill = NULL) +
  scale_fill_manual(values = c('#01665E', '#BF822E'),
                    name = "Microsite",
                    aesthetics = c("colour", "fill")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20)) +
  labs(x = "MDS1",
       y = "MDS2",
       caption = paste0( "Stress = ", round(grass_mds$stress*100, digits = 2), "%"))

forb_nmds_plot + 
  labs(tag = "A") + 
  grass_nmds_plot + 
  labs(tag = "B") +
  plot_layout(ncol = 1,
              guides = 'collect') +
  plot_annotation(theme = theme(
    legend.position = 'bottom'))

ggsave("figures/forb_grass_pca.png",
       width = 11,
       height = 14)

####Analysis - Functional traits####

ft_mds <- metaMDS(FTarray[4:ncol(FTarray)], distance = "bray")

####Plot####

ft_traits_nmds <- as.data.frame(ft_mds$species) %>% 
  mutate(traits = row.names(.))

ft_species_nmds <- as.data.frame(ft_mds$points) %>%
  mutate(Microsite = FTarray$Habitat,
         Site = str_extract(row.names(.), "^.{1}"),
         species = FTarray$Species)

ft_spp_colours = tibble(
  species = c("Commelina", "Helichrysum","Miscanthus", "Oxalis",
              "Pseudognaphalium", "Themeda", "Tristachya"),
  colour = c('#003C2F','#01665E', '#543006', '#369890',
             '#7FCEC2', '#BF822E',  '#8D5108'),
  fgroup = c("Forb", "Forb", "Grass", "Forb",
             "Forb", "Grass", "Grass")
)

####Plot####

ggplot(ft_species_nmds,
       aes(x = MDS1,
           y = MDS2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.8,
             alpha = 0.6,
             colour = "grey50") +
  stat_density_2d(geom = "polygon",
                  aes(fill = species,
                      alpha = after_stat(nlevel)),
                  contour_var = "ndensity",
                  breaks = c(0.5, 0.9)) +
  stat_density_2d(geom = "polygon",
                  aes(colour = species),
                  contour_var = "ndensity",
                  fill = NA,
                  breaks = c(0.5)) +
  guides(alpha = FALSE,
         fill = NULL) +
  facet_grid(cols = vars(Microsite)) +
  scale_fill_manual(values = ft_spp_colours$colour,
                    aesthetics = c("colour", "fill")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20)) +
  xlim(-0.4,0.4) +
  ylim(-0.4,0.4) +
  labs(x = "MDS1",
       y = "MDS2",
       caption = paste0( "Stress = ", round(ft_mds$stress*100, digits = 2), "%"))

ggsave("figures/FT_pca.png",
       width = 11,
       height = 8)
