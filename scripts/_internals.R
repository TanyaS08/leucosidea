####Libraries####
library(tidyverse)

####Assign species as forb/grass####

dat <- read.table("data/spp_richness.txt", header = T, row.names = 1)

#create species list from column names
spp_list = colnames(dat[2:ncol(dat)])

# grasses species list
grass_sp = c("Andropogon_appendiculatus", 
             "Andropogon_schirensis",
             "Aristida_adscensionis",
             "Aristida_junciformes", 
             "Brachiaria_serrata", 
             "Cymbopogon_pospischilii", 
             "Digitaria_monodactyla", 
             "Elionurus_muticus", 
             "Eragrostis_capensis",
             "Eragrostis_chloromelas",
             "Eragrostis_curvula",
             "Eragrostis_lehmanniana",
             "Eragrostis_plana",
             "Eragrostis_purple",
             "Eragrostis_trichophora",
             "Harpochloa_flax",
             "Helictotrichon_turgidulum",
             "Heteropogon_contortus",
             "Hyparrhenia_tamba",
             "Loudetia_simplex",
             "Miscanthus_capensis",
             "Setaria_sphacelata",
             "Themeda_triandra",
             "Tristachya_leucothrix")

growth_forms = tibble(growth_form = ifelse(spp_list %in% grass_sp,
                                           "grass",
                                           "forb"), 
                      species = spp_list)

####Functional group cover/richness####

# richness at each plot for grasses and forbs

fg_cover <-
  dat %>% 
  # get microsite from ID
  mutate(Microsite = str_extract(row.names(.), ".{1}$"),
         Pair = rep(1:135, 2)) %>%
  # standardise naming
  pivot_longer(-c(Microsite, Pair, Site),
               names_to = "species",
               values_to = "cover") %>% 
  left_join(.,
            growth_forms) %>% 
  filter(species != "Dead") %>% 
  select(-species) %>% 
  group_by(Site, Pair, Microsite, growth_form) %>% 
  mutate(richness = ifelse(cover > 0, 1, 0)) %>% 
  summarise(richness = sum(richness),
            cover = sum(cover)) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = c(growth_form, Microsite),
    values_from = c(cover, richness)
  ) %>% 
  ungroup()


####Site fidelity fo species####

# basically see if a species occurs in only one or both microsites

spp_fidelity <- dat %>% 
  # get microsite from ID
  mutate(Microsite = str_extract(row.names(.), ".{1}$")) %>%
  # standardise naming
  mutate(Microsite = case_when(Microsite == "U" ~ "Under",
                               TRUE ~ "Away")) %>% 
  group_by(Microsite) %>%
  summarise(across(Ajuga_ophrydis:Tristachya_leucothrix, 
                   ~ sum(.x))) %>% 
  pivot_longer(-Microsite,
               names_to = "species",
               values_to = "cover")  %>% 
  filter(!species == "Dead")

#remove for cleaner working environment
rm(dat, spp_list, grass_sp)