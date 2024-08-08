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
COMMarray <- read.table("data/spp_richness.txt", header=T, row.names = 1) %>%
  #remove unneeded cols
  select(-c(Site, Dead)) %>% 
  # remove problem site
  filter(!row.names(.) %in% "A11C") %>% 
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

comm_mds <- metaMDS(COMMarray[2:ncol(COMMarray)], distance = "bray") 

####Plot####

comm_sp_nmds <- as.data.frame(comm_mds$species) %>% 
  mutate(species = row.names(.)) %>% 
  left_join(.,
            growth_forms)  %>% 
  right_join(.,
             spp_fidelity) %>% 
  filter(cover > 0)

comm_site_nmds <- as.data.frame(comm_mds$points) %>% 
  mutate(Microsite = COMMarray$Habitat) 

ggplot(comm_site_nmds,
       aes(x = MDS1,
           y = MDS2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.2,
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
  xlim(-1,1) +
  ylim(-1,1) +
  labs(x = "MDS1",
       y = "MDS2",
       caption = paste0( "Stress = ", round(comm_mds$stress*100, digits = 2), "%"))

ggsave("figures/community_pca.png",
       width = 11,
       height = 7)


####Analysis - Functional traits####

dat = FTarray[4:ncol(FTarray)]
Blob_sp = FTarray$Species
Blob_microsite = FTarray$Habitat

PCA <- ade4::dudi.pca(dat, center = T, scale = T, scannf = F, nf = ncol(dat))
means<- PCA$cent
sds <- PCA$norm
eigen <- round(PCA$eig / sum(PCA$eig) * 100, 1)

signx <- ifelse(sign(sum(sign(PCA$c1[, "CS1"]))) < 0, -1 ,1)
signy <- ifelse(sign(sum(sign(PCA$c1[, "CS2"]))) < 0, -1, 1)
x = signx * PCA$li[, 'Axis1']
y = signy * PCA$li[, 'Axis2']
Axis1 <- signx * PCA$li[, 'Axis1']
Axis2 <- signy * PCA$li[, 'Axis2']

titre <- c(Chlorophyll = "Chlorophyll~(mg/m^2)", PHeight = "Plant~height~(m)", SLA = "SLA~(cm^2%.%g)",
           LDMC = "LDMC", Toughness = "Toughness~(N)")[colnames(dat)]

mult <- 4
arrows_all = tibble(
  Trait = titre, 
  xend = signx * PCA$co[, 'Comp1'] * mult,
  yend = signy * PCA$co[, 'Comp2'] * mult,
  x = rep(0, length(titre)),
  y = rep(0, length(titre)))

ft_spp_colours = tibble(
  species = c("Commelina", "Helichrysum","Miscanthus", "Oxalis",
              "Pseudognaphalium", "Themeda", "Tristachya"),
  colour = c('#003C2F','#01665E', '#543006', '#369890',
             '#7FCEC2', '#BF822E',  '#8D5108'),
  fgroup = c("Forb", "Forb", "Grass", "Forb",
             "Forb", "Grass", "Grass")
)

plot_data_all = tibble(Axis1 = x, 
                       Axis2 = y,
                       Species = Blob_sp,
                       Microsite = Blob_microsite,
                       fgroup = ifelse(Blob_sp %in% c("Miscanthus", "Themeda", "Tristachya"),
                                       "Grass",
                                       "Forb"))

####Plot####

ggplot(plot_data_all,
         aes(x = Axis1,
             y = Axis2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.2,
             alpha = 0.6,
             colour = "grey50") +
  stat_density_2d(geom = "polygon",
                  aes(fill = Species,
                      alpha = ..nlevel..),
                  contour_var = "ndensity",
                  breaks = c(0.5, 0.9))+
  stat_density_2d(geom = "polygon",
                  aes(colour = Species),
                  contour_var = "ndensity",
                  fill = NA,
                  breaks = c(0.5)) +
  guides(alpha = FALSE,
         fill = NULL) +
  facet_grid(cols = vars(Microsite)) +
  geom_segment(data = arrows_all,
               aes(xend = xend, 
                   yend = yend,
                   x = x,
                   y = y),
               arrow = arrow(length=unit(0.2,"cm")),
               size = 0.3) +
  geom_label_repel(data = arrows_all,
                   aes(label = Trait,
                       x = xend,
                       y = yend),
                   hjust = 0.1,
                   parse = TRUE, 
                   size = 3,
                   check_overlap = FALSE) +
  scale_fill_manual(values = ft_spp_colours$colour,
                    aesthetics = c("colour", "fill")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20)) +
  xlim(-4.3,4.3) +
  ylim(-4.3,4.3) +
  labs(x = "PC1",
       y = "PC2")

ggsave("figures/FT_pca.png",
       width = 11,
       height = 8)
