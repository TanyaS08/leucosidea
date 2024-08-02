####Libraries####
library(tidyverse)
library(ggbeeswarm)
library(ggforce)
library(ggtext)
library(vegan)
library(ade4)
library(ggrepel)

####Import data####

FTarray = read.table("data/FT_NMDS.txt",header=T, row.names = 1)
#remove NA
FTarray <- na.omit(FTarray)
FTarray$Reproductive <- NULL
FTarray =
  FTarray %>%
  mutate(Habitat = case_when(Habitat == "Under" ~ "Under",
                             TRUE ~ "Away"))

dat = FTarray[5:ncol(FTarray)]
Blob_sp = FTarray[2]
Blob_microsite = FTarray[3]

titre <- c(Chlorophyll = "Chlorophyll~(mg/m^2)", PHeight = "Plant~height~(m)", SLA = "SLA~(cm^2%.%g)",
           LDMC = "LDMC", Toughness = "Toughness~(N)")[colnames(dat)]

####Analysis####

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

mult <- 4
arrows_all = tibble(
  Trait = titre, 
  xend = signx * PCA$co[, 'Comp1'] * mult,
  yend = signy * PCA$co[, 'Comp2'] * mult,
  x = rep(0, length(titre)),
  y = rep(0, length(titre)))

####Plot####

plot_data_all = tibble(Axis1 = x, 
                       Axis2 = y,
                       FGroup = Blob_sp$Species,
                       Microsite = Blob_microsite$Habitat)

ggplot(plot_data_all,
       aes(x = Axis1,
           y = Axis2)) +
  scale_alpha(range = c(0.3, 0.7)) +
  geom_point(size = 0.2,
             alpha = 0.6,
             colour = "grey50") +
  stat_density_2d(geom = "polygon",
                  aes(fill = FGroup,
                      alpha = ..nlevel..),
                  contour_var = "ndensity",
                  breaks = c(0.5, 0.9))+
  stat_density_2d(geom = "polygon",
                  aes(colour = FGroup),
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
               size = 0.2) +
  geom_label_repel(data = arrows_all,
                   aes(label = Trait,
                       x = xend,
                       y = yend),
                   hjust = 0.1,
                   parse = TRUE, 
                   size = 2,
                   check_overlap = FALSE) +
  scale_fill_manual(values = c('#003C2F','#01665E','#543006','#369890',
                               '#7FCEC2','#8D5108','#BF822E')) +
  scale_colour_manual(values = c('#003C2F','#01665E','#543006','#369890',
                                 '#7FCEC2','#8D5108','#BF822E')) +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_text(size = 20)) +
  xlim(-5,5) +
  ylim(-5,5) +
  labs(x = "PC1",
       y = "PC2") +
  theme(legend.position = 'none')

ggsave("figures/FT_pca.png",
       width = 11,
       height = 6)