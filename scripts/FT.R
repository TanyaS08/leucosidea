####Libraries####
library(tidyverse)
library(ggbeeswarm)
library(ggforce)
library(ggtext)
library(vegan)
library(ade4)
library(ggrepel)

#Import dataset
FT = read.table("data/FT.txt",header=T) %>%
  mutate(Pair = rep(1:286, each = 2), #unique pair code
         LDMC = Dmass/WMass, 
         SLA = SA/LDMC,
         microsite = case_when(Habitat == "Under" ~ "Under",
                               TRUE ~ "Away"))

####PLOTS per SPP####

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


# Quick n dirty t-test

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
  geom_boxplot(aes(colour = microsite)) +
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


#### PCA ####

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

####COEFFICIENT OF VARIATION####
#libraries
library(plyr)
library(reshape2)
#now manipulate
melted <- melt(FT, id.vars=c("Plot", "Species"), 
               measure.vars=c('PHeight', 'Chlorophyll', 'Toughness',
                              'LDMC', 'SLA'))
melted2 <- na.omit(melted)
melted2$value <- as.numeric(melted2$value)


CV <- ddply(melted2, c("Plot", "Species", "variable"), summarise,
            mean = mean(value), sd = sd(value),
            sem = sd(value)/sqrt(length(value)))

CV$CV <- CV$sd/CV$mean
CV$CV2 <- CV$mean/CV$sd
write.csv(CV, "output/CV.csv")

CV2 <- na.omit(CV)

library(ggplot2)
ggplot(CV2, aes(Species, CV2)) +   
  geom_bar(aes(fill = Plot), position = "dodge", stat="identity")

###

####QUICK N DIRTY GLZs####
FT <- na.omit(FT)
FT$LDMC <- as.numeric(FT$LDMC)
FT$PHeight <- as.numeric(FT$PHeight)

PlotModelC <- glm(FT$Chlorophyll ~ FT$Plot, family = 'gaussian')
PlotModelSLA <- glm(FT$SLA ~ FT$Plot, family = 'gaussian')
PlotModelT <- glm(FT$Toughness ~ FT$Plot, family = 'gaussian')
PlotModelPH <- glm(FT$PHeight ~ FT$Plot, family = 'gaussian')
PlotModelLDMC <- glm(FT$LDMC ~ FT$Plot, family = 'gaussian')


