# Leucosidea and veg community data
# GGHNP, January 2017
#
# Authors: Tanya Strydom
#          Peter le Roux
# Contact: tanya.strydom@icloud.com

#' ------------------------------------------------------------------#
#'  ISSUES
#'  Traitsstrap only does community weighted - i.e. would still have 
#'  to bootstrap trait values if doing by species analyses
#' ------------------------------------------------------------------#

#' ------------------------------------------------------------------#
#'   TO DO:
#'   Species level bootstrapping
#'   Include pariwise analysis for species
#' ------------------------------------------------------------------#
#' 
#' ### 0) Preamble ----
### >> a) Dependencies ----
#install.packages("devtools")
if(!require(traitstrap)){ # for bootstrapping trait data 
  devtools::install_github("richardjtelford/traitstrap")
  library(traitstrap)
}
library(tidyverse)
library(tidylog)
library(FactoMineR)
library(factoextra)

### >> b) Dataframes ----
#environmental data
traits <- read.table(file.path("data", "FT.txt"),
                 header = T)

### 1) Bootstrap trait values using CWM ----

### >> a) traits df into long format ---- 
traits_long <-
  traits %>%
  #select relevant variables
  ##STILL NEED TO ADD ALL TRAIT VARIABLES
  select(.,
         Site,
         Plot,
         Species,
         SA,
         Chlorophyll,
         Toughness,
         PHeight,
         LDMC,
         SLA) %>%
  #pivot into long format
  
  pivot_longer(.,
               #columns not to be pivoted but are 'retained' i.e. not the traits
               cols = -c(Site,
                         Plot,
                         Species),
               #traits column
               names_to  = "Trait",
               #traits value column
               values_to = "Value") %>%
  
  #remove NA values to avoid sampling
  na.omit()

### >> b) species cover in long format ----

#since we have no cover data assigning same cover vals for all species

species_long <-
  tibble(
    #call sites
    Site = traits_long$Site,
    #call plots
    Plot = traits_long$Plot,
    #call species
    Species = traits_long$Species,
    #assign same cover value
    Cover = rep(0.3,
                nrow(traits_long))
  )

### >> c) Bootstrapping calculations ----
## We need to decide what out 'levels' or heirarchy are.
##For now using Site > Treatment > PlotID

trait_bootstrap <-
  #This imputes the  trait data 
  trait_fill(
    species_long,
    traits_long,
    scale_hierarchy = c("Plot",
                        "Site"),
    taxon_col = "Species",
    trait_col = "Trait",
    value_col = "Value",
    abundance = "Cover",
    keep_all = TRUE) %>%
  #This bootstraps the imputed data
  trait_np_bootstrap(.,
                     nrep = 1000) #TBD


#Summary of different moments (\mu, lower CI and upper CI)
trait_bootstrap_summary <-
  trait_summarise_boot_moments(trait_bootstrap)

### 2) Plotting Boot Moments ----

### >> a) Density plots ----
##Mean trait value
trait_bootstrap %>%
  ggplot(aes(mean, fill = Plot)) +
  facet_wrap(vars(Trait),
             scales = 'free') +
  geom_density(alpha = .5, kernel = "gaussian") +
  scale_fill_brewer(palette = "RdYlGn")  +
  theme_bw() +
  labs(y = "density")

##Variance
trait_bootstrap %>%
  ggplot(aes(variance, fill = Plot)) +
  facet_wrap(vars(Trait),
             scales = 'free') +
  geom_density(alpha = .5, kernel = "gaussian") +
  scale_fill_brewer(palette = "RdYlGn")  +
  theme_bw() +
  labs(y = "density")

## This code can be modified and tweaked for the other
# moments if desired

### >> b) Whisker plots ----

##Mean trait value
trait_bootstrap_summary %>%
  ggplot() +
  #split into panels by trait type and site
  facet_grid(cols = vars(Trait),
             scales = 'free_x',) +
  geom_pointrange(aes(x = Plot,
                      y = mean,
                      ymin = ci_low_mean,
                      ymax  = ci_high_mean,
                      colour = Plot)) +
  scale_colour_brewer(palette = "RdYlGn")  +
  theme_bw()

## This code can probably be modified and tweaked for the other
# moments if so desired

### 3) Mulivariate analysis of bootstrapped data ----

### >> a) Manipulate df onto wide format ----

traits_wide  <-
  trait_bootstrap %>%
  ungroup() %>%
  #remove all moments except for mean
  select(-c(variance,
            skewness,
            kurtosis)) %>%
  #change into wide form i.e. each trait becomes a column
  pivot_wider(.,
              names_from = Trait,
              values_from = mean) %>%
  #need to create a unique ID for each entry - make row name
  mutate(row_ID = paste0(n, "_", Site, "_", Treatment, "_", PlotID)) %>%
  #set new column as row name
  column_to_rownames("row_ID")

### >> b) PCA ----

trait.pca <- PCA(
  #remove identification columns
  traits_wide[,-c(1:5)],
  graph = FALSE,
  scale = TRUE)

fviz_pca(trait.pca,
         label = "var",
         #create a grouping var consisting of treatment and site
         habillage = as.factor(paste0(traits_wide$Site,
                                      "_",
                                      traits_wide$Treatment)),
         addEllipses = TRUE,
         #using the selected colour scheme for consistency
         #NEED TO UPDATE COLOURS WHEN ALL SITES ARE INCLUDED
         palette = c(acj_bb,
                     que_bb,
                     tre_bb,
                     tre_c))

# End of script ----
