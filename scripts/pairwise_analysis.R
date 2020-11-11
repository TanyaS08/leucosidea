# Leucosidea and veg community data
# GGHNP, January 2017
#
# Authors: Tanya Strydom
#          Peter le Roux
# Contact: tanya.strydom@icloud.com

#' ------------------------------------------------------------------#
#'   TODO:
#'   
#' ------------------------------------------------------------------#
#' 
#' ### 0) Preamble ----
### >> a) Dependencies ----
#install.packages("devtools")
library(tidyverse)
library(tidylog)
library(lmerTest)

### >> b) Dataframes ----
#environmental data
env <- read.table(file.path("data", "LeucoENV.txt"), 
                       header = T)
GrassForb <- read.table(file.path("data", "GFrichcov.txt"), 
                        header = T)

### 1) Mixed effect models ----

x <- 
  rbind(
    env %>%
      #remove unwanted variables
      select(-c(Wpt, Alt, TimeH, TrHei, Circ, LeuCov, ContHei, nn1, nn2, TempA, nna)) %>%
      #pivor into long format
      pivot_longer(cols = -c(Site, Pair),
                   names_to = "variable") %>%
      #create microsite variable
      mutate(Microsite = ifelse(str_detect(variable, "U"),
                                "Under",
                                "Away"),
             #remove U and C from env variable
             variable = as.factor(str_replace(variable,
                                              "[[U, C]]*$",
                                              ""))),
    GrassForb %>%
      mutate(Site = rep(rep(LETTERS[1:9], each=15), 4),
             Pair = rep(rep(1:15, 9), 4),
             Microsite = ifelse(Microsite == "Control",
                                "Away",
                                "Under")) %>%
      pivot_wider(names_from = Type,
                  values_from = c(Richness, Cover)) %>%
      pivot_longer(-c(Site, Pair, Microsite),
                   names_to = "variable") %>%
      select(Site, Pair, variable, value, Microsite))
  
#loop over the different environmental variables
for (i in levels(x$variable)) {
  #subset of df by environ variable
  dat = x[x$variable==i, ]
  
  print(summary(lmer(value ~ Microsite + (1 | Site),
      data = dat)))
  
}
  

# End of script ----