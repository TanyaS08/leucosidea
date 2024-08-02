####Libraries####
library(lme4)
library(MuMIn)
library(tidyverse)

####Import data####
env <- read.table("data/microclimate.txt", 
                  header = T) %>%
  # remove some columns
  dplyr::select(!c(Wpt, Alt, TimeH, nn1, nn2, TempA)) %>%
  # data in long format
  pivot_longer(.,
               cols = TempU:RichC,
               names_to = "variable") %>%
  # create variable for microsite (and drop from variable)
  mutate(microsite = as.factor(str_extract(variable, ".{1}$")),
         variable = str_replace(variable, ".{1}$", ""))

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
                       R_conditional = NA)

for (i in 1:length(vars)) {
  
  dat <- env %>%
    filter(variable == vars[i])
  
  # model
  lmm <- lmer(formula = value ~ microsite + (1 | Site), 
       dat)
  summ = summary(lmm)$coefficients
  
  # extract relevant summary stats
  model_results[i, 2] <- summ[1,1]
  model_results[i, 3] <- summ[1,2]
  model_results[i, 4] <- summ[1,3]
  model_results[i, 5] <- summ[2,1]
  model_results[i, 6] <- summ[2,2]
  model_results[i, 7] <- summ[2,3]
  model_results[i, 8] <- r.squaredGLMM(lmm)[1]
  model_results[i, 9] <- r.squaredGLMM(lmm)[2]
  
}

write.csv(model_results,
          "outputs/mixed_effect.csv")