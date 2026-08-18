##### GLMMs for per-species occurrence models ####


####Libraries####
library(tidyverse)
library(glmmTMB)
library(dplyr)
library(performance)

spp <- read.table("data/spp_richness.txt", 
                  header = T) %>% 
  select(-c(Dead, Seedling)) %>%
  # set abundance to 1 if cover > 0
  mutate(across(Ajuga_ophrydis:Tristachya_leucothrix, 
                ~ ifelse(.x == 0, 0 ,1))) %>%
  # create variable for microsite (and drop from variable)
  mutate(pair = str_extract(Plot, "(?<=^[A-I])\\d+(?=[CU]$)") |> as.integer(),
         microsite = str_extract(Plot, ".{1}$")) %>%
  relocate(pair, microsite, .after = Site)
spp <- spp[,-1]

cover <- read.table("spp_richness.txt", 
                  header = T) %>% 
  select(-c(Dead, Seedling)) %>%
  # create variable for microsite (and drop from variable)
  mutate(pair = str_extract(Plot, "(?<=^[A-I])\\d+(?=[CU]$)") |> as.integer(),
         microsite = str_extract(Plot, ".{1}$")) %>%
  relocate(pair, microsite, .after = Site)
cover <- cover[,-1]


# Identify and subset species that meet the criteria for being common enough to include in the analyses
cover <- cover %>%
  select(
    1:3,
    names(spp)[-(1:3)][
      colSums(spp[-(1:3)], na.rm = TRUE) >= 26
    ]
  )

spp <- spp %>%
  select(
    1:3,
    names(spp)[-(1:3)][
      colSums(spp[-(1:3)], na.rm = TRUE) >= 26
    ]
  )

# Defining:
# cover = original species cover data
# spp   = presence/absence version of cover

species_cols <- names(spp)[-(1:3)]
spp[species_cols] <- lapply(
  spp[species_cols],
  function(x) as.integer(x > 0)
)

# Species names
spp_names <- species_cols

# Results dataframe
model_results <- data.frame(
  species = spp_names,
  n_C = NA_integer_,
  n_U = NA_integer_,
  mean_cover_C = NA_real_,
  min_cover_C  = NA_real_,
  max_cover_C  = NA_real_,
  mean_cover_U = NA_real_,
  min_cover_U  = NA_real_,
  max_cover_U  = NA_real_,
  intercept = NA_real_,
  intercept_stderror = NA_real_,
  intercept_tval = NA_real_,
  micrositeU = NA_real_,
  micrositeU_stderror = NA_real_,
  micrositeU_tval = NA_real_,
  R_marginal = NA_real_,
  R_conditional = NA_real_,
  p_val = NA_real_
)


# Loop through species
for(i in seq_along(spp_names)){
  sp <- spp_names[i]
  ## Cover statistics
  cover_C <- cover[[sp]][cover$microsite == "C" & cover[[sp]] > 0]
  cover_U <- cover[[sp]][cover$microsite == "U" & cover[[sp]] > 0]
  model_results$n_C[i] <- length(cover_C)
  model_results$n_U[i] <- length(cover_U)
  if(length(cover_C) > 0){
    model_results$mean_cover_C[i] <- mean(cover_C)
    model_results$min_cover_C[i]  <- min(cover_C)
    model_results$max_cover_C[i]  <- max(cover_C)
  }
  if(length(cover_U) > 0){
    model_results$mean_cover_U[i] <- mean(cover_U)
    model_results$min_cover_U[i]  <- min(cover_U)
    model_results$max_cover_U[i]  <- max(cover_U)
  }
  ## Fit GLMM to presence/absence data
  form <- as.formula(paste(sp, "~ microsite + (1|Site)"))
  mod <- tryCatch(
    glmmTMB(
      form,
      data = spp,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  if(is.null(mod)) next
  coefs <- summary(mod)$coefficients$cond
  model_results$intercept[i]          <- coefs["(Intercept)", "Estimate"]
  model_results$intercept_stderror[i] <- coefs["(Intercept)", "Std. Error"]
  model_results$intercept_tval[i]     <- coefs["(Intercept)", "z value"]
  if("micrositeU" %in% rownames(coefs)){
    model_results$micrositeU[i]          <- coefs["micrositeU", "Estimate"]
    model_results$micrositeU_stderror[i] <- coefs["micrositeU", "Std. Error"]
    model_results$micrositeU_tval[i]     <- coefs["micrositeU", "z value"]
  }
  R2 <- tryCatch(performance::r2_nakagawa(mod), error = function(e) NULL)
  if(!is.null(R2)){
    model_results$R_marginal[i]    <- R2$R2_marginal
    model_results$R_conditional[i] <- R2$R2_conditional
  }
  
  ## Likelihood-ratio test for microsite
  mod0 <- update(mod, . ~ . - microsite)
  model_results$p_val[i] <- anova(mod0, mod)$"Pr(>Chisq)"[2]
}

write.csv(model_results, "mixed_effect_occurrence_site_only.csv")

