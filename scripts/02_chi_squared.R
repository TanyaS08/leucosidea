####Libraries####


####Import data####
spp <- read.table("data/spp_richness.txt",
                    header = T) %>%
  # set abundance to 1 if cover > 0
  mutate(across(Ajuga_ophrydis:Tristachya_leucothrix, 
                ~ ifelse(.x == 0, 0 ,1))) %>%
  # create variable for microsite (and drop from variable)
  mutate(microsite = str_extract(Plot, ".{1}$")) %>%
  group_by(microsite) %>%
  summarise(across(Ajuga_ophrydis:Tristachya_leucothrix, 
                   ~ sum(.x))) %>%
  ungroup()
  
# filter out the low abundance species (i.e. not occurring in more than 27 plots)
spp <- spp[,-1][colSums(spp[,-1]) > 26] %>%
  mutate(microsite = c("away", "under")) %>%
  column_to_rownames(var = "microsite")

####Chi2 test####

chi <- chisq.test(spp,
                  correct = F)
