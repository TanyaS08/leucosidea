####Libraries####


####Import data####
spp <- read.table("data/spp_richness.txt",
                    header = T) %>% 
  select(-c(Dead, Seedling)) %>%
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

ChiSq <- t(spp)
ChiSq <- as.data.frame(ChiSq)
ChiSq$Control_ab <- 135 - ChiSq$away
ChiSq$Under_ab <- 135 - ChiSq$under
tempabun <- which(ChiSq$away + ChiSq$under > 26.5)
length(tempabun)
dat <- ChiSq[tempabun,]

#Empty Matrix to store results
Results <- matrix(data = NA, ncol = 4, nrow = dim(dat)[1])
colnames(Results) <- c("Spp", "Chi2", "df", "p")

#loop for chi test
for (i in seq_len(nrow(dat))) {
  temp.tab <- t(structure(array(as.numeric(dat[i, c(1:4)])), dim = c(2, 2)))
  colnames(temp.tab) <- c("Control", "Under")
  rownames(temp.tab) <- c("Present", "Absent")
  temp.res <- chisq.test(temp.tab) #, simulate.p.value = TRUE)
  Results[i, 1] <- as.character(row.names(dat)[i])
  Results[i, 2] <- round(temp.res$statistic, 4)
  Results[i, 3] <- as.numeric(temp.res$parameter)
  Results[i, 4] <- round(temp.res$p.value, 4)
}

Results

Results <- as.data.frame(Results)
Results$Control <- dat$Control
Results$Under <- dat$Under
Results <- Results[c(1,5,6,2,3,4)]

write.csv(Results, "outputs/chisq_results.csv")