# Leucosidea and veg community data
# GGHNP, January 2017
#
# Authors: Tanya Strydom
#          Peter le Roux
# Contact: tanya.strydom@icloud.com

#' ------------------------------------------------------------------#
#'  DATA IMPORTING AND DATAFRAMES
#'  
#' ------------------------------------------------------------------#

#' ------------------------------------------------------------------#
#'   TO DO:
#'   
#' ------------------------------------------------------------------#
#' 
#' ### 0) Preamble ----
### >> a) Dependencies ----
install.packages("devtools")
if(!require(traitstrap)){ # for bootstrapping trait data 
  devtools::install_github("richardjtelford/traitstrap")
  library(traitstrap)
}
library(tidyverse)
library(tidylog)

### >> b) Dataframes ----
#environmental data
FT <- read.table(file.path("data", "FT.txt"),
                 header = T)



# End of script ----