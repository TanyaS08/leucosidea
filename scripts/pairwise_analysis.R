# Leucosidea and veg community data
# GGHNp, January 2017
#
# Authors: Tanya Strydom
#          Peter le Roux
# Contact: tanya.strydom@icloud.com

### IMPORTANT NOTICE ----
#' ------------------------------------------------------------------#
#' Please avoid working in section 3) 
#'   Mulivariate analysis of bootstrapped data
#'   Tanya Strydom is working in that section at the moment
#' ------------------------------------------------------------------#

#' ------------------------------------------------------------------#
#'  DATA IMPORTING AND DATAFRAMES
#'  - Import and clean data using Gr3_data_import_checking.R
#'    Can also be found in PFTC5_Gr3 Repo at: 
#'    (https://github.com/TanyaS08/PFTC5_Gr3/blob/master/scripts/Gr3_data_import_checking.R)
#'  - 'species' df = community cover data
#'  - 'traits' df = traits data
#' ------------------------------------------------------------------#

#' ------------------------------------------------------------------#
#'   TO DO:
#'  - Need to update trait variables (when available) that are selected 
#'    when converting to long format -> ca. l. 50
#'  - TBD: hierarchy when imputing/bootstrapping
#'    for now using: Site > Treatment > PlotID
#'        - changing this requires modifying select and gather when 
#'          making the long dfs
#'  - TBD: number of reps and sample size for bootstrapping
#'  - Rank sites from high to low as opposed to alphabetical
#'    for plotting
#'  - if we do decide to plot outputs maybe set better colour scheme
#'    manually
#'  - might be worth exchanging gather for pivot_longer() as gather()
#'    is a depreciated function  
#' ------------------------------------------------------------------#