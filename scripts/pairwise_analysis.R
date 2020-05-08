# Leucosidea and veg community data
# GGHNp, January 2017
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
#install.packages("devtools")
library(tidyverse)
library(tidylog)

### >> b) Dataframes ----
#environmental data
env <- read.table(file.path("data", "LeucoENV.txt"), 
                       header = T)




### 1) Bootstrap trait values using CWM ----

env_RII <- (env[,c(13:17)]-env[,c(18:22)])/
  (env[,c(13:17)]+env[,c(18:22)])

#Bootstrapping
library(boot)
bootmean <- function(x, i) mean(x[i]) #function for bootstrap
bootstrap_env <- matrix(NA, nrow = ncol(env_RII), ncol = 3) #empty matrix

for(i in 1:ncol(env_RII)) {
  x <- na.omit(env_RII[[i]]) #remove NA from col
  bss <- boot(x, bootmean, R=10000) #bootstrap
  bc <- boot.ci(bss, type="perc") #CI's for Bootstrap
  t <- bss$t #save simulated means as vector
  bootstrap_env[i,1] <- mean(t) #simulated grand mean
  bootstrap_env[i,2] <- bc$percent[1,4] #lower CI's
  bootstrap_env[i,3] <- bc$percent[1,5] #upper CI's
}

te <- as.data.frame(bootstrap_env)
row.names(te) <- c("Temp","STemp","SMoist","Cover",
                   "Richness")
colnames(te) <- c("Mean", "LCI", "UCI")
y <- ifelse(te$LCI <=0 & te$UCI >=0, F, T) #ID if CI interval crosses zero

GF <- read.delim("GFrichcov.txt", header = T)
GF[GF$Microsite == 'Under',]
GF = cbind(GF[GF$Microsite == 'Under',],GF[GF$Microsite == 'Control',])

GFRII <- (GF[,c(3,4)]-GF[,c(7,8)])/
  (GF[,c(3,4)]+GF[,c(7,8)])
GFRII$type = GF$Type

GFF = GFRII[GFRII$type == 'Forb',c(1,2)]
GFG = GFRII[GFRII$type == 'Grass',c(1,2)]

#Bootstrapping
library(boot)
bootmean <- function(x, i) mean(x[i]) #function for bootstrap
bootstrap_GFF <- matrix(NA, nrow = ncol(GFF), ncol = 3) #empty matrix

for(i in 1:ncol(GFF)) {
  x <- na.omit(GFF[[i]]) #remove NA from col
  bss <- boot(x, bootmean, R=10000) #bootstrap
  bc <- boot.ci(bss, type="perc") #CI's for Bootstrap
  t <- bss$t #save simulated means as vector
  bootstrap_GFF[i,1] <- mean(t) #simulated grand mean
  bootstrap_GFF[i,2] <- bc$percent[1,4] #lower CI's
  bootstrap_GFF[i,3] <- bc$percent[1,5] #upper CI's
}

tf <- as.data.frame(bootstrap_GFF)
colnames(tf) <- c("Mean", "LCI", "UCI")
row.names(tf) <- c("Richness", "Cover")

bootstrap_GFG <- matrix(NA, nrow = ncol(GFG), ncol = 3) #empty matrix

for(i in 1:ncol(GFG)) {
  x <- na.omit(GFG[[i]]) #remove NA from col
  bss <- boot(x, bootmean, R=10000) #bootstrap
  bc <- boot.ci(bss, type="perc") #CI's for Bootstrap
  t <- bss$t #save simulated means as vector
  bootstrap_GFG[i,1] <- mean(t) #simulated grand mean
  bootstrap_GFG[i,2] <- bc$percent[1,4] #lower CI's
  bootstrap_GFG[i,3] <- bc$percent[1,5] #upper CI's
}

tg <- as.data.frame(bootstrap_GFG)
colnames(tg) <- c("Mean", "LCI", "UCI")
row.names(tg) <- c("Richness", "Cover")

other_boot <- rbind(te, tg, tf)
other_boot$var <- c("Air Temperature","Soil Temperature","Soil Moisture",
                    "Species Cover", "Species Richness",
                    "Species Cover", "Species Richness",
                    "Species Cover", "Species Richness")
other_boot$grp <- c(rep("Microclimate",3),
                    rep("Vegetation",2),
                    rep("Grass",2),
                    rep("Forb",2))




# End of script ----