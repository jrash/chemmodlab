## code to prepare `pparg` dataset goes here

library(dplyr, warn.conflicts = FALSE)
library(rvest)
library(tidyr)
library(readr)
library(usethis)

# data file obtained from http://stats.drugdesign.fr on July 19, 2020

df <- read.csv("~/RESEARCH/JeremyAsh/RecallPaper/2022.01.13/screen_explorer_data.csv")

# Compile dataframe where records are ordered by compound ID for all 3 docking 
#   scores (surf, icm, vina).
# Rescale docking scores so that larger values suggest active compounds. 

pparg <- data.frame(df[order(df$surf_id),1:3], df[order(df$icm_id),4:6], df[order(df$vina_id),7:9] )
pparg$icm_scores <- -pparg$icm_scores
pparg$vina_scores <- -pparg$vina_scores

# Create minimum rank consensus of surf & icm 

order_surf <- order(pparg$surf_scores, decreasing = TRUE)
order_icm <- order(pparg$icm_scores, decreasing = TRUE)

n <- nrow(pparg)

minrank.mat <- matrix(ncol = 2, nrow = n)
minrank.mat[order_surf, 1] <- 1:n
minrank.mat[order_icm, 2] <- 1:n

min.rank <- apply(minrank.mat, 1, min)  
pparg$minr_id <- pparg$surf_id
pparg$minr_scores <- min.rank
pparg$minr_actives <- pparg$surf_actives

# Create maximum z-score consensus of surf & icm

z.mat <- matrix(ncol = 2, nrow = n)
z.mat[, 1] <- scale(pparg$surf_scores, center = T, scale = T)
z.mat[, 2] <- scale(pparg$icm_scores, center = T, scale = T)

z.score <- apply(z.mat, 1, max)    
pparg$maxz_id <- pparg$surf_id
pparg$maxz_scores <- z.score
pparg$maxz_actives <- pparg$surf_actives

# Output data files
write_csv(pparg, "data-raw/pparg.csv")
usethis::use_data(pparg, overwrite = TRUE)
