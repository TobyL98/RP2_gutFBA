#########################
# spear_rank.R
#########################

library(dplyr)
library(tidyr)

setwd('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research Project 2/development/RP2_gutFBA/first_model')
getwd()

abundances_df <- read.csv('Outputs/healthy_df_out.csv', header = T)

abundances_df <- abundances_df %>%
  pivot_longer(!c(Genus, species),
               names_to = "Patient",
               values_to = "Abundance") %>%
  select(!species) %>%
  arrange(Patient)