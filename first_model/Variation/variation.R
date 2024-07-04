# Variation.R

library(dplyr)

setwd('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research Project 2/development/RP2_gutFBA/first_model')
getwd()

abundances_df <- read.csv('Outputs/healthy_df_out.csv', header = T)

# calculates the mean and sd for each row where every column is numeric
abundances_df <- abundances_df %>%
  rowwise() %>%
  mutate(
    mean = mean(c_across(where(is.numeric))),
    sd = sd(c_across(where(is.numeric)))
  )


colnames(abundances_df)
