################
# individual_variation
################

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(grid)
library(gridExtra)

####
# Functions
####
colour_function <- function(df) {
  # function assigns black to all the samples and red
  # to the average
  df <- df %>%
    mutate(colour = if_else(sample == 'average', 'average', 'sample'))
  return(df)
}

plot_function <- function(healthy_df, CRC_df, title_name, tag_name) {
  # plots the point plot to compare healthy and CRC
  ggplot(healthy_df,
        aes(x = flux, 
            y = 'Healthy',
            color = colour)) +
    geom_point(shape = 'X', size = 5) +
    geom_point(data = CRC_df, aes(x = flux,
                                  y = ' CRC',
                                  color = colour),
               size = 5) +
    scale_color_manual(values = c("red", "black")) +
    labs(title = title_name,
         tag = tag_name) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(hjust = 0.7),
          title = element_text(size = 18),
          legend.position = 'none'
    )
}

organise_function <- function(fp) {
  # reads in and organises datframe before plotting
  flux_df <- read.csv(fp, header = T) %>%
    colour_function() %>%
    mutate(flux = flux * -1)
  return(flux_df)
  }

# plot community flux for healthy and CRC
healthy_community_flux_df <- 
  read.csv('variation_healthy/obj_values.csv'
           , header = T) %>%
  colour_function()

CRC_community_flux_df <-
  read.csv('variation_CRC/objective_values_crc.csv'
           , header = T) %>%
  colour_function()
community_flux_plot <- plot_function(healthy_community_flux_df,
                                     CRC_community_flux_df,
                                     title_name = "Community Biomass",
                                     tag_name = 'A')

# plot acetate flux for healthy and CRc
ac_healthy_flux_path <- 'variation_healthy/variation_EX_ac_medium_healthy.csv'
ac_healthy_flux_df <- organise_function(ac_healthy_flux_path)

ac_CRC_flux_path <- 'variation_crc/variation_EX_ac_medium_crc.csv'
ac_CRC_flux_df <- organise_function(ac_CRC_flux_path)
acetate_flux_plot <- plot_function(ac_healthy_flux_df,
                                   ac_CRC_flux_df,
                                   title_name = 'Acetate',
                                   tag_name = 'B')

# plot alanine flux for healthy and CRC
ala_healthy_flux_path <- 'variation_healthy/variation_EX_ala_L_medium_healthy.csv'
ala_healthy_flux_df <- organise_function(ala_healthy_flux_path)
ala_CRC_flux_path <- 'variation_crc/variation_EX_ala_L_medium_crc.csv'
ala_CRC_flux_df <- organise_function(ala_CRC_flux_path)
alanine_flux_plot <- plot_function(ala_healthy_flux_df,
                                   ala_CRC_flux_df,
                                   title_name = 'L-alanine',
                                   tag_name = 'C')

# plot propionate flux for healthy and CRC
ppa_healthy_flux_path <- 'variation_healthy/variation_EX_ppa_medium_healthy.csv'
ppa_healthy_flux_df <- organise_function(ppa_healthy_flux_path)
ppa_CRC_flux_path <- 'variation_crc/variation_EX_ppa_medium_crc.csv'
ppa_CRC_flux_df <- organise_function(ppa_CRC_flux_path)
propionate_flux_plot <- plot_function(ppa_healthy_flux_df,
                                      ppa_CRC_flux_df,
                                      title_name = 'Propionate',
                                      tag_name = 'D')

rbflvrd_healthy_flux_path <- 'variation_healthy/variation_EX_rbflvrd_medium_healthy.csv'
rbflvrd_healthy_flux_df <- organise_function(rbflvrd_healthy_flux_path)
rbflvrd_CRC_flux_path <- 'variation_CRC/variation_EX_rbflvrd_medium_crc.csv'
rbflvrd_CRC_flux_df <- organise_function(rbflvrd_CRC_flux_path)
rbflvrd_flux_plot <- plot_function(rbflvrd_healthy_flux_df,
                                   rbflvrd_CRC_flux_df,
                                   title_name = 'Reduced Riboflavin',
                                   tag_name = 'E')

layout <- matrix(c(1, 1, 
                   2, 3,
                   4, 5),
                 ncol = 2,
                 byrow = TRUE)
plot <- grid.arrange(community_flux_plot,
                     acetate_flux_plot,
                     alanine_flux_plot,
                     propionate_flux_plot,
                     rbflvrd_flux_plot,
                     layout_matrix = layout,
                     left = textGrob("Condition", rot = 90, gp = gpar(fontsize = 22)
                     ), 
                     bottom = textGrob("Flux (mmol/d)", gp = gpar(fontsize = 22)
                     ))
plot





