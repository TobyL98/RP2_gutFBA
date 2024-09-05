# Variation.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)

#########
# Functions
#########

filter_function <- function(df, genus_name) {
  # correctly filters the data
  # to obtain just one bacterias abundances
  df %>%
  filter(Genus == genus_name) %>%
  pivot_longer(!c(Genus, species),
               names_to = "Patient",
               values_to = "Abundance")
}

plot_function <- function(df, colour, genus_name, tag_name, mean_value) {
  # Plots the histogram showing the variation
  # distribution
  ggplot(df,
         aes(x = Abundance)) +
    geom_histogram(fill = colour,
                   color = "white",
                   binwidth = 5,
                   boundary = 0) +
    labs(title = genus_name,
         tag = tag_name) +
    geom_vline(aes(xintercept = mean_value), bacteroides_abundance_stats,
               color = "green", linewidth = 2) +
    #scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_blank(),
          title = element_text(size = 16))
}

stats_function <- function(abundance_df) {
  # calculates the mean and SD for abundances of
  # a specific bacteria
  abundance_df %>%
  summarise(mean_abundance = mean(Abundance),
            SD_abundance = sd(Abundance),
            n = n())
}

setwd('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research Project 2/development/RP2_gutFBA/first_model')
getwd()

abundances_df <- read.csv('Outputs/healthy_df_out.csv', header = T)


# picking out just bacteroides
bacteroides_abundances_df <- filter_function(abundances_df,
                                             genus_name = 'Bacteroides')

# bacteroides summary statistics
bacteroides_abundance_stats <- stats_function(bacteroides_abundances_df)
# bacteroides plot
bacteroides_abundances_plot <- plot_function(bacteroides_abundances_df,
                                             colour = '#1B9E77',
                                             genus_name = 'Bacteroides - Healthy',
                                             tag_name = 'A',
                                             bacteroides_abundance_stats$mean_abundance)

# picking out just prevotella
prevo_abundances_df <- filter_function(abundances_df,
                                       genus_name = 'Prevotella')

# prevotella_summary_statistics
prevo_abundance_stats <- stats_function(prevo_abundances_df)

# prevotella plot
prevo_abundances_plot <- plot_function(prevo_abundances_df,
                                       colour = '#1B9E77',
                                       genus_name = 'Prevotella - Healthy',
                                       tag_name = 'B',
                                       prevo_abundance_stats$mean_abundance)


# picking out just eubacterium
eubac_abundances_df <- filter_function(abundances_df,
                                       genus_name = 'Eubacterium')

# prevotella_summary_statistics
eubac_abundance_stats <- stats_function(eubac_abundances_df)
# prevotella plot
eubac_abundances_plot <- plot_function(eubac_abundances_df,
                                       colour = '#1B9E77',
                                       genus_name = 'Eubacterium - Healthy',
                                       tag_name = 'C',
                                       eubac_abundance_stats$mean_abundance)

# read in data for CRC
crc_abundances_df <- read.csv('Outputs/Stage_I_II_df_out.csv', header = T)

crc_bacteroides_df <- filter_function(crc_abundances_df,
                                      genus_name = 'Bacteroides')
crc_bacteroides_stats <- stats_function(crc_bacteroides_df)
crc_bacteroides_plot <- plot_function(crc_bacteroides_df,
                                      colour = '#D95F02',
                                      genus_name = 'Bacteroides - CRC',
                                      tag_name = 'D',
                                      crc_bacteroides_stats$mean_abundance)

crc_prevo_df <- filter_function(crc_abundances_df,
                                genus_name = 'Prevotella')
crc_prevo_stats <- stats_function(crc_prevo_df)
crc_prevo_plot <- plot_function(crc_prevo_df,
                                colour = '#D95F02',
                                genus_name = 'Prevotella - CRC',
                                tag_name = 'E',
                                crc_prevo_stats$mean_abundance)

crc_eubac_df <- filter_function(crc_abundances_df,
                                genus_name = 'Eubacterium')
crc_eubac_stats <- stats_function(crc_eubac_df)
crc_eubac_plot <- plot_function(crc_eubac_df,
                                  colour = '#D95F02',
                                  genus_name = 'Eubacterium - CRC',
                                tag_name = 'F',
                                crc_eubac_stats$mean_abundance)
crc_eubac_plot

layout <- matrix(c(1, 4,
                          2, 5,
                          3, 6), 
                        nrow = 3,
                        byrow = TRUE)
plot <- grid.arrange(bacteroides_abundances_plot,
                     prevo_abundances_plot,
                     eubac_abundances_plot,
                     crc_bacteroides_plot,
                     crc_prevo_plot,
                     crc_eubac_plot,
                     layout_matrix = layout,
                     left = textGrob("Count", rot = 90, gp = gpar(fontsize = 20)),
                     bottom = textGrob("Abundance (%)", gp = gpar(fontsize = 20)))

plot
