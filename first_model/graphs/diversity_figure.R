################
# diversity_figure.R
################

# plots the shannon index
# and richness for healthy and CRC

library(dplyr)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)

### FUNCTIONS
shannon_boxplot_function <- function(df, title_name, colours, tag_name) {
  # plots the boxplots for both richness and shannon_index
  ggplot(df, 
         aes(x= source, 
             y = shannon_index,
             fill = source)) +
    labs(title = title_name,
         tag = tag_name) +
    scale_fill_manual(values = colours) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_blank(),
          legend.position = "none")
  
}

richness_boxplot_function <- function(df, title_name, colours, tag_name) {
  # plots the boxplots for both richness and shannon_index
  ggplot(df, 
         aes(x= source,
             y = richness,
             fill = source)) +
    labs(title = title_name,
         tag = tag_name) +
    scale_fill_manual(values = colours) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_blank(),
          legend.position = "none")
  
}

### MAIN CODE
genus_CRC_df <- read_csv('diversity/genus_CRC.csv') %>%
  mutate(source = 'CRC')
genus_healthy_df <- read_csv('diversity/species_healthy.csv') %>%
  mutate(source = 'Healthy')
species_CRC_df <- read_csv('diversity/species_CRC.csv') %>%
  mutate(source = 'CRC')
species_healthy_df <- read_csv('diversity/species_healthy.csv') %>%
  mutate(source = 'Healthy')

# to check if the data is normally distributed, it isn't
histo_plot <- ggplot(species_CRC_df, 
                     aes(x= shannon_index)) +
  geom_histogram()
histo_plot

qq_plot <- ggplot(species_CRC_df,
                  aes(sample = shannon_index)) +
  stat_qq() +
  stat_qq_line()
qq_plot

# merging datasets for plotting
combined_genus_df <- bind_rows(genus_CRC_df, genus_healthy_df)
combined_species_df <- bind_rows(species_CRC_df, species_healthy_df)

colour1 <- c("#ea5545", "#87bc45")
colour2 <- c("#ef9b20", "#27aeef")
# plotting box plots
genus_shannon_box <- shannon_boxplot_function(combined_genus_df,
                                              title = 'Genus - Shannon Index',
                                              colour1,
                                              tag_name = 'C')
species_shannon_box <- shannon_boxplot_function(combined_species_df,
                                                title = 'Species - Shannon Index',
                                                colour2,
                                                tag_name = 'D')
genus_richness_box <- richness_boxplot_function(combined_genus_df,
                                                title = 'Genus - Richness',
                                                colour1,
                                                tag_name = 'E')
species_richness_box <- richness_boxplot_function(combined_species_df,
                                                  title = 'Species - Richness',
                                                  colour2,
                                                  tag_name = 'F')

# plotting shannon and richness data for 95% average
healthy_95_df <- read_csv('diversity/healthy95.csv') %>%
  mutate(source = 'Healthy')
crc_95_df <- read_csv('diversity/crc95.csv') %>%
  mutate(source = 'CRC')
combined_95_df <- bind_rows(healthy_95_df, crc_95_df)

bar_95_shannon <- ggplot(combined_95_df, 
                         aes(x = source, y = shannon_index, fill = source)) +
  labs(title = 'Genus 95% - Shannon Index',
       tag = 'A') +
  geom_col() +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_blank(),
        legend.position = "none")
  
bar_95_richness <- ggplot(combined_95_df, 
                         aes(x = source, y = richness, fill = source)) +
  labs(title = 'Genus 95% - Richness',
       tag = 'B') +
  geom_col() +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_blank(),
        legend.position = "none")
bar_95_richness


plot <- grid.arrange(bar_95_shannon,
                     bar_95_richness,
                     genus_shannon_box, 
                     genus_richness_box,
                     species_shannon_box,
                     species_richness_box,
                     nrow = 3,
                     bottom = textGrob("Condition", gp = gpar(fontsize = 22)
                     ))
plot

# doing t-tests to see if there is a difference
var.test(genus_CRC_df$shannon_index, genus_healthy_df$shannon_index)
# genus shannon index, is a significant difference
t.test(genus_CRC_df$shannon_index, genus_healthy_df$shannon_index)

# genus richness, is a significant difference
t.test(genus_CRC_df$richness, genus_healthy_df$richness)

# species shannon index, no significant difference
t.test(species_CRC_df$shannon_index, species_healthy_df$shannon_index)

# species richness, is a significant difference
t.test(species_CRC_df$richness, species_healthy_df$richness)