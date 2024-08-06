# Variation.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)

setwd('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research Project 2/development/RP2_gutFBA/first_model')
getwd()

abundances_df <- read.csv('Outputs/healthy_df_out.csv', header = T)


# picking out just bacteroides
bacteroides_abundances_df <- abundances_df %>%
  filter(Genus == "Bacteroides") %>%
  pivot_longer(!c(Genus, species),
               names_to = "Patient",
               values_to = "Abundance")

# bacteroides summary statistics
bacteroides_abundance_stats <- bacteroides_abundances_df %>%
  summarise(mean_abundance = mean(Abundance),
            SD_abundance = sd(Abundance),
            n = n())

# bacteroides plot
bacteroides_abundances_plot <- ggplot(bacteroides_abundances_df,
                                      aes(x = Abundance)) +
  geom_histogram(fill = "#1B9E77",
                 color = "white",
                 binwidth = 5,
                 boundary = 0) +
  labs(title = "Bacteroides") +
  geom_vline(aes(xintercept = mean_abundance), bacteroides_abundance_stats,
             color = "green", linewidth = 2) +
  #scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        title = element_text(size = 16))


# picking out just prevotella
prevo_abundances_df <- abundances_df %>%
  filter(Genus == "Prevotella") %>%
  pivot_longer(!c(Genus, species),
               names_to = "Patient",
               values_to = "Abundance")

# prevotella_summary_statistics
prevo_abundance_stats <- prevo_abundances_df %>%
  summarise(mean_abundance = mean(Abundance),
            SD_abundance = sd(Abundance),
            n = n())

# prevotella plot
prevo_abundances_plot <- ggplot(prevo_abundances_df,
                                      aes(x = Abundance)) +
  labs(title = "Prevotella") +
  geom_histogram(fill = "#D95F02",
                 color = "white",
                 binwidth = 5,
                 boundary = 0) +
  geom_vline(aes(xintercept = mean_abundance), prevo_abundance_stats,
             color = "green", linewidth = 2) +
  #scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        title = element_text(size = 16))


# picking out just eubacterium
eubac_abundances_df <- abundances_df %>%
  filter(Genus == "Eubacterium") %>%
  pivot_longer(!c(Genus, species),
               names_to = "Patient",
               values_to = "Abundance")

# prevotella_summary_statistics
eubac_abundance_stats <- eubac_abundances_df %>%
  summarise(mean_abundance = mean(Abundance),
            SD_abundance = sd(Abundance),
            n = n())

# prevotella plot
eubac_abundances_plot <- ggplot(eubac_abundances_df,
                                aes(x = Abundance)) +
  geom_histogram(fill = "#7570B3",
                 color = "white",
                 binwidth = 5,
                 boundary = 0) +
  geom_vline(aes(xintercept = mean_abundance), eubac_abundance_stats,
             color = "green", linewidth = 2) +
  labs(title = "Eubacterium") +
  #scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        title = element_text(size = 16))

plot <- grid.arrange(bacteroides_abundances_plot,
                     prevo_abundances_plot,
                     eubac_abundances_plot,
                     left = textGrob("Count", rot = 90, gp = gpar(fontsize = 20)),
                     bottom = textGrob("Abundance (%)", gp = gpar(fontsize = 20)))

plot

