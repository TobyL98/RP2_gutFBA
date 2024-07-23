#########################
# secretion_graphs.R
#########################

# Create graphs to demonstrate the biggest
# secretion differences between healthy and CRC

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# REMEMBER TO SET WORKING DIRECTORY

sec_diff_df <- read.csv("new_results/secFlux_diff.csv")


# remove generic metabolites that I don't want to analyse
sorted_sec_diff_df <- sec_diff_df %>%
  filter(metabolite != "h_medium") %>%
  filter(metabolite != "h2_medium") %>%
  filter(metabolite != "co2_medium") %>%
  arrange(desc(flux_diff))

# select the top 10 from CRC
top_10_CRC <- sorted_sec_diff_df %>%
  slice_head(n = 10) %>%
  pivot_longer(c(flux_healthy, flux_CRC))

# plotting the CRC and helathy levels of metabolites for the top 10 greatest
# differences for CRC
x_metabolite_names = c("Glycine", "Phenylalanine", "Galactose",
                       "4-Aminobutanoate", "Acetate", "Putrescine",
                       "Indole", "Ammonia", "Alanine", "Succinate")
CRC_bar <- ggplot(top_10_CRC, aes(x = reorder(metabolite, flux_diff),
                                   y = -value
                                  ,fill = name)) +
  geom_col(position = "dodge") +
  labs(y = "", x = "Metabolite", fill = "Condition") +
  scale_fill_discrete(labels = c("CRC", "Healthy")) +
  scale_x_discrete(labels = x_metabolite_names) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = (c(0.2, 0.6)))

# Plotting the top 10 percentage differences towards CRC
top_10_perc_CRC <- sorted_sec_diff_df %>%
  arrange(perc_flux_diff) %>%
  slice_head(n = 10)

x_metabolite_names2 = c("Putrescine", "Thymidine", "Phenylalanine", 
                        "4-aminobutonoate", "Galactose", "Succinate","Glycine", 
                        "Indole", "Alanine", "Acetoacetate")
CRC_diff_bar <- ggplot(top_10_perc_CRC, aes(x = reorder(metabolite, -perc_flux_diff) ,
                                       y = -perc_flux_diff)) +
  geom_col(fill = "#00BA38") +
  labs(y = "", x = "Metabolite") +
  scale_x_discrete(labels = x_metabolite_names2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(angle = 90, vjust = 0.4, 
                                    hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 18))

# Plotting the top 10 differences towards healthy
top_10_healthy <- sorted_sec_diff_df %>%
  slice_tail(n = 10) %>%
  pivot_longer(c(flux_healthy, flux_CRC))


x_metabolite_names3 <- c("Choline", "MGlcn78", "Valine", "Tn antigen",
                         "N-acgal", "Glycerol", "Butyrate",
                         "Riboflavin", "Formate", "Propionate")
healthy_bar <- ggplot(top_10_healthy, aes(x = reorder(metabolite, -flux_diff),
                                  y = -value
                                  ,fill = name)) +
  geom_col(position = "dodge") +
  labs(y = "Secretion Flux", x = "Metabolite", fill = "Condition") +
  scale_fill_manual(values = c("#999999", "#D55E00"), 
                               labels = c("CRC", "Healthy")) +
  scale_x_discrete(labels = x_metabolite_names3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = -0.5),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = (c(0.2, 0.6)))


# Plotting the top 10 percentage differences towards healthy
# rmeoved MGlcn78 as it increased by 200% and was distorting the graph
# will need to include in figure
top_10_perc_healthy <- sorted_sec_diff_df %>%
  arrange(perc_flux_diff) %>%
  na.omit() %>%
  slice_tail(n = 10) %>%
  filter(metabolite != "MGlcn78_rl_medium")

x_metabolite_names4 <- c("Formate", "Butyrate", "Propionate", "Riboflavin",
                         "Valine", "Glycerol", "Tn antigen", "N-acgal",
                         "Choline")
healthy_diff_bar <- ggplot(top_10_perc_healthy, aes(x = reorder(metabolite, perc_flux_diff) ,
                                            y = perc_flux_diff)) +
  geom_col(fill = "#CC79A7") +
  labs(y = "Secretion Flux Difference (%)", x = "Metabolite") +
  scale_x_discrete(labels = x_metabolite_names4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(angle = 90, vjust = 0.5, 
                                    hjust = -0.5),
        axis.title = element_text(size = 18))

       
layout <- matrix(c(1, 3, 
                   1, 3,
                   1, 3,
                   2, 4,
                   2, 4), nrow = 2)
plot <- grid.arrange(CRC_bar, CRC_diff_bar, healthy_bar, healthy_diff_bar,
             layout_matrix = layout)

