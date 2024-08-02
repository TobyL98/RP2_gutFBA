##################
# species_graphs.R
#################

# Creates a grpah that shows biggest difference between
# species biomass flux in healthy and CRC samples

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

spec_diff_df <- read.csv("final_results/spec_flux_diff.csv")

# sorting the dataframe
spec_diff_df <- spec_diff_df %>%
  arrange(percentage_flux_diff) %>%
  rename(Species = Species_biomass) %>%
  mutate(Species = str_split_i(Species, "_", 1)) # fix this

# pick top 5 differences towards healthy
top5_healthy <- spec_diff_df %>%
  slice_tail(n = 5) %>%
  pivot_longer(c(percentage_flux_healthy, percentage_flux_CRC))

# plot top 5 healthy differences
top5_healthy_bar <- ggplot(top5_healthy, aes(x = reorder(Species, percentage_flux_diff),
                                     y = value,
                                     fill = name)) +
  geom_col(position = "dodge") +
  labs(y = "Percentage of Community Biomass Flux", 
       x = "Species",
       fill = "Condition") +
  scale_fill_discrete(labels = c("CRC", "Healthy")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = , vjust = 1, hjust=1),
        axis.text = element_text(size = 16),
        axis.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = (c(0.2, 0.7)))


# pick top 5 differences towards CRC
top5_CRC <- spec_diff_df %>%
  slice_head(n = 5) %>%
  pivot_longer(c(percentage_flux_healthy, percentage_flux_CRC))

# plot top 5 CRC
top5_CRC_bar <- ggplot(top5_CRC, aes(x = Species,
                                     y = value,
                                     fill = name)) +
  geom_col(position = "dodge") +
  labs(y = "Percentage of Biomass Flux", 
       x = "Species",
       fill = "Condition") +
  scale_fill_manual(values = c("#999999", "#D55E00"), 
                    labels = c("CRC", "Healthy")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = , vjust = 1, hjust=1),
    axis.text = element_text(size = 16),
    axis.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = (c(0.2, 0.7)))




plot <- grid.arrange(top5_CRC_bar, top5_healthy_bar,
                     left = textGrob("Percentage of Community Biomass", rot = 90, gp = gpar(fontsize = 20)
                     ), 
                     bottom = textGrob("Species", gp = gpar(fontsize = 20)
                     )
                     )