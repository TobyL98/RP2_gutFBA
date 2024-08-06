###################
# species_graph2
###################

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
  mutate(Species = str_split_i(Species, "_", 1)) %>%
  mutate(Health_or_CRC = percentage_flux_diff >= 0)

spec_diff_bar <- ggplot(spec_diff_df, 
                        aes(x = percentage_flux_diff,
                            y = reorder(Species, percentage_flux_diff),
                            fill = Health_or_CRC)) +
  geom_col() +
  labs(y = "Species",
       x = "Percentage Abundance Difference",
       fill = "Condition") +
  scale_fill_discrete(labels = c("CRC", "Healthy")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = (c(0.7, 0.4)))
spec_diff_bar
