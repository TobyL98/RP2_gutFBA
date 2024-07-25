##################
# species_graphs.R
#################

# Creates a grpah that shows biggest difference between
# species biomass flux in healthy and CRC samples

library(dplyr)
library(ggplot2)
library(patchwork)

spec_diff_df <- read.csv("final_results/spec_flux_diff.csv")

# sorting the dataframe
spec_diff_df <- spec_diff_df %>%
  arrange(percentage_flux_diff) %>%
  rename(Species = Species_biomass) %>%
  mutate(str_split(Species, "_", n = 1)) # fix this

# pick top 5 differences towards CRC
top5_CRC <- spec_diff_df %>%
  slice_tail(n = 5) %>%
  pivot_longer(c(percentage_flux_healthy, percentage_flux_CRC))

# plot top 5 CRC differences
top5_CRC_bar <- ggplot(top5_CRC, aes(x = Species_biomass,
                                     y = value,
                                     fill = name)) +
  geom_col(position = "dodge")

top5_CRC_bar