##################
# MES_graphs.R
##################

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

MES_diff_df <- read.csv("final_results/MES_final_results.csv")

MES_diff_df <- MES_diff_df %>%
  drop_na() %>%
  filter(str_detect(metabolite_id, '_rl_', negate = TRUE))

top10_MES_CRC <-  MES_diff_df %>%
  filter(metabolite_id != "o2_medium") %>%
  slice_head(n = 10)

metabolite_names <- c("Cholate", "Glutathione", "Adenosine", "Glycerol",
                      "Butyrate", "Urea", "Threonine", "Asparagine",
                      "Tryptophan", "Arabinose")
top10_MESCRC_bar <- ggplot(top10_MES_CRC, 
                        aes(reorder(x = metabolite_name_x, -MES_diff),
                                          y = -MES_diff)) +
  geom_col(fill = 'maroon') +
  labs(x = "Metabolite", 
       y = "MES Difference",
       title = "CRC") +
  scale_x_discrete(labels = metabolite_names) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
    axis.text = element_text(size = 16),
    axis.title = element_blank(),
    title = element_text(size = 18)
  )

top10_MESCRC_bar

top10_MES_healthy <- MES_diff_df %>%
  slice_tail(n = 10)

metabolite_names1 <- c("NANA", "MGlcn191", "MGlcn186", "GlcNAc", "Citrate", 
                       "Methionine", "MGlcn153", "Galactosamine", "2,6-DAP", 
                       "Deoxyguanosine")
top10_MEShealthy_bar <- ggplot(top10_MES_healthy, 
                               aes(reorder(x = metabolite_name_x, MES_diff),
                                   y = MES_diff)) +
  geom_col(fill = 'blue') +
  labs(x = "Metabolite", 
       y = "MES Difference",
       title = "Healthy") +
  scale_x_discrete(labels = metabolite_names1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        axis.text = element_text(size = 16),
        axis.title = element_blank(),
        title = element_text(size = 18)
  )
top10_MEShealthy_bar

plot <- grid.arrange(top10_MESCRC_bar, top10_MEShealthy_bar, 
                     left = textGrob("MES Difference", rot = 90, gp = gpar(fontsize = 20)), 
                     bottom = textGrob("Metabolite", gp = gpar(fontsize = 20)))
plot


###############
# density plot
##############

# plot the density of MES scores for
# CRC vs healthy
compare_MES_df <- MES_diff_df %>%
  pivot_longer(c(MES_health, MES_CRC))

density_plot <- ggplot(compare_MES_df, aes(x = value
                                           , color = name
                                           , fill = name)) +
  labs(x = "Metabolic Exchange Score (MES)",
       y = "Density",
       fill = "Condition",
       color = "Condition") +
  geom_density(alpha = 0.3) +
  scale_fill_discrete(labels = c("CRC", "Healthy")) +
  scale_color_discrete(labels = c("CRC", "Healthy")) +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = (c(0.8, 0.6)))

density_plot
  