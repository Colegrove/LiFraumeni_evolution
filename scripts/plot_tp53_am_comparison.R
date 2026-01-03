### Comparison of AlphaMissense pathogenicity between obverved and all possible missense variants

library(ggplot2)
library(dplyr)


######
## plot all observed together
######
obsv_am_blood <- ggplot(combine_df, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c("Not observed" = "grey60", "Observed" = "#3182bd")) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )

######
## separate large and small observed clones
######
combine_df_large_clones <- combine_df %>%
  mutate(group = case_when(
    group == "Observed" & t_alt_count == 1 ~ "Observed (1 read)",
    group == "Observed" & t_alt_count > 1  ~ "Observed (>1 read)",
    TRUE ~ group
  ))
combine_df_large_clones$group <- factor(combine_df_large_clones$group,
                                        levels = c("Not observed", "Observed (1 read)", "Observed (>1 read)"))

obsv_am_large_clones_blood <- ggplot(combine_df_large_clones, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c(
    "No nobserved" = "grey60",
    "Observed (1 read)" = "#6baed6", 
    "Observed (>1 read)" = "#08519c" 
  )) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) 

######
#### save plots
######

ggsave("results/am_blood.png", obsv_am_blood, width = 2.5, height = 1.5, units = "in", dpi = 300)
ggsave("results/am_large_clones_blood.png", obsv_am_large_clones_blood, width = 2.5, height = 1.5, units = "in", dpi = 300)
