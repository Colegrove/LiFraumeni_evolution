InitialVAF_pos_plot <- MAF_table %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = VAF, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "Initial (All mutations)") + 
  theme_bw() + 
  theme()
InitialDP_pos_plot <- MAF_table %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = t_depth, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "Initial depths (All mutations)") + 
  theme_bw() + 
  theme()
InitialNF_pos_plot <- MAF_table %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = NF, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "Initial NFs (All mutations)") + 
  scale_y_log10() + 
  theme_bw() + 
  theme()
InitialNF_pos_plot 

FirstFiltersVAF_pos_plot <- filt_maf_1 %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = VAF, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "First Filters (Sample depth, per position depth, NF, masking)") + 
  theme_bw() + 
  theme()
FirstFiltersDP_pos_plot <- filt_maf_1 %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = t_depth, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "First Filters depths(Sample depth, per position depth, NF, masking)") + 
  theme_bw() + 
  theme()

FinalVAF_pos_plot <- filt_maf %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = VAF, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "Final (First Filters + SNP and contamination filters)") + 
  theme_bw() + 
  theme()
FinalDP_pos_plot <- filt_maf %>% 
  ggplot(
    aes(
      x = Start_Position, 
      y = t_depth, 
      color = coding
    )
  ) + 
  geom_point() + 
  labs(title = "Final depths(First Filters + SNP and contamination filters)") + 
  theme_bw() + 
  theme()
combined_VAF_pos_plot <- 
  cowplot::plot_grid(
    InitialVAF_pos_plot,
    InitialDP_pos_plot,
    InitialNF_pos_plot,
    FirstFiltersVAF_pos_plot,
    # FirstFiltersDP_pos_plot,
    FinalVAF_pos_plot,
    # FinalDP_pos_plot,
    ncol = 1
  )
combined_VAF_pos_plot

ggsave(paste(out_prefix(),".VAF_pos_plots.pdf"),
       plot=combined_VAF_pos_plot,
       width = 11, 
       height = 8.5,
       units = "in")
