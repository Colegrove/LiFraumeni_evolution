

########################## tp53

am_groups <- am_groups %>% filter(gene_name %in% c("TP53"))

unique(am_groups$lfs_group)
wilcox.test(am_pathogenicity.y ~ lfs_group, 
            data = am_groups %>% filter(lfs_group %in% c("LFS\nno-CTX", "non-LFS\nCTX")), exact = FALSE)

obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = lfs_group, y = am_pathogenicity.y)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, stroke = 1.2, alpha = 1, aes(color = LFS, shape = CTX)) +
  scale_color_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +
  scale_shape_manual(values = c("CTX" = 16, "no-CTX" = 1)) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups

ggsave("results/AlphaMissense_tp53_supp.png", obsv_am_blood_groups, width = 3, height = 2, units = "in", dpi = 300)

