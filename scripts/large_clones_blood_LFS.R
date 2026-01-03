
filt_maf

LFS_filt_maf <- filt_maf %>%
  mutate(LFS = if_else(Subject %in% c("Patient", "Family member A", "Family member C"), "LFS", "non-LFS")) %>%
  filter(t_alt_count >= 2 & Hugo_Symbol == "TP53") %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) #%>%
  #group_by(LFS) %>%
  #summarise(count = n())
LFS_filt_maf

large_clone_binom <- binom.test(13,14,p=0.5)

large_clone_p <- large_clone_binom$p.value

large_clones_LFS <- ggplot(LFS_filt_maf, aes(x = LFS, fill = LFS)) +
  geom_bar() + 
  #theme_classic() +
  labs(y="Number of large clones") +
  ylim(0,15) +
  annotate("text",
           x = 1.5, y = 15,
           label = deparse(bquote(italic(p) == .(round(large_clone_p, 3)))),
           parse = TRUE,
           size = 8*25.4/72.27,
           color = "black"
  ) +
  scale_fill_manual(values = c("LFS" = "#882255", "non-LFS" = "#44AA99")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size=8), 
        axis.title.y = element_text(size=8),
        legend.position = "none")
large_clones_LFS
ggsave("results/LFS_blood_large_clones.png", large_clones_LFS, width = 1.5, height = 1.5, units= "in", dpi = 300)
