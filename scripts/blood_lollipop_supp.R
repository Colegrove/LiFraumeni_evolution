## Lollipop plot

# default plot size
source("scripts/01.1.05_makeLolipopComparison.R")
# plot for slide
source("scripts/makeLolipopComparison_slide.R")

# TP53 domains
filteredDomains <- read_delim(inputs$domains_file)
filteredDomains
subjects_blood <- c("UW volunteer 1","UW volunteer 2","UW volunteer 3","UW volunteer 4",
                                 "UW volunteer 5","UW volunteer 6","UW volunteer 7",
                                 "Patient","Family member A","Family member C","Family member B")
subjects_LFS <- c("Patient","Family member A","Family member C")


## input cosmic data
cosmic_data <- read_delim("inputs/COSMIC/TP53_mutations_07OCT25_23_08_31.csv")
cosmic_counts <- cosmic_data %>% 
  group_by(Position) %>% 
  summarise(count = sum(Count))


blood_filtering <- maf_masked_coding %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP")


#### Slide size

tissue <- "All tissues"
lollipop_src <- lolipopComparison_slide(
  blood_filtering, 
  tissue,
  COSMIC_table_forLolipop, 
  filteredDomains, 
  data_name = "Sample", 
  comp_data_name = "COSMIC")

lollipop_src
# if not tissue specific
lollipop_plot_peices <- lollipop_src[["TP53"]]

ms_theme <- theme(
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size = 8),
  axis.title.y = element_text(size=8)
)


cosmic_counts <- cosmic_data %>%
  group_by(Position) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  dplyr::rename(prot.pos = Position)

test_plot <- lollipop_plot_peices[["test_plot"]] + ms_theme + theme(
  plot.margin = margin(2,2,2,2)
)
test_domains_plot <- lollipop_plot_peices[["test_domains_plot"]]
#test_plot_2 <- lollipop_plot_peices[["test_plot_2"]] + ms_theme

test_plot_2 <- ggplot(cosmic_counts) +
  geom_linerange(aes(x = prot.pos, ymin = -Count, ymax = 0), size = 0.25) +
  scale_x_continuous(limits = c(0, prot_len)) +
  scale_y_continuous(labels = abs, name = "COSMIC\nCount") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(2, 2, -2, 2), 
    axis.title.y = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8)
  )

domain_legend <- get_legend(test_domains_plot +
                              theme(legend.key.size = unit(8, "pt"),
                                    legend.text = element_text(size = 8, margin= margin(0,0,0,1)), 
                                    legend.box.margin = margin(0, 0, 0, 0),
                                    legend.direction = "horizontal",
                                    legend.background = element_blank(),
                                    legend.margin = margin(0,0,0,0),
                                    legend.key.spacing.x = unit(3, "pt"))
)

#domain_legend <- get_legend(test_domains_plot)

lollipop_ms <- cowplot::plot_grid(
  test_plot,
  test_domains_plot + theme(legend.position="none"),
  test_plot_2,
  domain_legend,
  align='v', axis="lr", rel_heights = c(3.25,1,2,1.25), ncol = 1) +
  theme(plot.margin = margin(b=-5))

show(lollipop_ms)
ggsave("results/lollipop_blood_supp.png", lollipop_ms, width = 3, height = 1.5, units = "in", dpi = 300)
