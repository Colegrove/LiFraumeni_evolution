## Lollipop plot

# default plot size
source("scripts/01.1.05_makeLolipopComparison.R")
# plot for slide
source("scripts/makeLolipopComparison_slide.R")

# TP53 domains
filteredDomains <- read_delim(inputs$domains_file)

#tissue = "All MSTP volunteers"
#tissue = "Skin"
mstp <- c("All MSTP volunteers", "UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")

tissue_filtering <- filt_maf %>% 
  filter(!is.na(prot.pos)) %>% 
  filter(Subject %in% mstp) %>%
  filter(coding == "coding") %>% 
  filter(inMask == FALSE) %>% 
  filter(Variant_Type=="SNP")

## uncomment to plot all participants individually

# all_tissues_plot <- function(tissue_df){
#   lollipop_src <- lolipopComparison(
#     #tissue_filtering, 
#     tissue_df,
#     tissue,
#     COSMIC_table_forLolipop, 
#     filteredDomains, 
#     data_name = "Sample", 
#     comp_data_name = "COSMIC")
#   
#   # if tissue specific
#   lollipop_plot <- lollipop_src[["TP53"]]
#   #show(lollipop_plot[["aligned"]])
#   # if not tissue specific
#   lollipop_plot_peices <- lollipop_src[["TP53"]]
#   
#   slide_theme <- theme(
#     axis.text.x = element_text(size=10),
#     axis.text.y = element_text(size = 10),
#     axis.title.y = element_text(size=10)
#   )
# 
#   Gene <- "TP53"
#   plot_title <- ggdraw() +
#     draw_label(
#       glue("{Gene} - {volunteer}"),
#       #GeneIter,
#       fontface = 'bold.italic',
#       x = 0,
#       hjust = 0,
#       size = 12
#     ) +
#     theme(
#       # add margin on the left of the drawing canvas,
#       # so title is aligned with left edge of first plot
#       plot.margin = margin(0, 0, 0, 7),
#     )
#   
#   test_plot <- lollipop_plot_peices[["test_plot"]] + slide_theme
#   test_domains_plot <- lollipop_plot_peices[["test_domains_plot"]]
#   test_plot_2 <- lollipop_plot_peices[["test_plot_2"]] + slide_theme
#   domain_legend <- lollipop_plot_peices[["domain_legend"]] + slide_theme
#   
#   
#   domain_legend <- get_legend(
#     test_domains_plot +
#       theme(legend.box.margin = margin(7, 7, 7, 7),
#             legend.direction = "horizontal",
#             legend.background = element_blank(), 
#             legend.text = element_text(size = 10)) +
#       guides(
#         fill=guide_legend(
#           nrow=min(3, length(filteredDomains$Label))
#         )
#       )
#   )
#   
#   lollipop_slide <- cowplot::plot_grid(
#     plot_title,
#     test_plot,
#     test_domains_plot + theme(legend.position="none"),
#     test_plot_2,
#     domain_legend,
#     align='v', axis="lr", rel_heights = c(1,4,1,3,3), ncol = 1)
#   
#   show(lollipop_slide)
#   ggsave(glue("results/lollipop_plot_slide_{volunteer}.png"), lollipop_slide, width = 6, height = 3, units = "in", dpi = 300)
# }
# 
# 
# for(volunteer in mstp){
#   if(volunteer == "All MSTP volunteers"){
#     all_tissues_plot(tissue_filtering)
#   }
#   else{
#     tissue_df <- tissue_filtering %>% 
#       filter(Subject == volunteer)
#     all_tissues_plot(tissue_df)
#   }
# }




print(unique(tissue_filtering$prot.pos))
      
      
family_title = "All MSTP volunteers"
lollipop_src <- lolipopComparison(
  tissue_filtering, 
  family_title,
  COSMIC_table_forLolipop, 
  filteredDomains, 
  data_name = "Sample", 
  comp_data_name = "COSMIC")

# if tissue specific
lollipop_plot <- lollipop_src[["TP53"]]
#show(lollipop_plot[["aligned"]])
# if not tissue specific
lollipop_plot_peices <- lollipop_src[["TP53"]]

slide_theme <- theme(
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size=10)
)

Gene <- "TP53"
plot_title <- ggdraw() +
  draw_label(
    glue("{Gene} - {family_title}"),
    #GeneIter,
    fontface = 'bold.italic',
    x = 0,
    hjust = 0,
    size = 12
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7),
  )

test_plot <- lollipop_plot_peices[["test_plot"]] + slide_theme
test_domains_plot <- lollipop_plot_peices[["test_domains_plot"]]
test_plot_2 <- lollipop_plot_peices[["test_plot_2"]] + slide_theme
domain_legend <- lollipop_plot_peices[["domain_legend"]] + slide_theme


lollipop_slide <- cowplot::plot_grid(
  plot_title,
  test_plot,
  test_domains_plot + theme(legend.position="none"),
  test_plot_2,
  domain_legend,
  align='v', axis="lr", rel_heights = c(1,4,1,3,3), ncol = 1)


show(lollipop_slide)
ggsave("results/lollipop_plot_mstp_controls.png", lollipop_slide, width = 6, height = 3, units = "in", dpi = 300)
