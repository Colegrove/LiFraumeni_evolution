## Lolipop plot comparison with domains
require("tidyverse")
require("cowplot")
require("glue")




plot_lollipop_single_gene <- function(
    my_counts,
    cosmic_counts,
    domains,
    gene_name = "TP53",
    tissue = "LFS esophagus",
    data_name = "Sample",
    cosmic_name = "COSMIC",
    color_scale = c("#56B4DF", "#E69D00", "#009E74", "#D55C00", "#0071B2")
) {
  
  prot_len <- max(domains$aa.length)
  
  p_top <- ggplot(my_counts) +
    geom_linerange(aes(x = prot.pos, ymin = 0, ymax = Count), size = 0.25) +
    geom_point(aes(x = prot.pos, y = Count), size = 1) +
    scale_x_continuous(limits = c(0, prot_len)) +
    scale_y_continuous(name = paste0(data_name, "\nCount")) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      plot.margin = margin(0, 7, 0, 7)
    )
  
  p_domains <- ggplot() +
    geom_rect(
      aes(xmin = 0, xmax = prot_len, ymin = -0.25, ymax = 0.25),
      fill = "grey90"
    ) +
    geom_rect(
      data = domains,
      aes(xmin = Start, xmax = End, ymin = -0.35, ymax = 0.35, fill = Label2)
    ) +
    scale_fill_manual(values = color_scale) +
    scale_x_continuous(limits = c(0, prot_len)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      plot.margin = margin(0, 7, 0, 7)
    )
  
  domain_legend <- get_legend(p_domains)
  
  p_bottom <- ggplot(cosmic_counts) +
    geom_linerange(aes(x = prot.pos, ymin = -Count, ymax = 0), size = 0.25) +
    scale_x_continuous(limits = c(0, prot_len)) +
    scale_y_continuous(labels = abs, name = paste0(cosmic_name, "\nCount")) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(0, 7, 0, 7)
    )
  
  p_title <- ggdraw() +
    draw_label(
      glue("{gene_name} - {tissue}"),
      fontface = 'bold.italic',
      x = 0,
      hjust = 0,
      size = 10
    ) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  cowplot::plot_grid(
    p_title,
    p_top,
    p_domains + theme(legend.position = "none"),
    p_bottom,
    domain_legend,
    ncol = 1,
    rel_heights = c(1, 4, 1, 3, 2),
    align = 'v',
    axis = "lr"
  )
}

plot_lollipop_single_gene(my_counts, cosmic_counts, domains)