## Lolipop plot comparison with domains
require("tidyverse")
require("cowplot")
require("glue")


#' Make a lolipop plot comparing your data to data from a known mutations database, by gene.
#'
#' `lolipopComparison` returns a list of plots comparing data to known mutations, and with domains shown.
#'
#' @param in_data
#'   A tibble of coding mutations from a preprocessed MAF-based file.
#'   Must include column names:
#'     "Hugo_Symbol"
#'     "prot.pos"
#'     "VAF"
#'     "coding"
#'
#' @param comp_data:
#'   A file produced from a database of mutations.
#'   Should include a count of number of mutations per codon per gene.
#'   Must include column names:
#'     "Gene"
#'     "prot.pos"
#'      "Count"
#' @param in_domain_data:
#'   A tibble of domain annotations, including the following:
#'     "HGNC"
#'     "refseq.ID"
#'     "protein.ID"
#'     "aa.length"
#'     "Start"
#'     "End"
#'     "domain.source"
#'     "Label"
#'     "Label2"
#'     "domain.anno"
#'     "pfam"
#'     "Description"
#'   Label2 will be used to label the domains in the plot.
#' @param in_hotspot_data:
#'   The hotspots that were actually sequenced, if not whole-gene sequencing.
#'   column headers of:
#'     "Gene"
#'     "Start"
#'     "Stop"
#' @param data_name The label to give this data on the Y axis
#' @param comp_data_name The name of the database your comparison data came from.
#' @param comp_data_name The name of the database your comparison data came from.
 

lolipopComparison <- function(
  in_data, tissue, comp_data,
  in_domain_data, in_hotspot_data=tibble(),
  data_name = "Sample 1",
  comp_data_name = "Database",
  font_sizes = list("small" = 5, "medium" = 6, "large" = 7), 
  color_scale = c(
    "#000000",
    "#E69D00",
    "#56B4DF",
    "#009E74",
    "#F0E442",
    "#0071B2",
    "#D55C00",
    "#CC79A7"
  )
) {
  forPlot <- in_data %>%
    filter(!is.na(prot.pos)) %>%
    select(
      Gene=Hugo_Symbol,
      prot.pos, VAF, coding
    ) %>%
    group_by(Gene, prot.pos) %>%
    summarize(Count=n())

  domain_plot_list = list()

  for (GeneIter in unique(forPlot$Gene)) {
    plot_data <- forPlot %>%
      filter(Gene == GeneIter) %>% filter()
    max(plot_data$prot.pos[!is.na(plot_data$prot.pos)])
    plot_comp <- comp_data %>%
      filter(Gene == GeneIter)



    plot_domains <- in_domain_data %>%
      filter(HGNC == GeneIter)

    plot_genes <- plot_domains %>%
      select(HGNC, aa.length) %>%
      unique()

    plot_hotspots <- in_hotspot_data

    prot_len = max(plot_genes$aa.length)
    
    test_plot <- ggplot() +
      geom_linerange(
        data=plot_data,
        mapping=aes_string(
          x = "prot.pos", ymin=0,ymax="Count"
        ), size=0.25) +
      geom_point(
        data=plot_data,
        mapping=aes_string(x="prot.pos",y="Count"),
        size=.25) +
      expand_limits(y=c(3)) +
      scale_x_continuous(limits=c(0, prot_len))  +
      scale_y_continuous(
        name = paste0(data_name," Mutation\nCount")) +
      theme_bw() +
      theme(
        axis.text.x = element_text(family = "sans",
                                   size=font_sizes[["small"]],
                                   color = "black"),
        axis.text.y = element_text(family = "sans",
                                   size=font_sizes[["small"]],
                                   color = "black"),
        axis.title.y = element_text(family = "sans", 
                                    size = font_sizes[["medium"]], 
                                    color = "black"),
        axis.title.x = element_blank(),
        panel.border = element_rect(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.grid = element_line(size = 0.25),
        plot.margin = margin(0, 7, 0, 7)
      )
    test_domains_plot <- ggplot() +
      geom_rect(
        data = plot_genes,
        mapping = aes(xmin = 0,
                      xmax = aa.length,
                      ymin = -0.25,
                      ymax = 0.25)
        ) +
      geom_rect(
        data=plot_domains,
        mapping=aes(xmin=Start,
                    xmax=End,
                    ymin=-.35,
                    ymax=.35,
                    fill = Label2)
        ) +
      scale_fill_manual(values=color_scale)
    if (length(rownames(plot_hotspots))) {
      test_domains_plot = test_domains_plot +
        geom_rect(data=plot_hotspots,
                  mapping=aes(xmin=Start,
                              xmax=Stop,
                              ymin=-.75,
                              ymax=-.35),
                  fill="black", alpha=0.5)
    } else {
      test_domains_plot = test_domains_plot +
        geom_rect(data=plot_genes,
                  mapping=aes(xmin = 0,
                              xmax = aa.length,
                              ymin=-.85,
                              ymax=-.45),
                  fill="black", alpha=0.5)
    }

    test_domains_plot = test_domains_plot +
      scale_x_continuous(limits=c(0, prot_len)) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.text = element_text(family = "sans",
                                       size=font_sizes[["small"]],
                                       color = "black"),
            legend.key.size = unit(font_sizes[["medium"]] + 1, "pt"),
            legend.title = element_blank(),
            plot.margin = margin(0, 7, 0, 7)
      )
    test_plot_2 <- ggplot() +
      geom_linerange(
        data = plot_comp,
        mapping = aes(x=prot.pos, ymin=-Count, ymax = 0),
        size = 0.25) +
      scale_x_continuous(limits=c(0, prot_len)) +
      scale_y_continuous(
        labels=abs,
        name = paste0(comp_data_name, "\nCount")) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "sans",
                                   size=font_sizes[["small"]],
                                   color = "black"),
        axis.title.y = element_text(family = "sans", 
                                    size = font_sizes[["medium"]], 
                                    color = "black"),
        axis.title.x = element_blank(),
        panel.border = element_rect(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.grid = element_line(size = 0.25),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 7, 0, 7)
      )
    domain_legend <- get_legend(
      test_domains_plot +
      theme(legend.box.margin = margin(7, 7, 7, 7),
      legend.direction = "horizontal",
      legend.background = element_blank()) +
      guides(
        fill=guide_legend(
          nrow=min(3, length(plot_domains$Label))
          )
        )
      )
    plot_title <- ggdraw() +
      draw_label(
        glue("{GeneIter} - {tissue}"),
        #GeneIter,
        fontface = 'bold.italic',
        x = 0,
        hjust = 0,
        size = font_sizes[["medium"]]
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7),
      )
    aligned <- cowplot::plot_grid(
      plot_title,
      test_plot,
      test_domains_plot + theme(legend.position="none"),
      test_plot_2,
      domain_legend,
      align='v', axis="lr", rel_heights = c(1,4,1,3,3), ncol = 1)
    domain_plot_list[[GeneIter]] <- aligned
    
    
    # ggsave(plot=aligned, filename = paste(out_prefix(),".",GeneIter,".plot.pdf", sep = ""),
    #        width=8.5, height=8.5, units="in")
  }
  return(domain_plot_list)
}
