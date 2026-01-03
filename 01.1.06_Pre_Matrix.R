library("tidyverse")
library("cowplot")
ggMatrixLine = function(data, x_col, column, y_name = NA, fill_scale=NA,
                        facet_var=NA,
                        plot_text=FALSE, text_size = 8) {
  
  if (is.na(y_name)) {
    y_name = column
  }
  p <- ggplot(data, aes_string(x=x_col, fill=column, label = column)) +
    geom_tile(aes(y=y_name))
  if (! is.na(plot_text)) {
    if (plot_text) {
    p <- p + geom_text(aes(y=y_name), color = "black",
                       size = text_size *  25.4 / 72.27 )
  }}
  if (! is.na(fill_scale)) {
    p <- p + fill_scale
  }
  if (! is.na(facet_var)) {
    p <- p + facet_grid(reformulate(facet_var, "."), scales="free_x", space = "free_x")
  }
  p <- p +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      plot.margin = unit(c(-.1, .1, -.1, .1), "cm"),
      axis.title=element_blank(),
      axis.text.y=element_text(
        family = "sans",
        size = text_size,
        color = "black",margin = margin(0,0,0,0,unit = "pt")),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      panel.spacing = unit(0,"pt"),
      plot.background = element_blank(),
      
      
    )
  return (p)
}

gg_matrix_legend_extractor <- function(in_plot, text_size = 8) {
  return (
    get_legend(
      in_plot +
        theme(legend.box.margin = margin(7, 7, 7, 7),
              legend.background = element_blank(),
              legend.text = element_text(
                family = "sans",
                size = text_size,
                color = "black"),
              legend.direction = "horizontal",
              legend.title = element_text(
                family = "sans",
                size = text_size,
                color = "black"),
              legend.key = element_rect(color = "black"),
              legend.key.size = unit(text_size, "pt"),
              legend.key.width = unit(2*text_size, "pt"),
              legend.margin = margin(3,3,3,3)
        )
    )
  )
}

gg_matrix_preplot <- function(in_plot) {
  return(
    in_plot + theme(legend.position = "none")
  )
}

ggSampleMatrix = function(data=NULL,
                          x_col=NULL,
                          cols_to_plot=NULL,
                          y_names=NULL,
                          facet_var=NULL,
                          fill_scales=NULL,
                          plot_text=NULL,
                          legend_cols = 2,
                          show_legend = NULL,
                          text_size = 8,
                          return_type = "plot") {
  if (is.null(data) |
      is.null(x_col) |
      is.null(cols_to_plot) |
      ! is.vector(cols_to_plot) |
      ! x_col %in% colnames(data)) {
    stop("Err1")
  }
  if (! is.null(fill_scales) & ! is.vector(fill_scales) & length(fill_scales) != length(cols_to_plot)) {
    stop("Err2")
  }
  if (! is.null(plot_text) & ! is.vector(plot_text) & length(plot_text) != length(cols_to_plot)) {
    stop("Err3")
  }
  # create subplots
  subplots <- list()
  plot_sizes <- c()
  for (i in seq(length(cols_to_plot))) {
    print(i)
    print(cols_to_plot[i])
    subplots[[i]] <- ggMatrixLine(data = data,
                                  x_col = x_col,
                                  column = cols_to_plot[i],
                                  y_name = if (! is.null(y_names)) {y_names[i]} else {NA},
                                  fill_scale = if (! is.null(fill_scales)) {fill_scales[i]} else {NA},
                                  facet_var = if (! is.null(facet_var)) {facet_var} else {NA},
                                  plot_text = if (! is.null(plot_text)) {plot_text[i]} else {NA},
                                  text_size = text_size
    )
    plot_sizes <- c(plot_sizes, 1)
  }
  
  # extract legends
  if (is.null(show_legend)) {
    matrix_legends <- cowplot::plot_grid(
      plotlist = lapply(subplots, gg_matrix_legend_extractor,
                        text_size = text_size),
      align = 'v', ncol = legend_cols, axis = 'l'
    )
  }  else if (all(show_legend == F)) {
    matrix_legends = NULL
  } else {
    matrix_legends <- cowplot::plot_grid(
      plotlist = lapply(subplots[show_legend == TRUE], gg_matrix_legend_extractor,
                        text_size = text_size),
      align = 'v', ncol = legend_cols, axis = 'l'
    )
  }
  
  if (return_type == "plot") {
    subplots[[i + 1]] <- matrix_legends
    out_plot <- cowplot::plot_grid(
      plotlist = lapply(
        subplots, gg_matrix_preplot
      ),
      rel_heights = c(plot_sizes, length(cols_to_plot) / 2),
      ncol=1, align = 'v', axis = 'lr'
    )
    return(out_plot)
  } else if (return_type == "subplot") {
    return(subplots)
  }
}


# sample_matrix_data <- sample_data %>% 
#   left_join(
#     maf_table%>% 
#       select(PatientNum, Sample_type) %>% 
#       unique()
#   ) %>% 
#   mutate(Sample = paste(PatientNum, Sample_type, sep="_"),
#          Study_Group_Num=factor(Study_Group_Num), 
#          Prior_Pos_Pap = factor(Prior_Pos_Pap), 
#          StageSimple = factor(StageSimple))
# 
# matrix_colors_Sample <- sample_matrix_data %>% 
#   select(Sample) %>% 
#   mutate(color="white") %>% 
#   unique()
# matrix_colors_PatientNum <- sample_matrix_data %>% 
#   select(PatientNum) %>% 
#   mutate(color="white") %>% 
#   unique()
# test_matrix <-ggSampleMatrix(
#   data=sample_matrix_data,
#   x_col = "Sample",
#   cols_to_plot = c("Study_Group_Num",
#                    "Sample",
#                    "PatientNum",
#                    "Sample_type",
#                    "Age_At_Surg",
#                    "Cancer",
#                    "BRCA12",
#                    "Smoking_History",
#                    "Prior_Pos_Pap",
#                    "PrevChemo",
#                    "StageSimple"),
#   plot_text =   c(F,T,T,F,F,F,F,F,F,F,F),
#   show_legend = c(T,F,F,T,T,T,T,T,T,T,T),
#   facet_var = "Study_Group_Num",
#   text_size = 6,
#   fill_scales = c(
#     scale_fill_discrete(),
#     scale_fill_manual(
#       values=matrix_colors_Sample$color
#     ),
#     scale_fill_manual(
#       breaks=matrix_colors_PatientNum$PatientNum,
#       values=matrix_colors_PatientNum$color
#     ),
#     NA,
#     scale_fill_gradient(low="white", high="cyan"),
#     NA,NA,NA,NA,NA,NA)
# )
# test_matrix
