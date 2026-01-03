##### Theme & COLOR RScript (Temp) 
matrix_plot_line_color <- "#e8e8e8"

matrix_plots_theme <- theme(strip.text = element_blank(), 
                            strip.background = element_blank(), 
                            plot.margin = margin(0,0,0,0), 
                            axis.title=element_blank(), 
                            axis.text.y=element_markdown(
                              family = font_family, 
                              size = font.size,
                              color = "black",margin = margin(0,7,0,0)), 
                            panel.spacing = unit(0,"pt")
                            # panel.border = element_rect(color = "black", fill = NA)
)

# legend theme
matrix_legend_theme <- theme(legend.box.margin = margin(5, 5, 5, 5), 
                             legend.background = element_blank(), 
                             legend.text = font_small, 
                             legend.direction = "horizontal",
                             legend.title = font_small, 
                             legend.key = element_rect(color = "black", linewidth=0.125),
                             legend.key.size = unit(font.size, "pt"), 
                             legend.margin = margin(5,5,5,5))

matrix_legend_theme_no_title <- theme(legend.box.margin = margin(5, 5, 5, 5), 
                                      legend.background = element_blank(), 
                                      legend.text = font_small, 
                                      legend.direction = "horizontal",
                                      legend.title = element_blank(), 
                                      legend.key = element_rect(color = "black", linewidth=0.125),
                                      legend.key.size = unit(font.size, "pt"), 
                                      legend.margin = margin(5,5,5,5))

# Line formatting
matrix_v_lines <- geom_vline(
  data=tibble(
    Group=factor(
      c("1", 
        "1", 
        "2", 
        "2", 
        "3",
        "3"),
      levels = c("BRCA1","BRCA2","wt")), 
    x=c(-1,1,-1,1,-1,1)*Inf), 
  aes(xintercept=x), col="black", linewidth = .125)

matrix_h_lines_top <- geom_hline(
  data=tibble(
    Group=factor(
      c("BRCA1","BRCA2","wt"),
      levels = c("BRCA1","BRCA2","wt")), 
    x=c(1,1,1)*Inf), 
  aes(yintercept=x), col="black", linewidth = .125)

matrix_h_lines_bottom <- geom_hline(
  data=tibble(
    Group=factor(
      c("BRCA1","BRCA2","wt"),
      levels = c("BRCA1","BRCA2","wt")), 
    x=c(-1,-1,-1)*Inf), 
  aes(yintercept=x), col="black", linewidth = .125)


# AM score gradiant ####
amscore.breaks = c(
  0,
  (0.564-0.34)/2,
  1
)
amscore.values = c(
  "#a4f5ef",
  "#FAEEC1",
  "#ed211a"
)
scale.amscore.fill = scale_fill_gradientn(colors = amscore.values,
                                          values = amscore.breaks
)

Col.scale.tissue_exp.breaks <- c("COSMIC","PAP","BLO","No Selection")
Col.scale.tissue_exp.values <- c("#A74C56","#FFA892","#CED5E6","#BFBFBF")
Col.scale.tissue_exp.color <- scale_color_manual(
  breaks = Col.scale.tissue_exp.breaks, 
  values = Col.scale.tissue_exp.values,
  name = NULL)
Col.scale.tissue_exp.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.tissue_exp.breaks, 
  values = Col.scale.tissue_exp.values,
  name = NULL)


Col.amClass.breaks <- c("likely_benign","likely_benign_LC",
                        "ambiguous","ambiguous_LC",
                        "likely_pathogenic","likely_pathogenic_LC")
Col.amClass.values <- c(
  "#DAE1F1",
  "#92A8D6",
  "#FAEEC1",
  "#F4D870",
  "#EBB28A",
  "#E08143"
)

# Col.amClass.values <- c(
#   "#ABE2FD",
#   "#0588CA",
#   "#DDCC77",
#   "#A9942D",
#   "#E4A2AD",
#   "#CC6677"
# )

Col.amClas.labels <- c(
  "Likely benign",
  "Likely benign large clone",
  "Ambiguous",
  "Ambiguous large clone",
  "Likely pathogenic",
  "Likely pathogenic large clone"
)

# Col.amClass.fill <- scale_fill_manual(
#   breaks = Col.amClass.breaks,
#   values = Col.amClass.values,
#   labels = Col.amClas.labels,
#   na.value = NA,
#   name = "AlphaMissense (AM)\nPathogenicity\nScore"
# )
Col.amClass.fill <- scale_fill_manual(
  breaks = Col.amClass.breaks,
  values = Col.amClass.values,
  labels = Col.amClas.labels,
  na.value = NA,
  name = "AlphaMissense pathogenicity"
)



# Renaming Facet Labels (Sample Type and Germline Mut!) ####
samp_type.labs <- c("Paps", "Peritoneal Fluid", "Blood")
names(samp_type.labs) <- c("PAP", "PFL","BLO")

germ_mut.labs <- c("<i>BRCA1</i>", "<i>BRCA2</i>", "<i>BRCA wt</i>")
names(germ_mut.labs) <- c("BRCA1", "BRCA2","BRCAwt")

unambigousColorScale <- c(
  "#000000",
  "#E69D00",
  "#56B4DF",
  "#009E74",
  "#F0E442",
  "#0071B2",
  "#D55C00",
  "#CC79A7"
)

#### Coding vs Non-coding Colors ####
Col.scale.coding.breaks <- c("coding", "non-coding")
Col.scale.coding.values <- unambigousColorScale[0:2]
Col.scale.coding.labels <- c("Coding", "Non-coding")
Col.scale.coding.color <- scale_color_manual(
  breaks = Col.scale.coding.breaks, 
  values = Col.scale.coding.values,
  labels = Col.scale.coding.labels,
  name = NULL)
Col.scale.coding.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.coding.breaks, 
  values = Col.scale.coding.values,
  labels = Col.scale.coding.labels,
  name = NULL)


#### Germline Mut Colors ####
Col.scale.germline_mut.breaks <- c("BRCA1", "BRCA2", "BRCAwt")
Col.scale.germline_mut.values <- c("#F8A16D",
                                   "#45CBCE",
                                   "#AB88FB")
Col.scale.germline_mut.labels <- c("BRCA1", "BRCA2", "BRCAwt")
Col.scale.germline_mut.color <- scale_color_manual(
  breaks = Col.scale.germline_mut.breaks, 
  values = Col.scale.germline_mut.values,
  labels = Col.scale.germline_mut.labels,
  name = "Germline Mutation Status")
Col.scale.germline_mut.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.germline_mut.breaks, 
  values = Col.scale.germline_mut.values,
  labels = Col.scale.germline_mut.labels,
  name = "Germline Mutation Status")


#### Age Colors ####
Col.scale.Age.fill <- scale_fill_gradient(na.value = 'white', 
                                          low = "#F7F7E2", 
                                          high = "#020ff5", 
                                          name="Age")


# Breast Cancer History Colors ####
Col.scale.prev_breast.breaks <- c("no","yes")
Col.scale.prev_breast.values <- c("#fffee3", "#fc9e2b")
Col.scale.prev_breast.labels <- c("No", "Yes")
Col.scale.prev_breast.color <- scale_color_manual(
  breaks = Col.scale.prev_breast.breaks, 
  values = Col.scale.prev_breast.values,
  labels = Col.scale.prev_breast.labels,
  name = NULL)
Col.scale.prev_breast.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.prev_breast.breaks, 
  values = Col.scale.prev_breast.values,
  labels = Col.scale.prev_breast.labels,
  name = "Prior\nBreast\nCancer")

# Chemotherapy History Colors ####
Col.scale.prev_chemo.breaks <- c("no","yes", "unknown")
Col.scale.prev_chemo.values <- c("#fcd7f0", "#ff29b8","darkgrey")
Col.scale.prev_chemo.labels <- c("No", "Yes", "Unknown")
Col.scale.prev_chemo.color <- scale_color_manual(
  breaks = Col.scale.prev_chemo.breaks, 
  values = Col.scale.prev_chemo.values,
  labels = Col.scale.prev_chemo.labels,
  name = NULL)
Col.scale.prev_chemo.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.prev_chemo.breaks, 
  values = Col.scale.prev_chemo.values,
  labels = Col.scale.prev_chemo.labels,
  name = "Prior\nChemotherapy")


#### Ca125 Colors ####
Col.scale.Ca125.fill <- scale_fill_gradient(na.value = 'white', 
                                            low = "#e3f3fa", 
                                            high = "#07b1fa", 
                                            name="Pre-operative CA-125")

# Mut Count ####
col.scale.mutCount.fill <- scale_fill_gradient(
  na.value = 'white', 
  low = "#DEEBF7", 
  high = "#2E75B6", 
  limits = c(0, 20),
  expand = c(0,0), 
  breaks = c(0,5,10,15,20), 
  name="Coding Mutations")

# Depth Scale Color ####
Col.scale.dp.fill <- scale_fill_gradient(na.value = "white", 
                                         low = "#E2F0D9", 
                                         high = "#73b54a", 
                                         limits = c(0, 16000),
                                         expand = c(0,0), 
                                         breaks = c(0, 3000,6000,9000,12000,15000), 
                                         labels = c("0","3k","6k","9k","12k","15k"),
                                         name="Mean Depth")

# MBF Scale Color ####
Col.scale.MF.fill <- scale_fill_gradient(na.value = "white", 
                                         low = "#E7E6E6",
                                         high ="#808080",
                                         limits = c(5e-8, 5e-5), 
                                         trans = "log",
                                         expand = c(0,0), 
                                         breaks = c(5e-8,5e-7,5e-6,5e-5), 
                                         labels = c("5E-8","5E-7","5E-6","5E-5"), 
                                         name="Mutation Burden")

# Path Mut Count ####
col.scale.pathCount.fill <- scale_fill_gradient(
  na.value = 'white', 
  low = "#e3fcfc", 
  high = "#f74a4a", 
  limits = c(0, 20),
  expand = c(0,0), 
  breaks = c(0,5,10,15,20), 
  name="Pathogenic Mutations")

# Path MB scale color
Col.scale.path_MB.fill <- scale_fill_gradient(na.value = 'white', 
                                              low = "#ffeded",
                                              high ="#f55858",
                                              limits = c(5e-8, 5e-5), 
                                              trans = "log",
                                              expand = c(0,0), 
                                              breaks = c(5e-8,5e-7,5e-6,5e-5), 
                                              labels = c("5E-8","5E-7","5E-6","5E-5"), 
                                              name="Pathogenic Mutation Burden")


# Path MB scale color
Col.scale.path_MF.fill <- scale_fill_gradient(na.value = 'red', 
                                              low = "white",
                                              high ="black",
                                              limits = c(0, 3.6e-6),
                                              expand = c(0,0), 
                                              breaks = c(0,1.8e-6,3.6e-6), 
                                              labels = c("0","1.8E-6","3.6E-6"), 
                                              name="Pathogenic\nMutation\nFrequency")

# Pathogenic Colors ####
Col.scale.path.breaks <- c("likely_pathogenic",
                           "ambiguous",
                           "likely_benign")
Col.scale.path.values <- c("#ed211a",
                           "#ebc14d",
                           "#a4f5ef")
Col.scale.path.color <- scale_color_manual(
  breaks = Col.scale.path.breaks, 
  values = Col.scale.path.values,
  labels = Col.scale.path.breaks,
  name = NULL)
Col.scale.path.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.path.breaks, 
  values = Col.scale.path.values,
  labels = Col.scale.path.breaks,
  name = NULL)


# Path + LC coloring line ####
Col.scale.mut.breaks <- c("1","2", "3","4")
Col.scale.mut.values <- c("#a10202", "#f74a4a", "#f0d148" ,"#cfeffa")
Col.scale.mut.labels <- c("Pathogenic Large Clone", "Pathogenic",
                          "Non-Pathogenic Large Clone", "Non-Pathogenic")
Col.scale.mut.color <- scale_color_manual(
  breaks = Col.scale.mut.breaks, 
  values = Col.scale.mut.values,
  labels = Col.scale.mut.labels,
  name = NULL)
Col.scale.mut.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.mut.breaks, 
  values = Col.scale.mut.values,
  labels = Col.scale.mut.labels,
  name = "Legend")

# for mutational spectrum plot DON'T CHANGE COLORS####
Col.scale.spectrum.breaks <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
Col.scale.spectrum.values <- c("#5abdeb", "#050708", "#d43c32", "#cbcacb", "#aacb72", "#e7c9c6")

Col.scale.spectrum.fill <- scale_fill_manual(
  breaks=Col.scale.spectrum.breaks,
  values=Col.scale.spectrum.values, 
  name = "Mutation Type")

# for mutational type ####
Col.scale.mut_type.fill <- scale_fill_manual(
  breaks=c("Indel",
           "Splice",
           "Nonsense_Mutation",
           "Missense_Mutation",
           "Silent"), 
  labels=c("Indel",
           "Splice",
           "Nonsense",
           "Missense",
           "Silent"),
  values = c("#320d0f",
             "#f6a3c9",
             "#f2d28b",
             "#887675",
             "#f3eae8"), 
  name = "Mutation Type")
Col.scale.mut_type.color <- scale_color_manual(
  breaks=c("Indel",
           "Splice",
           "Nonsense",
           "Missense",
           "Silent"), 
  values = c("#320d0f",
             "#f6a3c9",
             "#f2d28b",
             "#887675",
             "#f3eae8"), 
  name = "Mutation Type")


# hotspot proportion plot ####
Col.scale.hotspot.breaks <- c("COSMIC","BRCA1","BRCA2","BRCAwt","No Selection")
Col.scale.hotspot.values <- c("#A74C56","#F8A16D","#45CBCE","#AB88FB","#8C7F7D")
Col.scale.hotspot.labels <- c("COSMIC", "<i>BRCA1</i>", "<i>BRCA2</i>","<i>BRCA</i> WT",
                              "No Selection")
Col.scale.hotspot.color <- scale_color_manual(
  breaks = Col.scale.hotspot.breaks, 
  values = Col.scale.hotspot.values,
  labels = Col.scale.hotspot.labels,
  name = NULL)
Col.scale.hotspot.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.hotspot.breaks, 
  values = Col.scale.hotspot.values,
  labels = Col.scale.hotspot.labels,
  name = "Germline Mutation")

Col.scale.hotspot2.breaks <- c("Hotspot","Not Hotspot")
Col.scale.hotspot2.values <- c("#003870","gray")
Col.scale.hotspot2.textvalues <- c("white","gray")

Col.scale.hotspot2.color <- scale_color_manual(
  breaks = Col.scale.hotspot2.breaks, 
  values = Col.scale.hotspot2.textvalues,
  name = NULL,
  guide = "none")
Col.scale.hotspot2.fill <- scale_fill_manual(
  na.value = NA,
  breaks = Col.scale.hotspot2.breaks, 
  values = Col.scale.hotspot2.values,
  name = "Germline Mutation")

Col.scale.is_neoplasia.breaks = c("N","Y")
Col.scale.is_neoplasia.values = c("#E5F4F8","#2F74B5")
Col.scale.is_neoplasia.labels = c("No","Yes")
Col.scale.is_neoplasia.fill = scale_fill_manual(
  breaks = Col.scale.is_neoplasia.breaks,
  values = Col.scale.is_neoplasia.values,
  labels = Col.scale.is_neoplasia.labels,
  name = "Benign\nNeoplasia"
)

ageRaceColors = c(
  "No Selection"="#DDCC77",
  "Normal 0-29 White"="#117733",
  "Normal 0-29 Black"="#44AA99",
  "Normal 30-65 White"="#CC6677",
  "Normal 30-65 Black"="#AA4499",
  "Cancer"="#882235",
  "COSMIC"="#332288"
)

scale.ageRace.fill <-
  scale_fill_manual(
    values = ageRaceColors,
    name = ""
  )
scale.ageRace.color <-
  scale_color_manual(
    values = ageRaceColors,
    name = ""
  )
ageRaceColors2 = c(
  "No Selection"="#DDCC77",
  "No Cancer 0-29 White"="#117733",
  "No Cancer 0-29 Black"="#44AA99",
  "No Cancer 30-65 White"="#CC6677",
  "No Cancer 30-65 Black"="#AA4499",
  "Cancer"="#882235",
  "COSMIC"="#332288"
)
scale.ageRace2.fill <-
  scale_fill_manual(
    values = ageRaceColors2,
    name = ""
  )
scale.ageRace2.color <-
  scale_color_manual(
    values = ageRaceColors2,
    name = ""
  )

ageBidecadeColors = c(
  "No Selection"="#DDCC77",
  "0-19" = "#BBBBBB",
  "20-39" = "#999999",
  "40-59" = "#777777",
  "60+" = "#555555",
  "Cancer"="#882235",
  "COSMIC"="#332288"
)
scale.ageBidecade.fill <-
  scale_fill_manual(
    values = ageBidecadeColors,
    name = ""
  )

