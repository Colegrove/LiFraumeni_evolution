# Hunter Colegrove
# 29 Apr 2025
# Perform a PCA on the 96 mutation types across LF samples

SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/output/SBS/li_fraumeni.SBS96.all"
SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/output/SBS/li_fraumeni.SBS96.region"
#SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels/output/SBS/li_fraumeni.SBS96.region"
df <- read_delim(SigPath, delim = "\t")
df

sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sampleID_map, delim = "\t") %>%
  mutate(across(where(is.character), trimws)) %>%
  rename(sample = sample)
sample_map
df
new_names <- sample_map$subject[ match(colnames(df)[-1], sample_map$sample) ]
new_names
colnames(df)[-1] <- new_names

mat       <- df %>% column_to_rownames("MutationType") %>% as.matrix() %>% t()
mat_prop  <- sweep(mat, 1, rowSums(mat), FUN = "/")
pca   <- prcomp(mat_prop, center = TRUE, scale. = TRUE)

scores <- as.data.frame(pca$x)
scores$Sample <- rownames(scores)

var_exp <- pca$sdev^2 / sum(pca$sdev^2)
pc1_pct <- round(var_exp[1] * 100, 1)
pc2_pct <- round(var_exp[2] * 100, 1)
LFS_samples <- c("Patient:", "Family member A:", "Family member C:")
scores <- scores %>% 
  mutate(LFS = if_else(Sample %in% LFS_samples, "LFS", "non-LFS"))
scores

custom_colors_age <- c(
  "#fcce25","#fcce25","#fcce25","#fcce25","cyan", #"#fca636",
  "#f2844b","#b12a90","#8f0da4","#6a00a8","#0d0887","#0d0887")

custom_colors_age <- c(
  "#fcce25","#fcce25","#fcce25","#fcce25","cyan", #"#fca636",
  "#f2844b","#b12a90","#8f0da4","#6a00a8","#0d0887","#0d0444")

manual_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:", "Patient:")
legend_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Patient:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")
scores
scores <- scores %>%
  mutate(Sample = factor(Sample, levels = manual_subject_order)) #%>%
  #mutate(subject_age_label = paste0(Sample, " (", age, ")")) 

scores
pca_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = Sample, shape = LFS)) +
  geom_point(size = 8.5, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  scale_color_manual(values = custom_colors_age, breaks=legend_subject_order) +
  theme_minimal() + 
  theme(text = element_text(size = 20))
show(pca_plot)

path1 <- "results/PCA_blood_allPanels.png"
path2 <- "results/PCA_blood_bedPanel.png"
ggsave(path1, pca_plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(path2, pca_plot, width = 12, height = 6, units = "in", dpi = 300)
#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)


## 3 PCs
# 1. Extract scores and % variance explained
# scores <- as.data.frame(pca$x)
# scores$Sample <- rownames(scores)
# scores <- scores %>%
#   mutate(LFS = if_else(Sample %in% LFS_samples, "LFS", "non-LFS"))
# var_exp <- pca$sdev^2 / sum(pca$sdev^2)
# pc1_pct <- round(var_exp[1] * 100, 1)
# pc2_pct <- round(var_exp[2] * 100, 1)
# pc3_pct <- round(var_exp[3] * 100, 1)
# scores
# # 2. 3D plot
# plot_ly(
#   data = scores,
#   x    = ~PC1, y = ~PC2, z = ~PC3,
#   color = ~Sample,        # remove or change if you don't have a 'tissue' column
#   text  = ~Sample,        # hover shows sample name
#   marker = list(size = 4),
#   symbol = ~LFS,
#   symbols = c("circle", "x")
# ) %>%
#   add_markers() %>%
#   layout(
#     scene = list(
#       xaxis = list(title = paste0("PC1 (", pc1_pct, "%)")),
#       yaxis = list(title = paste0("PC2 (", pc2_pct, "%)")),
#       zaxis = list(title = paste0("PC3 (", pc3_pct, "%)"))
#     ),
#     title = "3D PCA of Mutation-Type Proportions"
#   )




