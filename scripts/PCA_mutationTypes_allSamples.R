# Hunter Colegrove
# 29 Apr 2025
# Perform a PCA on the 96 mutation types across LF samples

SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/output/SBS/li_fraumeni.SBS96.all"
SigPath2 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/output/SBS/li_fraumeni.SBS96.all"
df <- read_delim(SigPath, delim = "\t") %>% select(-DNA03975_S12.1)
df2 <- read_delim(SigPath2, delim = "\t")
df <- df %>%
  left_join(df2, by = "MutationType")
df %>% print(width = Inf)

sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sampleID_map, delim = "\t") %>%
  mutate(across(where(is.character), trimws)) %>%
  rename(sample = sample)
sample_map
df
new_names <- sample_map$subject[ match(colnames(df)[-1], sample_map$sample) ]
new_names
colnames(df)[-1] <- new_names
df
mat       <- df %>% column_to_rownames("MutationType") %>% as.matrix() %>% t()
mat_prop  <- sweep(mat, 1, rowSums(mat), FUN = "/")
pca   <- prcomp(mat_prop, center = TRUE, scale. = TRUE)
as.data.frame(pca$x)

scores <- as_tibble(pca$x, rownames = "Sample")
scores

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

manual_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:", "Patient:")
legend_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Patient:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")
scores
scores <- scores %>%
  mutate(Sample = factor(Sample, levels = manual_subject_order)) #%>%
  #mutate(subject_age_label = paste0(Sample, " (", age, ")")) 

## color by participant
pca_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = Sample, shape = LFS)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  scale_color_manual(values = custom_colors_age, breaks=legend_subject_order) +
  theme_minimal()
show(pca_plot)

path1 <- "results/PCA_allSamples.png"
ggsave(path1, pca_plot, width = 6, height = 6, units = "in", dpi = 300)
#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)


## color by tissue type

sample_colors <- c(
  "Buffy coat"   = "#67000d",
  "Plasma"       = "#ef3b2c",
  "Whole blood"  = "#cb181d",
  "Bone marrow"  = "#fc9272",
  "PBMC"         = "red",
  "Colon"        = "#6a51a3",
  "Esophagus 1"         = "#9ecae1",
  "Esophageal cancer 1" = "#6baed6",
  "Esophageal cancer 2" = "#2171b5",
  "Liver"               = "#fdb863",
  "Liver metastasis 1"  = "#e08214",
  "Liver metastasis 2"  = "#b35806",
  "Lung metastasis"         = "#636363",
  "Mediastinal metastasis"  = "#252525",
  "Skin"                   = "#a1d99b",
  "Skin, non-sun-exposed"  = "#31a354"
)
tissue_order <- c("Whole blood", "PBMC", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")

SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/output/SBS/li_fraumeni.SBS96.all"
SigPath2 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/output/SBS/li_fraumeni.SBS96.all"
df <- read_delim(SigPath, delim = "\t") %>% select(-DNA03975_S12.1)
df2 <- read_delim(SigPath2, delim = "\t")
df <- df %>%
  left_join(df2, by = "MutationType")
df %>% print(width = Inf)
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sampleID_map, delim = "\t") %>%
  mutate(across(where(is.character), trimws)) %>%
  rename(sample = sample)
sample_map
df
new_names <- sample_map$tissue[ match(colnames(df)[-1], sample_map$sample) ]
new_names
colnames(df)[-1] <- new_names
df
mat       <- df %>% column_to_rownames("MutationType") %>% as.matrix() %>% t()
mat_prop  <- sweep(mat, 1, rowSums(mat), FUN = "/")
pca   <- prcomp(mat_prop, center = TRUE, scale. = TRUE)
as.data.frame(pca$x)

scores <- as_tibble(pca$x, rownames = "Sample")
scores

var_exp <- pca$sdev^2 / sum(pca$sdev^2)
pc1_pct <- round(var_exp[1] * 100, 1)
pc2_pct <- round(var_exp[2] * 100, 1)
LFS_samples <- c("Patient:", "Family member A:", "Family member C:")
scores <- scores %>% 
  mutate(LFS = if_else(Sample %in% LFS_samples, "LFS", "non-LFS"))
scores

manual_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:", "Patient:")
legend_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Patient:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")
scores
scores <- scores %>%
  mutate(Sample = factor(Sample)) #%>%
#mutate(subject_age_label = paste0(Sample, " (", age, ")")) 



pca_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  scale_color_manual(values = sample_colors, breaks= tissue_order) +
  theme_minimal()
show(pca_plot)

path1 <- "results/PCA_allSamples_tissue.png"
ggsave(path1, pca_plot, width = 6, height = 6, units = "in", dpi = 300)
#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)





## 3 PCs
# # 1. Extract scores and % variance explained
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




