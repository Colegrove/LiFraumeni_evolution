# Hunter Colegrove
# 29 Apr 2025
# Perform a PCA on the 96 mutation types across LF samples

SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/output/SBS/li_fraumeni.SBS96.all"
df <- read_delim(SigPath, delim = "\t")
df

sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sampleID_map, delim = "\t") %>%
  mutate(across(where(is.character), trimws)) %>%
  rename(sample = sample)

new_names <- sample_map$tissue[ match(colnames(df)[-1], sample_map$sample) ]
colnames(df)[-1] <- new_names

mat       <- df %>% column_to_rownames("MutationType") %>% as.matrix() %>% t()
mat_prop  <- sweep(mat, 1, rowSums(mat), FUN = "/")
pca   <- prcomp(mat_prop, center = TRUE, scale. = TRUE)

scores <- as.data.frame(pca$x)
scores$Sample <- rownames(scores)

var_exp <- pca$sdev^2 / sum(pca$sdev^2)
pc1_pct <- round(var_exp[1] * 100, 1)
pc2_pct <- round(var_exp[2] * 100, 1)
cancer_samples <- c("Mediastinal metastasis", "Lung metastasis", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2")
scores <- scores %>% 
  mutate(cancer = if_else(Sample %in% cancer_samples, "cancer", "non-cancer"))
scores

sample_colors <- c(
  "Buffy coat"   = "#67000d",
  "Plasma"       = "#ef3b2c",
  "Whole blood"  = "#cb181d",
  "Bone marrow"  = "#fc9272",
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
tissue_order <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")

pca_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = Sample, shape = cancer)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  scale_shape_manual(values=c(15,17)) +
  scale_color_manual(values = sample_colors, breaks = tissue_order) +
  theme_minimal()
show(pca_plot)

path1 <- "results/PCA_patient_tissues_allPanels.png"
ggsave(path1, pca_plot, width = 6, height = 6, units = "in", dpi = 300)
#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)



library(plotly)

# 1. Extract scores and % variance explained
scores <- as.data.frame(pca$x)
scores$Sample <- rownames(scores)
scores <- scores %>% 
  mutate(cancer = if_else(Sample %in% cancer_samples, "cancer", "non-cancer"))
var_exp <- pca$sdev^2 / sum(pca$sdev^2)
pc1_pct <- round(var_exp[1] * 100, 1)
pc2_pct <- round(var_exp[2] * 100, 1)
pc3_pct <- round(var_exp[3] * 100, 1)
scores
# 2. 3D plot
plot_ly(
  data = scores,
  x    = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Sample,        # remove or change if you don't have a 'tissue' column
  text  = ~Sample,        # hover shows sample name
  marker = list(size = 4), 
  symbol = ~cancer, 
  symbols = c("circle", "x")
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1 (", pc1_pct, "%)")),
      yaxis = list(title = paste0("PC2 (", pc2_pct, "%)")),
      zaxis = list(title = paste0("PC3 (", pc3_pct, "%)"))
    ),
    title = "3D PCA of Mutation-Type Proportions"
  )

df_proportions <- df %>% 
  select("MutationType", "Patient:") %>%
  rename(Patient = "Patient:") %>%
  mutate(proportion = Patient/sum(Patient)) %>%
  #print(n=Inf) %>%
  filter(MutationType == "C[C>T]G" | MutationType == "C[C>T]G") %>%
  print()

