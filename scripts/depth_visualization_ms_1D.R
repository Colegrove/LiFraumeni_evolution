#### Figure 1D 
#### sequencing depth visualization
#### To generate depth panel file first run depth_by_panel_S1.R


## sequencing depth information
out_path <- "results/Manuscript_tables/"
out_file <- "supplemental_depth_panels.csv"
panel_depths <- read_delim(paste0(out_path,out_file))


metadata <- tribble(
  ~Subject_abbr, ~Subject, ~Tissue, ~GenePanel,
  "LFS01","Patient", "Whole blood", "TP53+MUT+CHIP",
  "LFS01","Patient", "Buffy coat", "TP53+MUT+CHIP",
  "LFS01","Patient", "Plasma", "TP53+MUT+CHIP",
  "LFS01","Patient", "Bone marrow", "TP53+MUT+CHIP",
  #"LFS01","Patient", "Buccal mucosa", "TP53", 16138,
  "LFS01","Patient", "Thyroid", "TP53",
  "LFS01","Patient", "Mainstem bronchus", "TP53",
  "LFS01","Patient", "Lung", "TP53",
  "LFS01","Patient", "Esophagus 1", "TP53+MUT",
  "LFS01","Patient", "Esophagus 2", "TP53",
  "LFS01","Patient", "Gastric 1", "TP53",
  "LFS01","Patient", "Gastric 2", "TP53",
  "LFS01","Patient", "Cardiac muscle", "TP53",
  "LFS01","Patient", "Spleen", "TP53",
  "LFS01","Patient", "Liver", "TP53+MUT",
  "LFS01","Patient", "Colon", "TP53+MUT",
  "LFS01","Patient", "Omentum", "TP53",
  "LFS01","Patient", "Peritoneum", "TP53",
  "LFS01","Patient", "Renal", "TP53",
  "LFS01","Patient", "Testis", "TP53",
  "LFS01","Patient", "Skeletal muscle", "TP53",
  "LFS01","Patient", "Skin", "TP53+MUT",
  "LFS01","Patient", "Skin, non-sun-exposed", "TP53+MUT",
  "LFS01","Patient", "Esophageal cancer 1", "TP53+MUT",
  "LFS01","Patient", "Esophageal cancer 2", "TP53+MUT",
  "LFS01","Patient", "Liver metastasis 1", "TP53+MUT",
  "LFS01","Patient", "Liver metastasis 2", "TP53+MUT",
  "LFS01","Patient", "Lung metastasis", "TP53+MUT",
  "LFS01","Patient", "Mediastinal metastasis", "TP53+MUT",
  "LFS02","Family member A", "PBMC", "TP53+MUT+CHIP",
  "REL01","Family member B", "PBMC", "TP53+MUT+CHIP",
  "LFS03","Family member C", "PBMC", "TP53+MUT+CHIP",
  "CON01","UW volunteer 1", "PBMC", "TP53+MUT+CHIP",
  "CON02","UW volunteer 2", "PBMC", "TP53+MUT+CHIP", 
  "CON03","UW volunteer 3", "PBMC", "TP53+MUT+CHIP", 
  "CON04","UW volunteer 4", "PBMC", "TP53+MUT+CHIP",
  "CON05","UW volunteer 5", "PBMC", "TP53+MUT+CHIP", 
  "CON06","UW volunteer 6", "PBMC", "TP53+MUT+CHIP", 
  "CON07","UW volunteer 7", "PBMC", "TP53+MUT+CHIP",
)

tissue_order <- c(
  "Bone marrow", "Whole blood","Plasma", "Buffy coat", "PBMC",
  "Buccal mucosa", "Thyroid", "Mainstem bronchus", "Lung",
  "Esophagus 1", "Esophagus 2", "Gastric 1", "Gastric 2",
  "Cardiac muscle", "Spleen", "Liver", "Colon",
  "Omentum", "Peritoneum", "Renal", "Testis",
  "Skeletal muscle", "Skin", "Skin, non-sun-exposed",
  "Esophageal cancer 1", "Esophageal cancer 2",
  "Liver metastasis 1", "Liver metastasis 2",
  "Lung metastasis", "Mediastinal metastasis"
)
## order option 2 to group seq panels better
tissue_order <- c(
  "Bone marrow", "Whole blood","Plasma", "Buffy coat", "PBMC",
  
  "Esophageal cancer 1", "Esophageal cancer 2",
  "Liver metastasis 1", "Liver metastasis 2",
  "Lung metastasis", "Mediastinal metastasis",
  "Liver", "Colon", "Skin", "Skin, non-sun-exposed", 
  "Esophagus 1", "Esophagus 2", "Gastric 1", "Gastric 2",
  "Thyroid", "Mainstem bronchus", "Lung",
  "Cardiac muscle", "Spleen",
  "Omentum", "Peritoneum", "Renal", "Testis",
  "Skeletal muscle"
)


tissue_abbreviations <- tribble(
  ~Tissue,                   ~Tissue_abbr,
  "Whole blood",             "WB",
  "Buffy coat",              "Buffy",
  "Plasma",                  "Plasma",
  "Bone marrow",             "BM",
  #"Buccal mucosa",           "Bucc",
  "Thyroid",                 "Thyroid",
  "Mainstem bronchus",       "Bronchus",
  "Lung",                    "Lung",
  "Esophagus 1",             "Esoph1",
  "Esophagus 2",             "Esoph2",
  "Gastric 1",               "Gast1",
  "Gastric 2",               "Gast2",
  "Cardiac muscle",          "Cardiac",
  "Spleen",                  "Spleen",
  "Liver",                   "Liver",
  "Colon",                   "Colon",
  "Omentum",                 "Omentum",
  "Peritoneum",              "Peritoneum",
  "Renal",                   "Renal",
  "Testis",                  "Testis",
  "Skeletal muscle",         "Skeletal",
  "Skin",                    "Skin",
  "Skin, non-sun-exposed",   "SkinNS",
  "Mediastinal metastasis",  "Med Met",
  "Lung metastasis",         "Lung Met",
  "Esophageal cancer 1",     "Esoph Ca1",
  "Esophageal cancer 2",     "Esoph Ca2",
  "Liver metastasis 1",      "Liver Met1",
  "Liver metastasis 2",      "Liver Met2",
  "PBMC",                    "PBMC"
)

metadata <- metadata %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))
# add abbreviations
metadata <- metadata %>%
  left_join(tissue_abbreviations, by="Tissue")
metadata

# add subject/tissue information to depth data
tissue_information <- filt_maf %>%
  distinct(Tumor_Sample_Barcode, Subject, Tissue) %>%
  filter(!(is.na(Subject))) %>%
  left_join(metadata, by = c("Subject", "Tissue")) %>%
  filter(Tissue != "Urine cells")

## summarize sequencing depth by individual and attach to metadata
subject_seq_depth <- panel_depths %>%
  group_by(Subject, Tissue) %>%
  summarise(DP = mean(MeanDepth))
tissue_information_depths <- tissue_information %>% 
  left_join(subject_seq_depth, by=c(Subject_abbr = "Subject", Tissue_abbr = "Tissue"))


tissue_order_df <- tibble(Tissue = tissue_order, 
                          OrderID = seq_along(tissue_order))
tissue_information_depths <- tissue_information_depths %>%
  left_join(tissue_order_df, by="Tissue")


## make cancer red
cancer_sample_abbr <- c("Med Met","Lung Met","Esoph Ca1","Esoph Ca2","Liver Met1","Liver Met2")
blood_samps <- c("PBMC", "Buffy coat", "Plasma", "Whole blood", "Bone marrow")
blood_samps_abbr <-c("Buffy", "Plasma", "WB", "BM")

tissue_information_depths_grouped <- tissue_information_depths %>%
  mutate(
    TissueGroup = if_else(Tissue %in% blood_samps | grepl("PBMC", Tissue, ignore.case = TRUE),
                          "Blood",
                          "Solid tissues"),
    SampleLabel = case_when(Tissue_abbr %in% cancer_sample_abbr ~ paste0("<span style='color:red;'>", Tissue_abbr), 
                            Tissue_abbr %in% blood_samps_abbr ~ Tissue_abbr,
                            Tissue_abbr == "PBMC" ~ paste0(Subject_abbr,"_", Tissue_abbr),
                            TRUE ~ Tissue_abbr),
    
    SampleLabel = fct_reorder(SampleLabel, OrderID, .fun = min) 
  ) %>% 
  mutate(GenePanel = factor(GenePanel, levels = c("TP53+MUT+CHIP","TP53+MUT","TP53")))

tissue_information_depths_grouped 

#######
## bar plot version
#######
depth_plot_split <- ggplot(tissue_information_depths_grouped,
                           aes(x = SampleLabel, y = DP, fill = GenePanel)) +
  geom_col() +
  scale_fill_manual(values = c( "TP53+MUT+CHIP" = "#000000", "TP53+MUT" = "#555555", "TP53" = "#AAAAAA")) +
  facet_grid(cols = vars(TissueGroup),
             scales = "free_x", space = "free_x", switch = "x") +
  labs(y = "Duplex depth", x = NULL, fill = "Gene panel") +
  scale_y_continuous(breaks = c(0, 10000, 20000)) +
  geom_hline(
    yintercept = 1000,
    #linetype = "dotdash",
    linewidth = 0.5,
    color = "#DDCC77"
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_blank(),
    axis.text.x.bottom = element_markdown(angle = 90, vjust = 0.5, hjust = 1, size = 8, margin = margin(t=8)),
    axis.text.y = element_text(size=8),
    axis.title.y = element_text(size=8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(0,1,-2,1),
  )

depth_plot_split
ggsave("results/depth_schematic.png", depth_plot_split, height = 1.9, width = 4.8)
ggsave("results/Manuscript_figures/Fig_1/depth_schematic.png", depth_plot_split, width = 4.8, height = 1.8, units = "in", dpi = 300)
