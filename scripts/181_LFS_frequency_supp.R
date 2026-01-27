## Frequency of LFS mutation by tissue type

cancer_samples = c("All tissue-types",
                   "Mediastinal metastasis",
                   "Lung metastasis",
                   "Esophageal cancer 1",
                   "Esophageal cancer 2",
                   "Liver metastasis 1",
                   "Liver metastasis 2")
tissue_order <- c("Whole blood", 
                  "Buffy coat", 
                  "Plasma", 
                  "Bone marrow", 
                  #"Buccal mucosa", 
                  "Thyroid", 
                  "Mainstem bronchus",
                  "Lung", 
                  "Esophagus 1", 
                  "Esophagus 2", 
                  "Gastric 1",
                  "Gastric 2",
                  "Cardiac muscle",
                  "Spleen",
                  "Liver",
                  "Colon",
                  "Omentum",
                  "Peritoneum",
                  "Renal",
                  "Testis",
                  "Skeletal muscle",
                  "Skin",
                  "Skin, non-sun-exposed",
                  "Mediastinal metastasis",
                  "Lung metastasis",
                  "Esophageal cancer 1",
                  "Esophageal cancer 2",
                  "Liver metastasis 1",
                  "Liver metastasis 2")
abbreviations <- c("WB", "Buffy", "Plas", "BM", "Thyr", "Bron", "Lung",
                   "Eso1", "Eso2", "Gast1", "Gast2", "CardM", "Spln", "Liver", "Colon",
                   "Omen", "Perit", "Renal", "Testis", "SkelM", "Skin", "SkinNS",
                   "MedMet", "LungMet", "EsoCa1", "EsoCa2", "LivMet1", "LivMet2")



LFS_frequency <- MAF_table %>% 
  filter(prot.pos == "181") %>% 
  filter(Subject == "Patient") %>% 
  filter(Hugo_Symbol == "TP53") %>%
  filter(Tissue != "Buccal mucosa") %>%
  dplyr::select(Hugo_Symbol, Start_Position, Variant_Classification, Reference_Allele, Tumor_Seq_Allele2, t_depth, protein_variant, Tissue, VAF) %>%
  arrange(desc(VAF)) %>%
  mutate(Tissue = reorder(Tissue, VAF, FUN = max)) %>%
  mutate(Tumor = ifelse(Tissue %in% cancer_samples, "Tumor", "Non-Tumor"))
tissue_labels <- abbreviations
names(tissue_labels) <- tissue_order


custom_label <- function(x) {
  sapply(x, function(t) {
    abbrev <- tissue_labels[[t]]
    if (t %in% cancer_samples) {
      paste0("<span style='color:red;'>", abbrev, "</span>")
    } else {
      abbrev
    }
  })
}

LFS_frequency
plot <- ggplot(LFS_frequency, aes(x=Tissue, y = VAF, color = Tumor)) + 
  geom_point() + 
  labs(y = "Germline p.R181H VAF") +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = "black") + 
  theme_minimal() + 
  scale_x_discrete(labels = custom_label) +
  scale_y_continuous(limits = c(0.4,1), breaks = seq(0,1,by = 0.1)) +
  theme(axis.text.x.bottom = element_markdown(angle = 90, size = 8, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =8),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = c(0.2,0.8)) + 
  scale_color_manual(name = NULL, values = c("Tumor" = "red", "Non-Tumor" = "black")) 
  
plot
ggsave("results/Manuscript_figures/Fig_S5/LFS_VAF_slide.png", plot, width = 4, height = 2.5, units = "in", dpi = 300)
