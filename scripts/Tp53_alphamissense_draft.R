### TP53 binding domain analysis


subjects_blood <- c("UW volunteer 1","UW volunteer 2","UW volunteer 3","UW volunteer 4",
                    "UW volunteer 5","UW volunteer 6","UW volunteer 7",
                    "Patient","Family member A","Family member C","Family member B")

all_samples = c("All tissue-types",
                "Whole blood", 
                "Buffy coat", 
                "Plasma", 
                "Bone marrow", 
                "Buccal mucosa", 
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

## all possible mutations
all_possible_muts_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/all_possible_sites_annotated.tsv"
all_possible_muts <- read_delim(all_possible_muts_path, delim="\t")

## assign alphamissense to all_possible missense mutations and exclude the rest
all_missense <- all_possible_muts %>%
  filter(GENE == "TP53") %>%
  left_join(alphamissense, by = c("POS" = "POS", "REF" = "REF", "ALT" = "ALT")) %>%
  filter(!is.na(am_pathogenicity))

################################################################################
### LFS01 Tissues
################################################################################
tissue_filtering <- maf_masked_coding %>% 
  filter(Subject == "Patient") %>%
  filter(Tissue %in% all_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP") %>% 
  filter(Variant_Classification == "Missense_Mutation") %>%
  filter(!is.na(am_pathogenicity))

tissue_filtering$group <- "Observed"
all_missense$group <- "Not observed"
combine_df <- bind_rows(tissue_filtering, all_missense)

combine_df %>% filter(group != "Observed") %>% print(width = Inf)
## plot all observed together
obsv_am_tissue <- ggplot(combine_df, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  scale_color_manual(values = c("Not observed" = "grey60", "Observed" = "#3182bd")) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_tissue
## separate large and small observed clones
combine_df_large_clones <- combine_df %>%
  mutate(group = case_when(
    group == "Observed" & t_alt_count == 1 ~ "Observed (1 read)",
    group == "Observed" & t_alt_count > 1  ~ "Observed (>1 read)",
    TRUE ~ group
  ))
combine_df_large_clones$group <- factor(combine_df_large_clones$group,
                           levels = c("Not observed", "Observed (1 read)", "Observed (>1 read)"))

obsv_am_large_clones_tissue <- ggplot(combine_df_large_clones, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c(
    "Not observed" = "grey60",
    "Observed (1 read)" = "#6baed6", 
    "Observed (>1 read)" = "#08519c" 
  )) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) 

obsv_am_large_clones_tissue

### proportion of SNVs in DNA-binding domain
### DNA-binding domain exons 5-8

tissue_filtering_SNV <- maf_masked_coding %>% 
  filter(Subject == "Patient") %>%
  filter(Tissue %in% all_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP") %>% 
  mutate(DBD = case_when(
    exon_number >= 5 & exon_number <= 8 ~ "DBD",
    TRUE ~ "non-DBD"
  )) %>%
  group_by(DBD) %>%
  summarise(mut_count = n(), .groups = "drop")

observed_dbd_propportion_tissues <- tissue_filtering_SNV %>%
  summarise(
    DBD_prop = mut_count[DBD == "DBD"] / sum(mut_count)
  )

tp53_total_length <- annotations %>% filter(gene_name == "TP53") %>% 
  filter(!is.na(exon_number)) %>% 
  summarise(cds_length = sum(width)) %>% 
  pull(cds_length) 

tp53_DBD_length <- annotations %>% filter(gene_name == "TP53")%>%
  filter(exon_number >= 5 & exon_number <= 8 ) %>% 
  summarise(dbd_length = sum(width)) %>% pull(dbd_length)

tp53_dbd_proportion <- tp53_DBD_length/tp53_total_length

dbd_comparison <- tibble(
  Category = c("Observed", "Expected"),
  Proportion = c(observed_dbd_propportion_tissues$DBD_prop, tp53_dbd_proportion)
)

dbd_prop_tissues <- ggplot(dbd_comparison, aes(x = Category, y = Proportion, fill = Category)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = sprintf("%.2f", Proportion)), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("grey60", "#3182bd")) +
  labs(
    x = NULL,
    y = "TP53 DBD proportion"
  ) +
  ylim(0,1) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )



################################################################################
### All blood
################################################################################
blood_filtering <- maf_masked_coding %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  #filter(Variant_Type=="SNP") %>% 
  #filter(Variant_Classification == "Missense_Mutation") %>%
  filter(!is.na(am_pathogenicity))

blood_filtering %>% print(width = Inf)
blood_filtering$group <- "Observed"
all_missense$group <- "All possible"
combine_df <- bind_rows(blood_filtering, all_missense)

lfs_subjects
ctx_subjects
am_groups <- combine_df %>% 
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>% 
  mutate(LFS = if_else(is.na(Subject), "All\npossible", LFS)) %>%
  mutate(CTX = if_else(Subject %in% ctx_subjects, "CTX", "non-CTX")) %>% 
  mutate(CTX = if_else(is.na(Subject), "", CTX)) %>%
  mutate(lfs_group = paste(LFS,CTX, sep='\n')) %>% 
  print(width = Inf)




obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = lfs_group, y = am_pathogenicity)) +
  # geom_violin(
  #   data = subset(am_groups, group == "All possible"),
  #   mapping = aes(x = lfs_group),
  #   fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
  #   position = position_nudge(x = 0.25)
  # ) +
  geom_boxplot() +
  # annotate("rect",
  #          xmin = 0.75, xmax = 1.25,
  #          ymin = -Inf, ymax = Inf,
  #          fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  #scale_color_manual(values = c("Not observed" = "grey60", "Observed" = "#3182bd")) +
  # scale_x_discrete(
  #   labels = c(
  #     "Not observed" = "Not\nobserved",
  #     "Observed (1 read)" = "Observed\n(1 read)",
  #     "Observed (>1 read)" = "Observed\n(>1 read)"
  #   )
  # ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups


obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = lfs_group, y = am_pathogenicity)) +
  geom_violin() +
  # annotate("rect",
  #          xmin = 0.75, xmax = 1.25,
  #          ymin = -Inf, ymax = Inf,
  #          fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  #scale_color_manual(values = c("Not observed" = "grey60", "Observed" = "#3182bd")) +
  # scale_x_discrete(
  #   labels = c(
  #     "Not observed" = "Not\nobserved",
  #     "Observed (1 read)" = "Observed\n(1 read)",
  #     "Observed (>1 read)" = "Observed\n(>1 read)"
  #   )
  # ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups

ggsave("results/am_blood_violin.png", obsv_am_blood_groups, width = 2.5, height = 1.5, units = "in", dpi = 300)

## plot all observed together
obsv_am_blood <- ggplot(combine_df, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c("Not observed" = "grey60", "Observed" = "#3182bd")) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood
## separate large and small observed clones
combine_df_large_clones <- combine_df %>%
  mutate(group = case_when(
    group == "Observed" & t_alt_count == 1 ~ "Observed (1 read)",
    group == "Observed" & t_alt_count > 1  ~ "Observed (>1 read)",
    TRUE ~ group
  ))
combine_df_large_clones$group <- factor(combine_df_large_clones$group,
                                        levels = c("Not observed", "Observed (1 read)", "Observed (>1 read)"))

obsv_am_large_clones_blood <- ggplot(combine_df_large_clones, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c(
    "No nobserved" = "grey60",
    "Observed (1 read)" = "#6baed6", 
    "Observed (>1 read)" = "#08519c" 
  )) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) 
obsv_am_large_clones_blood


### proportion of SNVs in DNA-binding domain
### DNA-binding domain exons 5-8

blood_filtering_SNV <- maf_masked_coding %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP") %>% 
  mutate(DBD = case_when(
    exon_number >= 5 & exon_number <= 8 ~ "DBD",
    TRUE ~ "non-DBD"
  )) %>%
  group_by(DBD) %>%
  summarise(mut_count = n(), .groups = "drop")

observed_dbd_propportion_blood <- blood_filtering_SNV %>%
  summarise(
    DBD_prop = mut_count[DBD == "DBD"] / sum(mut_count)
  )

tp53_total_length <- annotations %>% filter(gene_name == "TP53") %>% 
  filter(!is.na(exon_number)) %>% 
  summarise(cds_length = sum(width)) %>% 
  pull(cds_length) 

tp53_DBD_length <- annotations %>% filter(gene_name == "TP53")%>%
  filter(exon_number >= 5 & exon_number <= 8 ) %>% 
  summarise(dbd_length = sum(width)) %>% pull(dbd_length)

tp53_dbd_proportion <- tp53_DBD_length/tp53_total_length

dbd_comparison <- tibble(
  Category = c("Observed", "Expected"),
  Proportion = c(observed_dbd_propportion_blood$DBD_prop, tp53_dbd_proportion)
)


k <- blood_filtering_SNV %>% filter(DBD == "DBD") %>% pull(mut_count)    # observed mutations in DBD
n <- blood_filtering_SNV %>% summarise(total = sum(mut_count)) %>% pull(total)   # total observed mutations
p <- tp53_dbd_proportion # proportion of coding sites within the DBD
b_test <- binom.test(k, n, p = p)
b_test_p <- b_test$p.value
b_test_lower_ci <- b_test$conf.int[1]
b_test_upper_ci <- b_test$conf.int[2]


dbd_prop_blood <- ggplot(dbd_comparison, aes(x = Category, y = Proportion, fill = Category)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    data = dbd_comparison %>% filter(Category == "Observed") %>%
      mutate(ymin = b_test_lower_ci, ymax = b_test_upper_ci),
    aes(ymin = ymin, ymax = ymax),
    width = 0.1, size = 0.4
  ) +
  geom_text(aes(y = Proportion/2, label = sprintf("%.2f", Proportion)),
            size = 8/2.85) +
  scale_fill_manual(values = c("grey60", "#3182bd")) +
  labs(
    x = NULL,
    y = "TP53 DBD proportion"
  ) +
  ylim(0,1) +
  theme_classic(base_size = 8) +
  annotate(
    "text",
    x = 2, y = 0.95,
    label = paste0("italic(p) == ", signif(b_test_p, 2)),
    parse = TRUE,
    size = 8/2.85
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
dbd_prop_blood
##################
#### Save plots
##################

ggsave("results/am_blood.png", obsv_am_blood, width = 2.5, height = 1.5, units = "in", dpi = 300)
ggsave("results/am_large_clones_blood.png", obsv_am_large_clones_blood, width = 2.5, height = 1.5, units = "in", dpi = 300)
ggsave("results/dbd_proportion_blood.png", dbd_prop_blood, width = 2, height = 1.5, units = "in", dpi = 300)

ggsave("results/am_large_clones_tissues.png", obsv_am_large_clones_tissue, width = 2.5, height = 1.5, units = "in", dpi = 300)
ggsave("results/am_tissues.png", obsv_am_tissue, width = 2.5, height = 1.5, units = "in", dpi = 300)
ggsave("results/dbd_proportion_tissues.png", dbd_prop_tissues, width = 2, height = 1.5, units = "in", dpi = 300)











##################
#### compare among LFS indivudals
##################
maf_masked_coding %>% filter(Subject == "Family member C") %>% filter(Hugo_Symbol == "TP53") %>% print(width = Inf)
tissue_filtering <- maf_masked_coding %>% 
  filter(Subject %in% c("Patient", "Family member A","Family member B", "Family member C")) %>%
  #filter(Tissue %in% all_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  #filter(!is.na(exon_number)) %>%
  #filter(Variant_Type=="SNP") %>% 
  #filter(Variant_Classification == "Missense_Mutation") %>%
  filter(!is.na(am_pathogenicity))


tissue_filtering <- tissue_filtering %>%
  mutate(group = case_when(Subject == "Patient" ~ "LFS01", 
                           Subject == "Family member A" & Tissue == "PBMC" ~ "LFS02 blood",
                           Subject == "Family member A" & Tissue == "Urine cells" ~ "LFS02 urine", 
                           Subject == "Family member B" & Tissue == "PBMC" ~ "REL01 blood",
                           Subject == "Family member B" & Tissue == "Urine cells" ~ "REL01 urine",
                           Subject == "Family member C" & Tissue == "PBMC" ~ "LFS03 blood",
                           Subject == "Family member C" & Tissue == "Urine cells" ~ "LFS03 urine"))

tissue_filtering %>% print(width = Inf, n=Inf)

## plot all observed together

obsv_am_tissue <- ggplot(tissue_filtering, aes(x = group, y = am_pathogenicity, color = Variant_Classification)) +
  #geom_boxplot() +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +

  #scale_color_manual(values = c("LFS01" = "grey60", "LFS02 urine" = "#3182bd")) +
  
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    #legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_tissue


wilcox.test(
  am_pathogenicity ~ group,
  data = subset(tissue_filtering, group %in% c("LFS01", "LFS02 urine"))
)

## separate large and small observed clones
combine_df_large_clones <- combine_df %>%
  mutate(group = case_when(
    group == "Observed" & t_alt_count == 1 ~ "Observed (1 read)",
    group == "Observed" & t_alt_count > 1  ~ "Observed (>1 read)",
    TRUE ~ group
  ))
combine_df_large_clones$group <- factor(combine_df_large_clones$group,
                                        levels = c("Not observed", "Observed (1 read)", "Observed (>1 read)"))

obsv_am_large_clones_tissue <- ggplot(combine_df_large_clones, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  annotate("rect",
           xmin = 0.75, xmax = 1.25,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c(
    "Not observed" = "grey60",
    "Observed (1 read)" = "#6baed6", 
    "Observed (>1 read)" = "#08519c" 
  )) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "Not\nobserved",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) 

obsv_am_large_clones_tissue

### proportion of SNVs in DNA-binding domain
### DNA-binding domain exons 5-8

tissue_filtering_SNV <- maf_masked_coding %>% 
  filter(Subject == "Patient") %>%
  filter(Tissue %in% all_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP") %>% 
  mutate(DBD = case_when(
    exon_number >= 5 & exon_number <= 8 ~ "DBD",
    TRUE ~ "non-DBD"
  )) %>%
  group_by(DBD) %>%
  summarise(mut_count = n(), .groups = "drop")

observed_dbd_propportion_tissues <- tissue_filtering_SNV %>%
  summarise(
    DBD_prop = mut_count[DBD == "DBD"] / sum(mut_count)
  )

tp53_total_length <- annotations %>% filter(gene_name == "TP53") %>% 
  filter(!is.na(exon_number)) %>% 
  summarise(cds_length = sum(width)) %>% 
  pull(cds_length) 

tp53_DBD_length <- annotations %>% filter(gene_name == "TP53")%>%
  filter(exon_number >= 5 & exon_number <= 8 ) %>% 
  summarise(dbd_length = sum(width)) %>% pull(dbd_length)

tp53_dbd_proportion <- tp53_DBD_length/tp53_total_length

dbd_comparison <- tibble(
  Category = c("Observed", "Expected"),
  Proportion = c(observed_dbd_propportion_tissues$DBD_prop, tp53_dbd_proportion)
)

dbd_prop_tissues <- ggplot(dbd_comparison, aes(x = Category, y = Proportion, fill = Category)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = sprintf("%.2f", Proportion)), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("grey60", "#3182bd")) +
  labs(
    x = NULL,
    y = "TP53 DBD proportion"
  ) +
  ylim(0,1) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )

dbd_prop_tissues




#######################
######### All genes alphamissense supp figure
########################


all_genes <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-11-25-alphamissense/all_genes_am.tsv")

all_genes %>% print(width = Inf)
#write_delim(maf_masked_coding, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-11-25-alphamissense/maf_masked_coding.tsv")


blood_filtering <- all_genes %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  #filter(Variant_Type=="SNP") %>% 
  #filter(Variant_Classification == "Missense_Mutation") %>%
  filter(!is.na(am_pathogenicity.x) | !is.na(am_pathogenicity.y)) %>%
  mutate(am_pathogenicity.y = if_else(is.na(am_pathogenicity.y), am_pathogenicity.x, am_pathogenicity.y))


blood_filtering %>% print(width = Inf)
blood_filtering$group <- "Observed"
all_missense$group <- "All possible"
combine_df <- bind_rows(blood_filtering, all_missense)

lfs_subjects
ctx_subjects
am_groups <- combine_df %>% 
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>% 
  mutate(LFS = if_else(is.na(Subject), "All\npossible", LFS)) %>%
  mutate(CTX = if_else(Subject %in% ctx_subjects, "CTX", "no-CTX")) %>% 
  mutate(CTX = if_else(is.na(Subject), "", CTX)) %>%
  mutate(lfs_group = paste(LFS,CTX, sep='\n')) %>% 
  print(width = Inf)



am_groups <- am_groups %>% filter(gene_name %in% c("ASXL1", "BRINP3", "DNMT3A", "GATA2", "RAD21", "TET2", "TP53"))

mw_results <- am_groups %>%
  filter(group != "All possible") %>%        # same filtering as your plot
  filter(gene_name %in% c("ASXL1", "BRINP3", "DNMT3A", "GATA2", "RAD21", "TET2", "TP53")) %>%
  group_by(gene_name) %>%
  summarise(
    test = list(
      wilcox.test(
        am_pathogenicity.y ~ LFS, exact = FALSE
      )
    )
  ) %>%
  mutate(tidy = map(test, tidy)) %>%
  unnest(tidy) %>%
  dplyr::select(gene_name, statistic, p.value, method)
mw_results

obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = LFS, y = am_pathogenicity.y, color = LFS)) +
  facet_wrap(~ gene_name, nrow = 1) +
  # geom_violin(
  #   data = subset(am_groups, group == "All possible"),
  #   mapping = aes(x = lfs_group),
  #   fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
  #   position = position_nudge(x = 0.25)
  # ) +
  geom_boxplot(outlier.shape = NA) +
  # annotate("rect",
  #          xmin = 0.75, xmax = 1.25,
  #          ymin = -Inf, ymax = Inf,
  #          fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +
  # scale_x_discrete(
  #   labels = c(
  #     "Not observed" = "Not\nobserved",
  #     "Observed (1 read)" = "Observed\n(1 read)",
  #     "Observed (>1 read)" = "Observed\n(>1 read)"
  #   )
  # ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups

ggsave("results/AlphaMissense_chip_supp.png", obsv_am_blood_groups, width = 6, height = 2, units = "in", dpi = 300)


pathogenic_proportion <- am_groups %>% mutate(pathogenic = if_else(am_pathogenicity.y > 0.564, "Pathogenic", "Non-Pathogenic")) %>%
  group_by(gene_name, LFS) %>%
  summarise(n_total = n(),
            n_path = sum(pathogenic == "Pathogenic"),
            n_nonpath = n_total-n_path,
            prop_path = n_path/n_total,
            .groups = "drop")
pathogenic_proportion

fishers_path <- pathogenic_proportion %>% 
  dplyr::select(gene_name, LFS, n_total, n_path, n_nonpath) %>%
  pivot_wider(
    names_from = LFS,
    values_from = c(n_total, n_path, n_nonpath),
    names_sep = "_"
  )



prop_plot <- ggplot(pathogenic_proportion, aes(x = LFS, y = prop_path, fill = LFS)) +
  facet_wrap(~ gene_name, nrow = 1) +
  geom_col(width = 0.7, color = "black") +
  scale_fill_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    y = "Pathogenic\nproportion",
    x = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2,2,2,2, "pt")
  )

prop_plot
ggsave("results/AlphaMissense_proportion_supp.png", prop_plot, width = 6, height = 2, units = "in", dpi = 300)


combined <- obsv_am_blood_groups / prop_plot +
  plot_layout(heights = c(2, 1))

combined
ggsave("results/AlphaMissense_dist_proportion_supp.png", combined, width = 6, height = 3, units = "in", dpi = 300)





########################## tp53

am_groups <- am_groups %>% filter(gene_name %in% c("TP53"))
am_groups %>% print(n=Inf)
unique(am_groups$lfs_group)
wilcox.test(am_pathogenicity.y ~ lfs_group, 
            data = am_groups %>% filter(lfs_group %in% c("LFS\nno-CTX", "non-LFS\nCTX")), exact = FALSE)

obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = lfs_group, y = am_pathogenicity.y)) +
  #facet_wrap(~ gene_name, nrow = 1) +
  # geom_violin(
  #   data = subset(am_groups, group == "All possible"),
  #   mapping = aes(x = lfs_group),
  #   fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
  #   position = position_nudge(x = 0.25)
  # ) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  # annotate("rect",
  #          xmin = 0.75, xmax = 1.25,
  #          ymin = -Inf, ymax = Inf,
  #          fill = "white", color = NA) +
  geom_jitter(width = 0.15, size = 1.5, stroke = 1.2, alpha = 1, aes(color = LFS, shape = CTX)) +
  scale_color_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +
  scale_shape_manual(values = c("CTX" = 16, "no-CTX" = 1)) +
  # scale_x_discrete(
  #   labels = c(
  #     "Not observed" = "Not\nobserved",
  #     "Observed (1 read)" = "Observed\n(1 read)",
  #     "Observed (>1 read)" = "Observed\n(>1 read)"
  #   )
  # ) +
  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups

ggsave("results/AlphaMissense_tp53_supp.png", obsv_am_blood_groups, width = 3, height = 2, units = "in", dpi = 300)

