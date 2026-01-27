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
all_possible_muts_path <- "inputs/all_possible_sites_annotated.tsv.gz"
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
      "Not observed" = "All possible\nSNVs",
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


mask_rect_df <- data.frame(
  xmin = c(0.75, 1.75, 2.75),
  xmax = c(1.25, 2.25, 3.25),
  ymin = -Inf,
  ymax = Inf
)
combine_df_large_clones %>% print(width = Inf)
n_counts <- combine_df_large_clones %>% group_by(group) %>% summarise(n=n())
n_all <- n_counts %>% filter(group == "Not observed") %>% pull(n)
n_1 <- n_counts %>% filter(group == "Observed (1 read)") %>% pull(n)
n_large <- n_counts %>% filter(group == "Observed (>1 read)") %>% pull(n)

n_counts
n_all
n_1
n_large

obsv_am_large_clones_tissue <- ggplot(combine_df_large_clones, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Observed (1 read)"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Observed (>1 read)"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  # annotate("rect",
  #          xmin = 0.75, xmax = 1.25,
  #          ymin = -Inf, ymax = Inf,
  #          fill = "white", color = NA) +
  geom_rect(
    data = mask_rect_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white", color = NA
  ) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c(
    "Not observed" = "grey60",
    "Observed (1 read)" = "#6baed6", 
    "Observed (>1 read)" = "#08519c" 
  )) +
  scale_x_discrete(
    labels = c(
      "Not observed" = "All possible\nSNVs",
      "Observed (1 read)" = "Observed\n(1 read)",
      "Observed (>1 read)" = "Observed\n(>1 read)"
    )
  ) +
  scale_y_continuous(limit = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  #labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  labs(y = "Pathogenicity score", x = NULL) +
  scale_x_discrete(
    labels = c(
      "Not observed"       = paste0("All possible\nSNVs\nn = ", n_all),
      "Observed (1 read)"  = paste0("Observed\n(1 read)\nn = ", n_1),
      "Observed (>1 read)" = paste0("Observed\n(>1 read)\nn = ", n_large)
    )
  ) +
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

tissue_filtering_SNV

k <- tissue_filtering_SNV %>% filter(DBD == "DBD") %>% pull(mut_count)    # observed mutations in DBD
n <- tissue_filtering_SNV %>% summarise(total = sum(mut_count)) %>% pull(total)   # total observed mutations
p <- tp53_dbd_proportion # proportion of coding sites within the DBD
b_test <- binom.test(k, n, p = p)
b_test_p <- b_test$p.value
b_test_lower_ci <- b_test$conf.int[1]
b_test_upper_ci <- b_test$conf.int[2]




dbd_prop_tissues <- ggplot(dbd_comparison, aes(x = Category, y = Proportion, fill = Category)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    data = dbd_comparison %>% filter(Category == "Observed") %>%
      mutate(ymin = b_test_lower_ci, ymax = b_test_upper_ci),
    aes(ymin = ymin, ymax = ymax),
    width = 0.1, size = 0.4
  ) +
  geom_text(aes(y = Proportion/2, label = sprintf("%.2f", Proportion)),
            size = 8/2.85) +
  annotate(
    "text",
    x = 2, y = 0.95,
    label = paste0("italic(p) == ", signif(b_test_p, 2)),
    parse = TRUE,
    size = 8/2.85
  ) +
  scale_fill_manual(values = c("grey60", "#3182bd")) +
  scale_x_discrete(labels = function(x){ifelse(x == "Observed", paste0("Observed\nn = ", k), x)}) +
  labs(
    x = NULL,
    y = "*TP53* DBD proportion"
  ) +
  #ylim(0,1) +
  scale_y_continuous(limit = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_markdown(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 0, "pt")
  )
dbd_prop_tissues



combine_dbd_tissues <- dbd_prop_tissues + obsv_am_large_clones_tissue +
  plot_layout(widths = c(0.61, 1)) &
  theme(plot.margin = margin(2,2,2,1))
combine_dbd_tissues

#ggsave("results/dbd_snvs_combined_ms_4.png", combine_dbd_tissues, width = 3.75, height = 1.5, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_4/dbd_snvs_combined_ms_4.png", combine_dbd_tissues, width = 3.75, height = 1.5, units = "in", dpi = 300)


################################################################################
### All blood
################################################################################
blood_filtering <- maf_masked_coding %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP") %>% 
  filter(Variant_Classification == "Missense_Mutation") %>%
  filter(!is.na(am_pathogenicity))

blood_filtering %>% print(width = Inf)
blood_filtering$group <- "Observed"
all_missense$group <- "Not observed"
combine_df <- bind_rows(blood_filtering, all_missense)

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

n_counts <- combine_df_large_clones %>% group_by(group) %>% summarise(n=n())
n_all <- n_counts %>% filter(group == "Not observed") %>% pull(n)
n_1 <- n_counts %>% filter(group == "Observed (1 read)") %>% pull(n)
n_large <- n_counts %>% filter(group == "Observed (>1 read)") %>% pull(n)

obsv_am_large_clones_blood <- ggplot(combine_df_large_clones, aes(x = group, y = am_pathogenicity, color = group)) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Not observed"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Observed (1 read)"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  geom_violin(
    data = subset(combine_df_large_clones, group == "Observed (>1 read)"),
    mapping = aes(x = group),
    fill = "grey60", color = NA, alpha = 0.6, width = 0.8, trim = TRUE,
    position = position_nudge(x = 0.25)
  ) +
  geom_rect(
    data = mask_rect_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white", color = NA
  ) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c(
    "No nobserved" = "grey60",
    "Observed (1 read)" = "#6baed6", 
    "Observed (>1 read)" = "#08519c" 
  )) +
  scale_y_continuous(limit = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_x_discrete(
    labels = c(
      "Not observed"       = paste0("All possible\nSNVs\nn = ", n_all),
      "Observed (1 read)"  = paste0("Observed\n(1 read)\nn = ", n_1),
      "Observed (>1 read)" = paste0("Observed\n(>1 read)\nn = ", n_large)
    )
  ) +
  labs(y = "Pathogenicity score", x = NULL) +
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
    y = "*TP53* DBD proportion"
  ) +
  #ylim(0,1) +
  scale_y_continuous(limit = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme_classic(base_size = 8) +
  scale_x_discrete(labels = function(x){ifelse(x == "Observed", paste0("Observed\nn = ", k), x)}) +
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
    axis.title.y = element_markdown(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
dbd_prop_blood


##################
#### Save plots
##################


combine_dbd_blood <- dbd_prop_blood + obsv_am_large_clones_blood +
  plot_layout(widths = c(0.61, 1)) &
  theme(plot.margin = margin(2,2,2,1))
combine_dbd_blood

#ggsave("results/dbd_snvs_combined_ms_3.png", combine_dbd_blood, width = 3.75, height = 1.5, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/dbd_snvs_combined_ms_3.png", combine_dbd_blood, width = 3.75, height = 1.5, units = "in", dpi = 300)

