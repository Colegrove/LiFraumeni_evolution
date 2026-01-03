#### Hunter Colegrove
#### 28 Jul 2025



### input all possible mutations file from omega
### mutations in chip bed file regions

all_possible_muts_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/all_possible_sites_annotated.tsv"
all_possible_muts <- read_delim(all_possible_muts_path, delim="\t")


############## ALL genes ######################

expected_all <- all_possible_muts %>%
  filter(IMPACT != "intron_variant") %>%
  group_by(GENE) %>%
  summarise(
    synonymous = sum(IMPACT == "synonymous"),
    nonsynonymous = sum(IMPACT == "missense" | IMPACT == "nonsense" | IMPACT == "splice_region" | IMPACT == "essential_splice"),
    .groups = "drop"
  ) %>% 
  filter(synonymous >0 | nonsynonymous > 0)

expected_all %>% print(n = Inf)

observed_all <- filt_maf_CHIP %>%
  filter(Tumor_Sample_Barcode != "DevDNA1_S1.1", InBed == TRUE) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  mutate(type = case_when(
    Variant_Classification %in% c("Silent") ~ "synonymous",
    Variant_Classification %in% c(
      "Missense_Mutation", "Nonsense_Mutation", 
      "Frame_Shift_Del", "Frame_Shift_Ins", 
      "In_Frame_Ins", "In_Frame_Del", 
      "Nonstop_Mutation", "Splice_Site", "Splice_Region"
    ) ~ "nonsynonymous",
    TRUE ~ "exclude"
  )) %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    syn_obs = sum(type == "synonymous"),
    nonsyn_obs = sum(type == "nonsynonymous"),
    .groups = "drop"
  )

dnds_all <- expected_all %>%
  inner_join(observed_all, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    dnds = (nonsyn_obs / nonsynonymous) / (syn_obs / synonymous)
  )


dnds_all <- dnds_all %>% filter(dnds >0 & syn_obs > 0) %>% print(n = Inf)
dnds_classic <- ggplot(dnds_all, aes(x = GENE, y = dnds)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 8)) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dnds_classic


############### group by LFS
observed_all_LFS <- filt_maf_CHIP %>%
  filter(
    Tumor_Sample_Barcode != "DevDNA1_S1.1",
    InBed == TRUE,
    Tissue %in% c("PBMC", "Buffy coat")
  ) %>%
  mutate(
    # Create new grouping variable
    Group = if_else(
      Subject %in% c("Patient", "Family member A", "Family member C"),
      "LFS",  # LFS group
      "non-LFS"   # Non-LFS group
    ),
    # Classify variants
    type = case_when(
      Variant_Classification %in% c("Silent") ~ "synonymous",
      Variant_Classification %in% c(
        "Missense_Mutation", "Nonsense_Mutation", 
        "Frame_Shift_Del", "Frame_Shift_Ins", 
        "In_Frame_Ins", "In_Frame_Del", 
        "Nonstop_Mutation", "Splice_Site", "Splice_Region"
      ) ~ "nonsynonymous",
      TRUE ~ "exclude"
    )
  ) %>%
  group_by(Hugo_Symbol, Group) %>%
  summarise(
    syn_obs = sum(type == "synonymous"),
    nonsyn_obs = sum(type == "nonsynonymous"),
    .groups = "drop"
  )
observed_all_LFS %>% print(n=Inf)


dnds_all_LFS <- expected_all %>%
  inner_join(observed_all_LFS, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    dnds = (nonsyn_obs / nonsynonymous) / (syn_obs / synonymous)
  )


dnds_all_LFS <- dnds_all_LFS %>%
  group_by(GENE) %>%
  filter(all(syn_obs > 0 & dnds > 0)) %>%
  ungroup() %>%
  print(n = Inf)

dnds_all_LFS <- dnds_all_LFS %>%
  mutate(total_mut = syn_obs + nonsyn_obs) %>%
  mutate(label_fraction = paste0(nonsyn_obs, "\n―\n", syn_obs))

dnds_classic <- ggplot(dnds_all_LFS, aes(x = GENE, y = dnds, fill = Group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_text(
    aes(y = 0, label = label_fraction),
    position = position_dodge(width = 0.9),
    vjust = +1.1,
    size = 3, 
    lineheight = 0.5
  ) +
  scale_y_continuous(limits = c(-0.1, 5)) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(), 
    axis.title.x = element_blank(),
    legend.position = c(0.8,0.8)
  )
dnds_classic


ggsave("results/dnds_naive.png", dnds_classic, width = 6, height = 5, units = "in", dpi = 300)




############### group by LFS
observed_all_LFS <- filt_maf_CHIP %>%
  filter(
    Tumor_Sample_Barcode != "DevDNA1_S1.1",
    InBed == TRUE,
    Tissue %in% c("PBMC", "Buffy coat")
  ) %>%
  mutate(
    # Create new grouping variable
    Group = if_else(
      Subject %in% c("Patient", "Family member A", "Family member C"),
      "LFS",  # LFS group
      "non-LFS"   # Non-LFS group
    ),
    # Classify variants
    type = case_when(
      Variant_Classification %in% c("Silent") ~ "synonymous",
      Variant_Classification %in% c(
        "Missense_Mutation", "Nonsense_Mutation", 
        "Frame_Shift_Del", "Frame_Shift_Ins", 
        "In_Frame_Ins", "In_Frame_Del", 
        "Nonstop_Mutation", "Splice_Site", "Splice_Region"
      ) ~ "nonsynonymous",
      TRUE ~ "exclude"
    )
  ) %>%
  group_by(Hugo_Symbol, Group) %>%
  summarise(
    syn_obs = sum(type == "synonymous"),
    nonsyn_obs = sum(type == "nonsynonymous"),
    .groups = "drop"
  )
observed_all_LFS %>% print(n=Inf)


dnds_all_LFS <- expected_all %>%
  inner_join(observed_all_LFS, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    dnds = (nonsyn_obs / nonsynonymous) / (syn_obs / synonymous)
  )


dnds_all_LFS <- dnds_all_LFS %>%
  group_by(GENE) %>%
  filter(all(syn_obs > 0 & dnds > 0)) %>%   # ensures all rows in group have syn_obs > 0
  ungroup() %>%
  print(n = Inf)

dnds_classic <- ggplot(dnds_all_LFS, aes(x = GENE, y = dnds, fill = Group)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 8)) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dnds_classic




############## Tissues in patient ######################



expected_all <- all_possible_muts %>%
  filter(IMPACT != "intron_variant") %>%
  group_by(GENE) %>%
  summarise(
    synonymous = sum(IMPACT == "synonymous"),
    nonsynonymous = sum(IMPACT == "missense" | IMPACT == "nonsense" | IMPACT == "splice_region" | IMPACT == "essential_splice"),
    .groups = "drop"
  ) %>% 
  filter(synonymous >0 | nonsynonymous > 0)

expected_all %>% print(n = Inf)

unique(filt_maf$Subject)
observed_all_patient <- filt_maf %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>%
  mutate(type = case_when(
    Variant_Classification %in% c("Silent") ~ "synonymous",
    Variant_Classification %in% c(
      "Missense_Mutation", "Nonsense_Mutation", 
      "Frame_Shift_Del", "Frame_Shift_Ins", 
      "In_Frame_Ins", "In_Frame_Del", 
      "Nonstop_Mutation", "Splice_Site", "Splice_Region"
    ) ~ "nonsynonymous",
    TRUE ~ "exclude"
  )) %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    syn_obs = sum(type == "synonymous"),
    nonsyn_obs = sum(type == "nonsynonymous"),
    .groups = "drop"
  )

dnds_all_patient <- expected_all %>%
  inner_join(observed_all, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    dnds = (nonsyn_obs / nonsynonymous) / (syn_obs / synonymous)
  )


dnds_all <- dnds_all %>% filter(dnds >0 & syn_obs > 0) %>% print(n = Inf)
dnds_classic <- ggplot(dnds_all, aes(x = GENE, y = dnds)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 8)) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dnds_classic






