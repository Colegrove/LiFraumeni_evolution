#### Hunter Colegrove
#### 28 Jul 2025
#### dN/dS of chip genes using classic method


### input all possible mutations file from omega
### mutations in chip bed file regions

all_possible_muts_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/all_possible_sites_annotated.tsv"
all_possible_muts <- read_delim(all_possible_muts_path, delim="\t")

## calculate expected mutations
expected_all <- all_possible_muts %>%
  filter(IMPACT != "intron_variant") %>%
  group_by(GENE) %>%
  summarise(
    synonymous = sum(IMPACT == "synonymous"),
    nonsynonymous = sum(IMPACT == "missense" | IMPACT == "nonsense" | IMPACT == "splice_region" | IMPACT == "essential_splice"),
    .groups = "drop"
  ) %>% 
  filter(synonymous >0 | nonsynonymous > 0)

## look at samples within LFS by subject
observed_all_LFS <- filt_maf_CHIP %>%
  filter(
    Tumor_Sample_Barcode != "DevDNA1_S1.1",
    InBed == TRUE,
    Tissue %in% c("PBMC", "Buffy coat"),
    Subject %in% c("Patient", "Family member A", "Family member C"),
  ) %>%
  mutate(Group = Subject,
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


# calculate dn/ds
dnds_all_LFS <- expected_all %>%
  inner_join(observed_all_LFS, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    pseudo_count_syn = (synonymous/(synonymous + nonsynonymous)),
    pseudo_count_non = (nonsynonymous/(synonymous + nonsynonymous)),
    dnds = ((nonsyn_obs + pseudo_count_non)/ nonsynonymous) / ((syn_obs + pseudo_count_syn)/ synonymous)
  )


## remove cases where no synonymous and no nonsynonymous
dnds_all_LFS <- dnds_all_LFS %>%
  group_by(GENE) %>%
  filter(
    all(syn_obs >0 | nonsyn_obs > 0),
    n_distinct(Group) == 3 
  ) %>%
  ungroup()

## filter genes based on LFS/non-LFS grouping analyzed
dnds_all_LFS <- dnds_all_LFS %>%
  filter(GENE %in% c("ASXL1","BRINP3","CEBPA","DNMT3A","EZH2","GATA2",
                     "HNRNPK","RAD21","RUNX1", "SMC1A","TET2","TP53","WT1"))

dnds_all_LFS %>% print(n=Inf)
set.seed(1)
n_reps <- 200
df <- dnds_all_LFS %>%
  select(GENE, Group, synonymous, nonsynonymous, syn_obs, nonsyn_obs,
         pseudo_count_syn, pseudo_count_non)

boot_list <- vector("list", nrow(df))

for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  n <- r$syn_obs + r$nonsyn_obs
  
  labels <- c(rep("syn", r$syn_obs), rep("non", r$nonsyn_obs))
  
  draws <- replicate(n_reps, sample(labels, size = n, replace = TRUE))
  if (is.vector(draws)) draws <- matrix(draws, nrow = n)
  
  syn_b <- colSums(draws == "syn")
  non_b <- n - syn_b
  
  syn_adj <- syn_b + r$pseudo_count_syn
  non_syn_adj <- non_b + r$pseudo_count_non
  non_sites <- r$nonsynonymous
  syn_sites <- r$synonymous
  
  dnds <- ((non_syn_adj / non_sites) /
             (syn_adj / syn_sites))
  
  boot_list[[i]] <- tibble(
    GENE = r$GENE,
    Group = r$Group,
    draw = seq_len(n_reps),
    syn_b = syn_b,
    non_b = non_b,
    dnds = dnds
  )
}

boot_res <- bind_rows(boot_list)

boot_iqr <- boot_res %>%
  group_by(GENE, Group) %>%
  summarise(
    dnds_q25 = quantile(dnds, 0.25, na.rm = TRUE),
    dnds_med = median(dnds, na.rm = TRUE),
    dnds_q75 = quantile(dnds, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


dnds_all_LFS <- dnds_all_LFS %>%
  mutate(total_mut = syn_obs + nonsyn_obs) %>%
  mutate(label_fraction = paste0(nonsyn_obs, "\n―\n", syn_obs))

group_colors <- c(
  "Family member A" = "#DDDDDD",
  "Family member C"     = "#AAAAAA",
  "Patient" = "#888888"
)

dnds_all_LFS <- dnds_all_LFS %>% left_join(boot_iqr, by=c("GENE", "Group"))


dnds_all_LFS_TP53 <- dnds_all_LFS %>% filter(GENE == "TP53")

patient_labels <- c("Family member A" = "LFS02", "Family member C" = "LFS03", "Patient" = "LFS01")
pd <- position_dodge(width = 1)
dnds_classic <- ggplot(dnds_all_LFS_TP53, aes(x = Group, y = dnds)) +
  geom_col(position = position_dodge(width = 1), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_errorbar(data = boot_iqr %>% filter(GENE == "TP53"),
                aes(x = Group, ymin = dnds_q25, ymax = dnds_q75),
                position = pd, width = 0.2, inherit.aes = FALSE, na.rm = TRUE) +
  geom_text(
    aes(y = 0, label = label_fraction),
    position = position_dodge(width = 1),
    vjust = +1.1,
    size = 3, 
    lineheight = 0.5
  ) +
  scale_y_continuous(limits = c(-2, 21)) +
  scale_fill_manual(values = group_colors, labels = patient_labels) +
  scale_x_discrete(labels = patient_labels) +
  labs(x = "Gene", y = "dN/dS (TP53)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.title = element_blank(), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size =8),
    axis.text.y = element_text(size=8),
    legend.position = c(0.5,0.8)
  ) 
  #coord_cartesian(ylim = c(-0.2, 30))
dnds_classic


ggsave("results/dnds_naive_LFS.png", dnds_classic, width = 3, height = 3, units = "in", dpi = 300)
