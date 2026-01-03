### Run AMSD


devtools::install_github("sfhart33/mutspecdist")
library(mutspecdist)


sbs_matrix <- read_delim("/Users/huntc10/Downloads/li_fraumeni.SBS96-2.all")
sbs_matrix

sbs_matrix_flipped <- sbs_matrix %>%
  pivot_longer(cols = -MutationType,
               names_to = "Sample",
               values_to = "Count") %>%
  pivot_wider(names_from = MutationType, 
              values_from = Count)

names(sbs_matrix_flipped) <- gsub("([A-Z])\\[([A-Z])>([A-Z])\\]([A-Z])", "\\1\\2\\4>\\1\\3\\4", names(sbs_matrix_flipped))

matrix_spectra <- sbs_matrix_flipped %>%
  rowwise() %>%  
  mutate(across(-Sample, ~ .x / sum(c_across(-Sample)))) %>% 
  ungroup()

# mean
spectra_mat <- matrix_spectra %>%
  dplyr::select(-Sample) %>%
  as.matrix()
rownames(spectra_mat) <- matrix_spectra$Sample

# sum
spectra_mat <- sbs_matrix_flipped %>%
  dplyr::select(-Sample) %>%
  as.matrix()
rownames(spectra_mat) <- matrix_spectra$Sample

################# LFS vs non-LFS
subset_mat_nonLFS <- spectra_mat[c("DNA03963_S3.1", "DNA03968_S5.1", 
                                   "DNA03969_S6.1", "DNA03970_S7.1", 
                                   "DNA03971_S8.1", "DNA03972_S9.1",
                                   "DNA03973_S10.1", "DNA03974_S11.1"), ]
subset_mat_LFS <- spectra_mat[c("DNA03975_S12.1", "DNA03962_S2.1", 
                                "DNA03964_S4.1"), ]

subset_mat_nonLFS_no_chemo <- spectra_mat[c( 
                                            "DNA03969_S6.1", "DNA03970_S7.1", 
                                            
                                            "DNA03973_S10.1"), ]
subset_mat_LFS_no_chemo <- spectra_mat[c("DNA03962_S2.1", 
                                         "DNA03964_S4.1"), ]


col_sums <- colSums(subset_mat_LFS_no_chemo)
col_sums_df <- data.frame(
  mutation = names(col_sums),
  count = as.numeric(col_sums)
) %>% mutate(group = "LFS")
col_sums_df$prop <- col_sums_df$count / sum(col_sums_df$count)
sum(col_sums_df$prop)

col_sums2 <- colSums(subset_mat_nonLFS_no_chemo)
col_sums_df2 <- data.frame(
  mutation = names(col_sums2),
  count = as.numeric(col_sums2)
) %>% mutate(group = "non-LFS")
col_sums_df2$prop <- col_sums_df2$count / sum(col_sums_df2$count)
sum(col_sums_df$prop)

col_sums_df <- rbind(col_sums_df, col_sums_df2)
col_sums_df
lfs_96 <- ggplot(col_sums_df, aes(x = mutation, y = prop, fill = group)) +
  geom_col(position = 'dodge') +
  labs(
    x = "Mutation type",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.title.x = element_blank(), 
    legend.position = c(0.5,0.9)
  )

ggsave("results/LFS_96_comparison.png", lfs_96, width = 8, height = 4, units = "in", dpi=300)




col_sums_LFS <- colSums(subset_mat_LFS_no_chemo)
col_sums_df_LFS <- data.frame(
  mutation = names(col_sums_LFS),
  count    = as.numeric(col_sums_LFS)
) %>%
  mutate(
    group = "LFS",
    prop  = count / sum(count)
  )

# non-LFS
col_sums_nonLFS <- colSums(subset_mat_nonLFS_no_chemo)
col_sums_df_nonLFS <- data.frame(
  mutation = names(col_sums_nonLFS),
  count    = as.numeric(col_sums_nonLFS)
) %>%
  mutate(
    group = "non-LFS",
    prop  = count / sum(count)
  )

# combine
col_sums_df <- bind_rows(col_sums_df_LFS, col_sums_df_nonLFS)

## 2. Add COSMIC-style substitution info & ordering ----

# helper to normalize substitutions to pyrimidine context (C or T)
normalize_sub <- function(ref, alt) {
  rc <- c(A = "T", C = "G", G = "C", T = "A")  # reverse complement map
  if (ref %in% c("C", "T")) {
    paste0(ref, ">", alt)
  } else {
    paste0(rc[ref], ">", rc[alt])
  }
}

subs_order <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

col_sums_df <- col_sums_df %>%
  mutate(
    context_ref = str_sub(mutation, 1, 3),   # e.g. "ACA"
    context_alt = str_sub(mutation, 5, 7),   # e.g. "AAA"
    ref         = str_sub(context_ref, 2, 2),
    alt         = str_sub(context_alt, 2, 2)
  ) %>%
  rowwise() %>%
  mutate(subs_norm = normalize_sub(ref, alt)) %>%
  ungroup() %>%
  mutate(
    subs_norm   = factor(subs_norm, levels = subs_order)
  ) %>%
  arrange(subs_norm, context_ref, context_alt) %>%
  mutate(
    # set factor levels in the desired plotted order
    mutation = factor(mutation, levels = unique(mutation))
  )

## 3. Plot: 96-channel style, ordered by COSMIC substitution class ----

lfs_96 <- ggplot(col_sums_df, aes(x = mutation, y = prop, fill = group)) +
  geom_col(position = "dodge") +
  labs(
    x = "Mutation type",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.title.x = element_blank(),
    legend.position = c(0.5, 0.9)
  )

lfs_96















exp_v_con <- amsd(subset_mat_LFS_no_chemo, subset_mat_nonLFS_no_chemo, mean_or_sum = "mean", seed = 123) 
#exp_v_con <- amsd(subset_mat_LFS, subset_mat_nonLFS, mean_or_sum = "sum", seed = 123) 
plot_amsd_histogram(exp_v_con)
exp_v_con$cosine
#exp_v_con

exp_v_con <- amsd(subset_mat_LFS, subset_mat_nonLFS, mean_or_sum = "mean", seed = 123) 
#exp_v_con <- amsd(subset_mat_LFS, subset_mat_nonLFS, mean_or_sum = "sum", seed = 123) 
plot_amsd_histogram(exp_v_con)
exp_v_con$cosine
#exp_v_con
