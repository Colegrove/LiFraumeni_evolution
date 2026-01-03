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

subset_mat_nonLFS_no_chemo <- spectra_mat[c("DNA03963_S3.1", "DNA03968_S5.1", 
                                   "DNA03969_S6.1", "DNA03970_S7.1", 
                                   "DNA03971_S8.1", "DNA03972_S9.1",
                                   "DNA03973_S10.1"), ]
subset_mat_LFS_no_chemo <- spectra_mat[c("DNA03962_S2.1", 
                                "DNA03964_S4.1"), ]

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



# ################### control example where we expect a difference
# subset_mat_SBS7 <- spectra_mat[c("DNA03969_S6.1", "DNA03970_S7.1", 
#                                 "DNA03963_S3.1", "DNA03962_S2.1", 
#                                 "DNA03964_S4.1"), ]
# subset_mat_noSBS7 <- spectra_mat[c("DNA03968_S5.1","DNA03971_S8.1","DNA03972_S9.1", 
#                                    "DNA03973_S10.1", "DNA03974_S11.1"), ]
# 
# 
# exp_v_con <- amsd(subset_mat_SBS7, subset_mat_noSBS7, mean_or_sum = "sum", seed = 123) 
# plot_amsd_histogram(exp_v_con)
# exp_v_con
# 
# 
# 
# maf_masked_coding %>% group_by(Tumor_Sample_Barcode) %>% filter(!inRepeatMask | is.na(inRepeatMask)) %>% summarise(count = n()) %>% print(n =Inf)
# col_sums <- colSums(sbs_matrix[ , -1])
# col_sums
