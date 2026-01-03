### multiple regression analysis all chip genes

library(broom)
library(patchwork)

## generate mutFreq_combined from mutation_frequencies_blood_CHIP_ms_2BC_v3.R
mutFreq_prep %>% print(width = Inf, n = Inf)
#mutFreq_combined
# helper to encode and fit
mutFreq_prep_CHIP_encode <- mutFreq_prep %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  )

mutFreq_prep_CHIP_encode
##########################################
# CHIP Frequency
##########################################

# mutFreq_subject <- mutFreq_prep_CHIP_encode %>%
#   group_by(Subject, age, LFS, CTx, LFS_n, CTx_n) %>%
#   summarise(
#     n_muts = sum(n_muts, na.rm = TRUE),
#     mutReads = sum(mutReads, na.rm = TRUE),
#     denominator_coding = sum(denominator_coding, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(age_decades = age / 10,
#          depth_scaled = denominator_coding / 10000000) %>%
#   mutate(MF = n_muts/denominator_coding,
#          MB = mutReads/denominator_coding,
#          MF_scale = n_muts/depth_scaled,
#          MB_scale = mutReads/depth_scaled)

mutFreq_subject <- mutFreq_prep_CHIP_encode %>%
  mutate(age_decades = age / 10,
         depth_scaled = denominator / 10000000) %>%
  mutate(MF_scale = n_muts/depth_scaled)


# With UW07 and no CTx
# model_freq_CHIP_all_samples_noCTx <- lm(
#   #n_muts ~ age_decades + depth_scaled + LFS_n,
#   #MF ~ age_decades + depth_scaled + LFS_n, 
#   MF_scale ~ age_decades + depth_scaled + LFS_n,
#   data = mutFreq_subject
# )
# summary(model_freq_CHIP_all_samples_noCTx)

#non_coding_set = 0 # non-coding-total
#non_coding_set = 1 # non-coding-MUT
non_coding_set = 2 # non-coding-CHIP

if(non_coding_set == 0){
  #lm_model <- lm(mutFreq ~ age, data = mutFreq_subject %>% filter(coding == "non-coding-total"))
  mutFreq_subject_non_coding = mutFreq_subject %>% filter(coding == "non-coding-total")
  xlimits = c(-1.2, 1.5)
  xbreaks = seq(-1, 1.5, 0.5)
  #ylabel = "Non-coding\nCHIP + MUT MF"
  file_out = "MF_multi_coding_non_coding_total_ms.png"
}
if(non_coding_set == 1){
  #lm_model <- lm(mutFreq ~ age, data = mutFreq_subject %>% filter(coding == "non-coding-MUT"))
  mutFreq_subject_non_coding = mutFreq_subject %>% filter(coding == "non-coding-MUT")
  xlimits = c(-1.6, 1.63)
  xbreaks = seq(-1.5, 1.5, 0.5)
  #ylabel = "MUT MF"
  file_out = "MF_multi_coding_non_coding_MUTonly_ms.png"
}
if(non_coding_set == 2){
  #lm_model <- lm(mutFreq ~ age, data = mutFreq_subject %>% filter(coding == "non-coding-CHIP"))
  mutFreq_subject_non_coding = mutFreq_subject %>% filter(coding == "non-coding-CHIP")
  xlimits = c(-0.5, 1.5)
  xbreaks = seq(-0.5, 1.5, 0.5)
  #ylabel = "non-coding\nCHIP MF"
  file_out = "MF_multi_coding_non_coding_CHIPonly_ms.png"
}

## coding
mutFreq_subject_coding <- mutFreq_subject %>%
  filter(coding == "coding")

mutFreq_subject_coding
model_freq_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_coding
)
summary(model_freq_coding)


model_freq_non_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ 1 + age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_non_coding
)
summary(model_freq_non_coding)
mutFreq_subject_non_coding

coef_freq_CHIP_coding <- tidy(model_freq_coding, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))


term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
CHIP_freq_plot_coding <- ggplot(coef_freq_CHIP_coding,
                    aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.2, size = 5, color = "black") +
  labs(x = expression("Effect size\n(mutations / "  * 10^7 * " bases)"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

CHIP_freq_plot_coding

ggsave("results/MF_multi_coding_ms.png", CHIP_freq_plot_coding, width = 1.5, height = 1.5, units = "in", dpi = 300)

#########
##non-coding

tidy(model_freq_non_coding)

coef_freq_CHIP_non_coding <- tidy(model_freq_non_coding, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))

term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")


coef_freq_CHIP_non_coding

CHIP_freq_plot_non_coding <- ggplot(coef_freq_CHIP_non_coding,
                                aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = xlimits, breaks = xbreaks) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.2, size = 5, color = "black") +
  labs(x = expression("Effect size\n(mutations / "  * 10^7 * " bases)"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

CHIP_freq_plot_non_coding

ggsave("results/MF_multi_noncoding_ms.png", CHIP_freq_plot_non_coding, width = 1.5, height = 1.5, units = "in", dpi = 300)




# ##########################################
# # CHIP Burden
# ##########################################
# 
# # With UW07 and no CTx
# model_burden_CHIP_all_samples_noCTx <- lm(
#   #mutReads ~ age_decades + depth_scaled + LFS_n,
#   #MB ~ age_decades + depth_scaled + LFS_n,
#   MB_scale ~ age_decades + depth_scaled + LFS_n,
#   data = mutFreq_subject
# )
# summary(model_burden_CHIP_all_samples_noCTx)
# 
# ## with ctx
# model_burden_CHIP_all_samples_CTx <- lm(
#   #mutReads ~ age_decades + depth_scaled + LFS_n,
#   #MB ~ age_decades + depth_scaled + LFS_n,
#   #MB_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
#   MB_scale ~ age_decades + LFS_n + CTx_n,
#   data = mutFreq_subject
# )
# summary(model_burden_CHIP_all_samples_CTx)
# 
# 
# coef_burden_CHIP <- tidy(model_burden_CHIP_all_samples_CTx, conf.int = TRUE) %>%
#   filter(term != "(Intercept)") %>%
#   mutate(sig = case_when(
#     p.value < 0.0001 ~ "***",
#     p.value < 0.01   ~ "**",
#     p.value < 0.05   ~ "*",
#     TRUE ~ ""
#   ))
# 
# term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
# coef_burden_CHIP
# CHIP_burden_plot <- ggplot(coef_burden_CHIP,
#                     aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
#   geom_point(color = "black", size = 1) +
#   geom_errorbarh(height = 0, color = "black") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_y_discrete(labels = term_labels) +
#   scale_color_manual(
#     name = "Significance",
#     values = c("black"), 
#     labels = c("*** p<0.0001", "** p<0.01", "* p<0.05")
#   ) +
#   scale_x_continuous(limits = c(-30, 50), breaks = seq(-30, 60, 15)) +
#   geom_text(aes(label = sig, x = estimate),
#             hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
#   labs(x = "Effect size", y = NULL) +
#   theme_classic(base_size = 8) +
#   theme(axis.text.y = element_markdown(size = 8),
#         axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
#         legend.position = 'right')
# 
# CHIP_burden_plot
# ggsave("results/covariates_burden_CHIP_ms.png", CHIP_burden_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)

####### freq + burden plot

CHIP_freq_plot_coding <- CHIP_freq_plot_coding + theme(
  plot.margin = margin(1,1,1,1)) + 
  coord_cartesian(clip = "off")
CHIP_freq_plot_non_coding <- CHIP_freq_plot_non_coding + theme(
  plot.margin = margin(1,1,1,1)) + 
  coord_cartesian(clip = "off")
# CHIP_burden_plot <- CHIP_burden_plot + theme(
#   plot.margin = margin(1,1,1,1)
# )


combined_plot <- CHIP_freq_plot_coding | CHIP_freq_plot_non_coding

combined_plot

ggsave(paste0("results/", file_out), combined_plot, width = 3.5, height = 1, units = "in", dpi = 300)
ggsave(paste0("/Users/huntc10/Library/CloudStorage/OneDrive-UW/Li-Fraumeni/Manuscript_figures/Fig_2/", file_out), combined_plot, width = 3.5, height = 1, units = "in", dpi = 300)


# save supp table

model_to_row <- function(model, gene, coding) {
  
  coding_tidy <- broom::tidy(model)
  
  coding_wide <- coding_tidy %>%
    dplyr::select(term, estimate, std.error, p.value) %>%
    dplyr::mutate(
      est_se = sprintf("%.4g (%.4g)", estimate, std.error),
      p      = p.value
    ) %>% 
    dplyr::select(term, est_se, p) %>% 
    dplyr::filter(term %in% c("(Intercept)", "age_decades", "LFS_n", "CTx_n")) %>% 
    tidyr::pivot_wider(
      names_from  = term,
      values_from = c(est_se, p),
      names_glue  = "{term}.{.value}"
    ) %>% 
    dplyr::mutate(
      model  = "MF ~ 1 + decades + LFS + CTx",
      gene   = gene,
      coding = coding,
      .before = 1
    ) %>% 
    dplyr::rename(
      `(Intercept).Estimate` = `(Intercept).est_se`,
      `decades.Estimate`     = `age_decades.est_se`,
      `LFS.Estimate`         = `LFS_n.est_se`,
      `CTx.Estimate`         = `CTx_n.est_se`,
      `decades.p`            = `age_decades.p`,
      `LFS.p`                = `LFS_n.p`,
      `CTx.p`                = `CTx_n.p`
    ) %>%
    dplyr::select(
      model, gene, coding,
      `(Intercept).Estimate`, `(Intercept).p`,
      `decades.Estimate`,     `decades.p`,
      `LFS.Estimate`,         `LFS.p`,
      `CTx.Estimate`,         `CTx.p`
    )
  
  coding_wide
}
row_coding <- model_to_row(model_freq_coding,      gene = "CHIP", coding = "coding")
row_noncoding <- model_to_row(model_freq_non_coding,   gene = "CHIP", coding = "non-coding")
model_table <- bind_rows(row_coding, row_noncoding)

write.table(
  model_table,
  file = "results/MF_multivariate_CHIP_summary.csv",
  sep = ",",
  row.names = FALSE
)

