### multiple regression analysis tp53 only

## generate mutFreq_combined from mutation_frequencies_blood_CHIP_ms_2BC_v2.R

mutFreq_combined_TP53 <- mutFreq_combined %>%
  filter(Hugo_Symbol == "TP53")

# helper to encode and fit
mutFreq_prep_TP53_encode <- mutFreq_combined_TP53 %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  )

mutFreq_prep_TP53_encode
##########################################
# TP53 Frequency
##########################################

mutFreq_subject <- mutFreq_prep_TP53_encode %>%
  group_by(Subject, age, LFS, CTx, LFS_n, CTx_n) %>%
  summarise(
    n_muts = sum(n_muts, na.rm = TRUE),
    mutReads = sum(mutReads, na.rm = TRUE),
    denominator_coding = sum(denominator_coding, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(age_decades = age / 10,
         depth_scaled = denominator_coding / 10000000) %>%
  mutate(MF = n_muts/denominator_coding,
         MB = mutReads/denominator_coding,
         MF_scale = n_muts/depth_scaled,
         MB_scale = mutReads/depth_scaled)

mutFreq_subject
# # With UW07 and no CTx
# model_freq_TP53_all_samples_noCTx <- lm(
#   #n_muts ~ age_decades + depth_scaled + LFS_n,
#   MF_scale ~ age_decades + depth_scaled + LFS_n,
#   data = mutFreq_subject
# )
# summary(model_freq_TP53_all_samples_noCTx)

## with ctx
model_freq_TP53_all_samples_CTx <- lm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ age_decades + LFS_n + CTx_n,
  data = mutFreq_subject
)
summary(model_freq_TP53_all_samples_CTx)

coef_freq_TP53 <- tidy(model_freq_TP53_all_samples_CTx, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))

term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
coef_freq_TP53
TP53_freq_plot <- ggplot(coef_freq_TP53,
                    aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-3, 9), breaks = seq(-3, 9, 3)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  labs(x = "Effect size", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

TP53_freq_plot

ggsave("results/covariates_frequency_TP53_ms.png", TP53_freq_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)




##########################################
# TP53 Burden
##########################################

# # With UW07 and no CTx
# model_burden_TP53_all_samples_noCTx <- lm(
#   #mutReads ~ age_decades + depth_scaled + LFS_n,
#   MB_scale ~ age_decades + depth_scaled + LFS_n,
#   data = mutFreq_subject
# )
# summary(model_burden_TP53_all_samples_noCTx)

# with ctx
model_burden_TP53_all_samples_CTx <- lm(
  #mutReads ~ age_decades + depth_scaled + LFS_n,
  #MB_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MB_scale ~ age_decades + LFS_n + CTx_n,
  data = mutFreq_subject
)
summary(model_burden_TP53_all_samples_CTx)


coef_burden_TP53 <- tidy(model_burden_TP53_all_samples_CTx, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))

term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
coef_burden_TP53
TP53_burden_plot <- ggplot(coef_burden_TP53,
                    aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_color_manual(
    name = "Significance",
    values = c("black"), 
    labels = c("*** p<0.0001", "** p<0.01", "* p<0.05")
  ) +
  scale_x_continuous(limits = c(-10, 20), breaks = seq(-10, 20, 10)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  labs(x = "Effect size", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
        legend.position = 'right')

TP53_burden_plot
ggsave("results/covariates_burden_TP53_ms.png", TP53_burden_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)

####### freq + burden plot
TP53_freq_plot <- TP53_freq_plot + theme(
  plot.margin = margin(1,1,1,1)
)
TP53_burden_plot <- TP53_burden_plot + theme(
  plot.margin = margin(1,1,1,1)
)
combined_plot <- TP53_freq_plot | TP53_burden_plot

ggsave("results/covariates_combined_TP53_ms.png", combined_plot, width = 3.5, height = 1.5, units = "in", dpi = 300)

