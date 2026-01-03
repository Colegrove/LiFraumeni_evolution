### linear models


##### CHIP
#### mutation frequency

### with UW07 and CTx

mutFreq_prep_CHIP_encode 
mutFreq_prep_TP53_encode %>% filter(coding == "coding") %>% dplyr::select(coding, Subject, denominator_coding, age, n_muts, mutFreq, LFS_n, CTx_n, shape_group) %>% print()
test <- mutFreq_prep_TP53_encode %>% filter(coding == "coding") %>% dplyr::select(coding, Subject, denominator_coding, age, n_muts, mutFreq, LFS_n, CTx_n, shape_group) %>% print()
test %>% dplyr::select(-denominator_noncoding)

write.csv(mutFreq_prep_CHIP_encode, "/Users/huntc10/Desktop/CHIP_freq.csv")
write.csv(mutFreq_prep_TP53_encode, "/Users/huntc10/Desktop/TP53_freq.csv")
write.csv(mutBurden_prep_CHIP_encode, "/Users/huntc10/Desktop/CHIP_burden.csv")
write.csv(mutBurden_prep_TP53_encode, "/Users/huntc10/Desktop/TP53_burden.csv")


mutFreq_prep_CHIP


mutFreq_prep_CHIP_encode <- mutFreq_prep_CHIP %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  )
mutFreq_prep_CHIP_encode
model_freq_CHIP_all_samples_covariates <- lm(n_muts ~ age + denominator_coding + LFS_n + CTx_n, data = mutFreq_prep_CHIP_encode)
summary(model_freq)

### with UW07 and without CTx
model_freq_CHIP_all_samples_noCTx <- lm(n_muts ~ age + denominator_coding + LFS_n, data = mutFreq_prep_CHIP_encode)
summary(model_freq_CHIP_all_samples_noCTx)

### without UW07 and without CTx
mutFreq_prep_CHIP_encode_exclude <- mutFreq_prep_CHIP_encode %>% filter(Subject != "UW volunteer 7")
model_freq_exclude_noCTx <- lm(n_muts ~ age + denominator_coding + LFS_n, data = mutFreq_prep_CHIP_encode_exclude)
summary(model_freq_exclude_noCTx)
confint(model_freq_exclude_noCTx)

### plot frequency covariates
coef_freq_CHIP <- broom::tidy(model_freq_CHIP_all_samples_noCTx, conf.int = TRUE)
coef_freq_CHIP <- coef_freq_CHIP %>% filter(term != "(Intercept)")

coef_freq_CHIP <- coef_freq_CHIP %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))
#term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS","CTx_n" = "CTx")
term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS")
coef_freq_CHIP
term_plot <- ggplot(coef_freq_CHIP, aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-9, 33), breaks = seq(-10, 35, 5)) +
  #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), limits = c(-3, 12)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  labs(x = "Effect size", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8))
term_plot
ggsave("results/linear_model_frequency_CHIP_ms.png", term_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)




##### mutation burden

## with UW07 and CTx
mutBurden_prep_CHIP_encode <- mutBurden_prep_CHIP %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  ) %>% 
  filter(coding == "coding")
mutBurden_prep_CHIP_encode
model_burden_CHIP_all_samples_covariates <- lm(mutReads ~ age + denominator_coding + LFS_n + CTx_n, data = mutBurden_prep_CHIP_encode)
summary(model_burden_CHIP_all_samples_covariates)

## with UW07 without CTx
model_burden_CHIP_all_samples_noCTx <- lm(mutReads ~ age + denominator_coding + LFS_n, data = mutBurden_prep_CHIP_encode)
summary(model_burden_CHIP_all_samples_noCTx)

## without UW07 without CTx
mutBurden_prep_CHIP_exclude <- mutBurden_prep_CHIP_encode %>% filter(Subject != "UW volunteer 7")
model_burden_CHIP_exclude_noCTx <- lm(mutReads ~ age + denominator_coding + LFS_n, data = mutBurden_prep_CHIP_exclude)
summary(model_burden_CHIP_exclude_noCTx)


### plot burden covariates
coef_burden_CHIP <- broom::tidy(model_burden_CHIP_all_samples_noCTx, conf.int = TRUE)

coef_burden_CHIP <- coef_burden_CHIP %>% filter(term != "(Intercept)")

coef_burden_CHIP <- coef_burden_CHIP %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))
#term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS","CTx_n" = "CTx")
term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS")
coef_burden_CHIP
term_plot <- ggplot(coef_burden_CHIP, aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  #scale_x_continuous(limits = c(-2000, 1000), breaks = seq(-3000, 1000, 1000)) +
  scale_x_continuous(
    #trans = pseudo_log_trans(base = 10),
    limits = c(-5000, 300),
    #breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100),
    breaks = c(-1000, -100, -10, 0, 10, 100),
    labels = label_math()
  ) +
  #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), limits = c(-2000, 100)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  labs(x = "Effect size", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))
term_plot
ggsave("results/linear_model_burden_CHIP_ms.png", term_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)






##### TP53 only

#### mutation frequency
## with UW07 and CTx
mutFreq_prep_TP53_encode <- mutFreq_prep_TP53 %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  ) %>%
  filter(coding=="coding")

model_freq_TP53 <- lm(n_muts ~ age + denominator_coding + LFS_n + CTx_n, data = mutFreq_prep_TP53_encode)
summary(model_freq_TP53)

### with UW07 and without CTx
model_freq_TP53_all_samples_noCTx <- lm(n_muts ~ age + denominator_coding + LFS_n, data = mutFreq_prep_TP53_encode)
summary(model_freq_TP53_all_samples_noCTx)


model

### plot frequency covariates
coef_freq_TP53 <- broom::tidy(model_freq_TP53_all_samples_noCTx)
coef_freq_TP53 <- coef_freq_TP53 %>% filter(term != "(Intercept)")

coef_freq_TP53 <- coef_freq_TP53 %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))
#term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS","CTx_n" = "CTx")
term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS")
coef_freq_TP53
term_plot <- ggplot(coef_freq_TP53, aes(y = term, x = estimate, xmin = estimate - std.error, xmax = estimate + std.error)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), limits = c(-3, 12)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  labs(x = "Effect size", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8))
term_plot
ggsave("results/linear_model_frequency_TP53_ms.png", term_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)


##### mutation burden
mutBurden_prep_TP53_encode <- mutBurden_prep_TP53 %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  ) %>%
  filter(coding == "coding")

model_burden_TP53_all_samples_covariates <- lm(mutReads ~ age + denominator_coding + LFS_n + CTx_n, data = mutBurden_prep_TP53_encode)
summary(model_burden_TP53_all_samples_covariates)

## with UW07 without CTx
model_burden_TP53_all_samples_noCTx <- lm(mutReads ~ age + denominator_coding + LFS_n, data = mutBurden_prep_TP53_encode)
summary(model_burden_TP53_all_samples_noCTx)


### plot burden covariates
coef_burden_TP53 <- broom::tidy(model_burden_TP53_all_samples_noCTx)
coef_burden_TP53 <- coef_burden_TP53 %>% filter(term != "(Intercept)")

coef_burden_TP53 <- coef_burden_TP53 %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))
#term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS","CTx_n" = "CTx")
term_labels <- c("age" = "Age","denominator_coding" = "Depth","LFS_n" = "LFS")
coef_burden_TP53
term_plot <- ggplot(coef_burden_TP53, aes(y = term, x = estimate, xmin = estimate - std.error, xmax = estimate + std.error)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), limits = c(-3, 12)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  labs(x = "Effect size", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8))
term_plot
ggsave("results/linear_model_burden_TP53_ms.png", term_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)
