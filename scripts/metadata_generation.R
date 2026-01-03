
## send metadata to henry



original_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
metadata <- read_delim(original_path)
metadata %>% print(n=Inf)


sample_list <- c(
  "DNA03962_S2.1", "DNA03963_S3.1", "DNA03964_S4.1", "DNA03968_S5.1",
  "DNA03969_S6.1", "DNA03970_S7.1", "DNA03971_S8.1", "DNA03972_S9.1",
  "DNA03973_S10.1", "DNA03974_S11.1", "DNA03975_S12.1", "DNA03976_S13.1",
  "DNA03977_S14.1", "DNA03978_S15.1", "DNA03982_S33.1", "DNA03983_S34.1",
  "DNA03985_S35.1", "DNA03986_S36.1", "DNA03989_S37.1", "DNA03990_S38.1",
  "DNA03991_S39.1", "DNA03992_S40.1", "DNA03994_S41.1", "DNA03999_S42.1",
  "DNA04000_S43.1"
)

filtered_df <- metadata %>%
  filter(sample %in% sample_list) %>%
  mutate(
    subject = str_remove(subject, ":$"),
    LFS_status = case_when(
      subject %in% c("Patient", "Family member A", "Family member C") ~ "LFS",
      TRUE ~ "non-LFS"
    ),
    CTx_status = case_when(
      subject %in% c("Patient", "UW volunteer 7") ~ "CTx",
      TRUE ~ "no-CTx"
    )
  )

filtered_df %>% print(n=Inf)

write.table(
  filtered_df,
  file = "results/li_fraumeni_metadata.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


