
#### Hunter Colegrove
## 24 NOV 2025
#### Takes mutation dataframe and saves as .maf file for use in SigProfiler algorithms


dir_path_all <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mut_sig_2025_11_24/all_muts"
dir_path_coding <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mut_sig_2025_11_24/coding"
dir_path_non_coding <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mut_sig_2025_11_24/non_coding"
outfile_all <- "maf_matrix_all_24NOV2025.maf"
outfile_coding <- "maf_matrix_coding_24NOV2025.maf"
outfile_non_coding <- "maf_matrix_non_coding_24NOV2025.maf"

# all muts
maf_format <- maf_masked_coding %>%
  filter(Tissue != "Urine cells") %>%
  add_column(placeholder1 = ".", .before = 2) %>%
  add_column(placeholder2 = ".", .before = 3) %>%
  add_column(placeholder3 = ".", .before = 8) %>%
  add_column(placeholder4 = ".", .before = 12) %>%
  add_column(placeholder5 = ".", .before = 15)

new_file_path <- file.path(dir_path_all,outfile_all)
cat(paste("#version 2.4\n"), file = new_file_path)
write_delim(maf_format, new_file_path, delim = "\t", append = TRUE, col_names = TRUE)
  

## coding only
maf_format_coding <- maf_format %>% 
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site")

new_file_path_coding <- file.path(dir_path_coding,outfile_coding)
cat(paste("#version 2.4\n"), file = new_file_path_coding)
write_delim(maf_format_coding, new_file_path_coding, delim = "\t", append = TRUE, col_names = TRUE)


## non-coding only
maf_format_non_coding <- maf_format %>%
  filter(is.na(exon_number) & Variant_Classification != "Splice_Site")

new_file_path_non_coding <- file.path(dir_path_non_coding,outfile_non_coding)
cat(paste("#version 2.4\n"), file = new_file_path_non_coding)
write_delim(maf_format_non_coding, new_file_path_non_coding, delim = "\t", append = TRUE, col_names = TRUE)

