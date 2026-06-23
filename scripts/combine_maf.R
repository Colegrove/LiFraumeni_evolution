### Hunter Colegrove
### 18 Sep 2024
### Li-Fraumeni
###

### combine all maf files with removed INFO contents and SV mutations 
### output a maf file containing all samples: allSamples_noINFO.maf


## path of files
dat_path <- "/inputs/mafs/"

if (! file.exists("inputs/mafs/allSamples_noINFO.maf") ) {
  ## find all relevant mafs
  maf_NoINFO_list <- list.files(dat_path, pattern = "\\.variant-calls\\.noINFO\\.maf$", recursive = TRUE, full.names = TRUE)
  
  ## combine the mafs - PHENO column has interfile discrepancy call string to not remove values
  combined_data <- maf_NoINFO_list %>% 
    map_dfr(~ read_delim(.x, delim = "\t", skip = 1, col_types = cols(PHENO = col_character())))
  
  print(combined_data)
  output_path <- "/inputs/mafs/"
  out_filename <- "allSamples_noINFO.maf"
  out_filepath <- glue("{output_path}{out_filename}")
  print(out_filepath)
  write_delim(combined_data, out_filepath, delim='\t')
}