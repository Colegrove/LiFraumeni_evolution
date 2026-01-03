COSMIC_MAF <- 
  loadMaf(inputs$COSMIC_MAF,
          inputs$COSMIC_Seshat,
          testBed,
          seshat_col_trans,
          seshat_path_trans,
          variant_clasification_table,
          alphamissense,
          intOGen
  ) %>%
  bind_cols(
    read_delim(inputs$COSMIC_VCF, "\t", 
               escape_double = FALSE, 
               comment = "##", 
               trim_ws = TRUE) %>% 
      # Rename #CHROM
      dplyr::select(CHROM = `#CHROM`, everything()) %>% 
      dplyr::rename(FILTER_VCF = FILTER)
  )

COSMIC_all <-   
  # list.files(inputs$COSMIC_dir,pattern=inputs$COSMIC_pattern, full.names = T, recursive = T) %>% 
  # Load all COSMIC files and concatenate
  read_csv(inputs$COSMIC_file, col_types = cols(.default = "c" )) %>%
  type_convert() %>% 
  filter(PRIMARY_HISTOLOGY=="carcinoma") %>% 
  filter(GENOME_WIDE_SCREEN=="y") %>%
  filter(SAMPLE_TYPE != "cell line") %>% 
  mutate(rownum = row_number()) %>% 
  # filter(rownum == 5285) %>% 
  left_join(COSMIC_MAF %>% dplyr::select(-FILTER, -INFO) %>% unique(), by=c("GENOMIC_MUTATION_ID" = "ID"
  )) %>% 
  filter(!is.na(Mutation_type))

COSMIC_file <- COSMIC_all %>% 
  filter(grepl("Substitution", MUTATION_DESCRIPTION))

