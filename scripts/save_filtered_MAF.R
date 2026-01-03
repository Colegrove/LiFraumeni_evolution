
#### Hunter Colegrove
## 25 APR 2025
#### Script filters and recreates .maf files as they are in R pipeline for use in SigProfiler algorithms

filt_maf %>% filter(Hugo_Symbol == 'TP53') %>% filter(Tissue != "Urine cells") %>% arrange(desc(NF)) %>% print(width = Inf)

filt_maf %>% filter(Tissue == "Skeletal muscle") %>% print(width = Inf)

patient_tissues_with_mutagenesis <- c("Whole blood", "Buffy coat", "Plasma", 
                                      "Bone marrow", "Colon", "Skin", 
                                      "Skin, non-sun-exposed", "Liver", 
                                      "Esophagus 1", "Esophageal cancer 1", 
                                      "Esophageal cancer 2", "Liver metastasis 1", 
                                      "Liver metastasis 2", "Lung metastasis", 
                                      "Mediastinal metastasis")
controls_patient_blood_samples <- c("PBMC", "Buffy coat")

#### Filtering criteria
VAF_FILT = 0.3

save_filtered_maf <- function(panel, tissues, out_path){
  
  if(panel == "ALL"){
    if(tissues == "blood"){
      barcodes_list = filt_maf %>%
        filter(Tissue %in% controls_patient_blood_samples) %>%
        distinct(Tumor_Sample_Barcode) %>%
        pull(Tumor_Sample_Barcode)
    }
    
    if(tissues == "tissues"){
      barcodes_list = filt_maf %>%
        filter(Tissue %in% patient_tissues_with_mutagenesis) %>%
        distinct(Tumor_Sample_Barcode) %>%
        pull(Tumor_Sample_Barcode)
    }
  }
  # }else if(panel == "MUT"){
  #   barcodes_list = filt_maf %>%
  #     filter(Tissue %in% controls_patient_blood_samples) %>%
  #     distinct(Tumor_Sample_Barcode) %>%
  #     pull(Tumor_Sample_Barcode)
  # }
  
  for(barcode in barcodes_list){
    sample_df <- filt_maf %>%
      filter(Tumor_Sample_Barcode == barcode) %>%
      filter(VAF < VAF_FILT) %>%
      add_column(placeholder1 = ".", .before = 2) %>%
      add_column(placeholder2 = ".", .before = 3) %>%
      add_column(placeholder3 = ".", .before = 8) %>%
      add_column(placeholder4 = ".", .before = 12) %>%
      add_column(placeholder5 = ".", .before = 15)
    new_file_name <- paste(barcode,".filtered.maf", sep="")
    new_file_path <- file.path(out_path,new_file_name)
    cat(paste("#version 2.4\n"), file = new_file_path)
    
    print(sample_df, width = Inf)
    write_delim(sample_df, new_file_path, delim = "\t", append = TRUE, col_names = TRUE)
  }
}

############## All available panels
panel = "ALL"

## blood only
tissues = "blood"
out_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered"
save_filtered_maf(panel, tissues, out_path)


## tissues only
tissues = "tissues"
out_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered"
save_filtered_maf(panel, tissues, out_path)





############## Mutagenesis panel only
panel = "MUT"

## tissues only
tissues = "tissues"
out_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_mutagenesis_filtered"
save_filtered_maf(panel, tissues, out_path)









file1 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/output/SBS/li_fraumeni.SBS96.all"
file2 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/output/SBS/li_fraumeni.SBS96.all"

mat1 <- read.table(file1, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
mat2 <- read.table(file2, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

if (!identical(rownames(mat1), rownames(mat2))) {
  stop("Row names (mutation contexts) do not match between the two matrices.")
}

combined <- cbind(mat1, mat2)
nrow(combined)
duplicated_columns <- colnames(combined_fixed)[duplicated(colnames(combined_fixed))]

dup_indices <- which(colnames(combined) == "DNA03975_S12.1")
combined_fixed <- combined[, -dup_indices[2]] 

combined_fixed
