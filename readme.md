# Li-Fraumeni evolution

This repository contains the analysis pipeline used to investigate somatic mutation and selection in normal and cancer tissues from individuals with Li-Fraumeni syndrome as part of the pre-print: https://doi.org/10.64898/2026.01.12.699071. 
The main branch of this repository contains the code used to generate all analysis within the publication using twinstrand's nexus-hosted variant calling pipeline. 
Due to potential unavailability to access the twinstrand analysis pipeline for general use, the <deepUMI_branch> of this repository contains code used to produce analysis using input data generated from the <deepUMI_tool> pipeline.

## Pipeline overview
This repository is organized as an RStudio project and is intended to be run in R.
Analyses are run using a main script (scripts/<main_script>), which sources the analysis scripts in a fixed order. 
Note: The order of script execution within the main script must be preserved, as later steps depend on objects, annotations, and intermediate results generated throughout the workflow. 

1. Open the Rproject file LiFraumeni_evolution.Rproj
2. Install required packages (see supplementary table for list of package versions: <link_to_supplement>)
3. Supply required reference files: see following section for list of necessary files and recommended locations.
4. Configure inputs: edit processing_config.txt file to specify input data paths (see sections below for input data types) and additional data filtering parameters.
5. Run scripts/<main_script>

## Required reference files
Large reference files are excluded from version control and must be provided locally prior to running pipeline:

1. Human reference genome, hg38.fa: to be placed in inputs/refs/
- <link_to_data>
2. All possible mutations table, all_possible_sites_annotated.tsv.gz: to be placed in inputs/
- A file containing all possible single-nucleotide variants must be generated and annotated with VEP
- <link_to_VEP_webtool>
3. Gene annotation file, gencode.v38.annotation.gff3.gz: to be placed in inputs/
- <link_to_download_file>
4. Alphamissense annotations file, <file_name>: to be placed in inputs/alphamissense/
- <link_to_download_file>

## Input clinical data
Input sequencing datasets required to reproduce analysis are available at: <link_to_data>
Sequence alignment and variant calling will need to be performed on the hosted .fastq files prior to the use of this pipeline. 
We recommend the deepUMI analysis pipeline: <link_to_pipeline>

## Input synthetic data
To test the pipeline operation without access to patient data, a .maf file containing randomly simulated single nucleotide variants can be generated. 
These mutations are not biologically observed and are intended for validating pipeline operation. 

1. Run: <synthetic_data_script> to generate .maf file
2. Assign synthetic data to config file. 
 
Note: 
- All reference files will still be required prior to pipeline operation.
- Some analyses such as dinucleotide variant and phasing analysis that require BAM files will not be produced.
