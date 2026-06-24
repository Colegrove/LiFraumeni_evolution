# Li-Fraumeni evolution
This repository contains the analysis pipeline used to generate analysis and figures as part of the pre-print:  
[Ultra-deep duplex sequencing reveals unique features of somatic evolution in the normal tissues of a family with Li-Fraumeni syndrome](https://doi.org/10.64898/2026.01.12.699071). 

Data to run pipeline and recreate figures will be available at dbGaP accession phs004484.v1.p1 upon publication.

Sample .BAM, .MAF, and .MUT files will be required to run analysis in entirety.

## Pipeline overview
This repository is organized as an RStudio project and is intended to be run in R with the exception of phasing analysis run as a python script.

Analyses are run using a main script (scripts/00.0_MainScript.R).  
Note: The order of script execution within the main script must be preserved, as later steps depend on objects, annotations, and intermediate results generated throughout the workflow. 

### Using the pipeline
1. Clone this repository: `git clone https://github.com/Colegrove/LiFraumeni_evolution.git`
2. Open `LiFraumeni_evolution.Rproj` in RStudio
3. Install required packages:
- Run `setup.R` (at the repo root) once in R to install all dependencies automatically.
- For exact package versions used in the published analysis, see [pre-print supplemental materials](https://www.biorxiv.org/content/10.64898/2026.01.12.699071v1.supplementary-material).
4. Supply required reference files:
- See following section for list of necessary files and recommended locations.
5. Configure inputs:
- Edit processing_config.txt file to specify input data paths and additional data filtering parameters.
6. Run scripts/00.0_MainScript.R
- Note: To generate the final figure, phasing_tp53_181_indels.py will need to be run in python between tissue_phasing_prepare.R and tissue_phasing.R. 

## Required input files
Large input files are excluded from version control and must be provided prior to running pipeline.
Below are suggested filepath locations for input files. Edit processing_config.txt if input files are located in different locations.

1. Human reference genome; inputs/refs/hg38.fa
2. All possible mutations table; inputs/all_possible_sites_annotated.tsv.gz:
- A file containing all possible single-nucleotide variants must be generated and annotated with VEP
- https://www.ensembl.org/vep
3. Gene annotation file, inputs/gencode.v38.annotation.gff3.gz:
- https://www.gencodegenes.org/human/release_48.html
4. Alphamissense annotations file, inputs/alphamissense/AlphaMissense_hg38.tsv.gz:
- https://zenodo.org/records/8208688
