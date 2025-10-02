# PDMR Landscape — Figure Scripts:

This folder contains R scripts used to generate the main and supplementary figures for the PDMR landscape paper.

## Input files:

Some large data required to run above scripts are located in Figshare (https://figshare.com/account/items/29853116/). The other input files are in DATA folder.

## For the scripts in sup_code folder, run them in sequence
compute_pca_model_level_expression.R
compute_PDX_TCGA_CCLE_data.R
do_pca_for_PDX_TCGA_CCLE.R
prepare_PDX_decoder_table.R
*run_decoder_PDX.sh
    *The decoder package is too large to include, but the provided code includes the command to download it.

## The generated figures are in Figures folder

