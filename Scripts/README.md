# PDMR Landscape — Figure Scripts

This folder contains R scripts used to generate the main and supplementary figures for the PDMR landscape paper.


## Inputs

- Large data files in Figshare https://figshare.com/s/a47ea87939355f7cc558 
- Other inputs in the repository’s `DATA/` folder.

## For the scripts in sup_code folder, run them in sequence

1. compute_pca_model_level_expression.R
2. compute_PDX_TCGA_CCLE_data.R
3. do_pca_for_PDX_TCGA_CCLE.R
4. prepare_PDX_decoder_table.R
5. run_decoder_PDX.sh
    * Note: The decoder package is not included due to size. The shell script contains the command to download it.

## The generated figures are in Figures folder

