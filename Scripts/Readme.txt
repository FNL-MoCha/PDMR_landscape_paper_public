Instructions:

0. Some large data required to run some scripts are located in Figshare (https://figshare.com/account/items/29853116/). The other required data are in DATA folder.

1. For sup_code, run these in sequence:
compute_pca_model_level_expression.R
compute_PDX_TCGA_CCLE_data.R
do_pca_for_PDX_TCGA_CCLE.R
prepare_PDX_decoder_table.R
run_decoder_PDX.sh # run this on HPC

*The decoder package is too large to include, but the provided code includes the command to download it.


2. Run each of the scripts named by figure number to get the figure pdfs.

3. The generated figures are in Figures folder.
