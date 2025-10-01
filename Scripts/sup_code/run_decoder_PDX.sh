#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --mem=128g

git clone https://github.com/laurapeng/decoder.git #clone the decoder software

mv HNSC_originator_expression.tsv decoder/data/
mv PAAD_originator_expression.tsv decoder/data/
mv COAD_originator_expression.tsv decoder/data/

mv config_de_novo_HNSC.tsv decoder
mv config_de_novo_PAAD.tsv decoder
mv config_de_novo_COAD.tsv decoder

module purge 
module load matlab
configfile="./decoder/config_de_novo_HNSC.tsv" 
matlab -batch "Decon_de_novo('$configfile')"

configfile="./decoder/config_de_novo_PAAD.tsv" 
matlab -batch "Decon_de_novo('$configfile')"

configfile="./decoder/config_de_novo_COAD.tsv" 
matlab -batch "Decon_de_novo('$configfile')"