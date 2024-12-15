#!/bin/bash

# Define input file paths and output directory
expression_matrix="../data/input/expression_matrix.csv"
metadata="../data/input/metadata.csv"
output_dir="../data/cache"

# Define output prefix and number of shuffles
output_prefix="The_compare"
n_shuffles=5

# Run the analysis script
python run_analysis.py \
--expression_matrix $expression_matrix \
--metadata $metadata \
--output $output_dir \
--output_prefix $output_prefix \
--n_shuffles $n_shuffles \
--metadata-filter sample GSM6886536_Sigmoid_03N GSM6886537_Sigmoid_03T \
--metadata-filter CellType "IgA+ Plasma cells"

# Run the R comparative analysis script
Rscript compare.R --comparison_list ../data/info.txt --output_dir ../data/report1




