#!/bin/bash

# Phenotypic data needs to be split by phenotype. This results in one file per phenotype.
filename="$1"

# Get the header row (column names)
n_col=4 #$(awk 'END {print NF;}' $filename)
col_names=$(head -n1 $filename)
# Loop through phenotype columns (starting from the first column containing phenotypic information, the first two contain genotype info)
for (( i=3; i<=$n_col; i++ )); do
    # Extract the column name based on its position
    out_name=$(echo $col_names | awk -v col="$i" '{print $col;}')
    # Extract the first two columns (IDs) and the phenotype column, and write to a new file
    awk -v col="$i" '{print $1, $2, $col}' $filename > Pheno_$out_name
done