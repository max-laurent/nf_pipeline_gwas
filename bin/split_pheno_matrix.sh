#!/bin/bash

filename="$1"
n_col=4 #$(awk 'END {print NF;}' $filename)
col_names=$(head -n1 $filename)

for (( i=3; i<=$n_col; i++ )); do
    out_name=$(echo $col_names | awk -v col="$i" '{print $col;}')
    awk -v col="$i" '{print $1, $2, $col}' $filename > Pheno_$out_name
done