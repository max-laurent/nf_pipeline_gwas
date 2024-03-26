process SPLIT_PHENO {

    input:
    path matrix

    output:
    path 'Pheno*'

    shell:
    '''
    #!/bin/bash

    filename="!{matrix}"
    n_col=$(awk 'END {print NF;}' $filename)
    col_names=$(head -n1 $filename)

    for (( i=2; i<=$n_col; i++ )); do
        out_name=$(echo $col_names | awk -v col="$i" '{print $col;}')
        awk -v col="$i" '{print $1, $col}' $filename > Pheno_$out_name.txt
    done
    '''
}

process FORMAT_PHENO {
    debug true
    input:
    path geno_file
    path pheno_file
    

    output:
    path 'FORMATED_*'

    script:
    // """
    // ls $pheno_file
    // """
    """
    #!/usr/bin/env R --vanilla
    library("tidyverse")
    pheno_file <- read.table(file = "$pheno_file", header = TRUE, sep = " ")

    pheno_file <- pheno_file %>% rename(FID = genotype)

    geno_file <- read.table("$geno_file", stringsAsFactors = FALSE, header = FALSE)

    frame <- geno_file[, c(1:2)]
    frame <- left_join(frame, pheno_file, by = join_by("V2" == "FID"))
    frame[is.na(frame)] <- -9

    colnames(frame)[1] <- "FID"
    colnames(frame)[2] <- "IID"

    write.table(frame, file = paste0("FORMATED_", "$pheno_file"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    """

}
workflow {
    Channel
        .fromPath('./data/expression_matrix.txt')
        .set{expr_matrix_ch}
    Channel
        .fromPath('./data/geno_data.fam')
        .set{fam_data}
    
    split_pheno_ch = SPLIT_PHENO(expr_matrix_ch)
    formated_pheno_ch = FORMAT_PHENO(fam_data.collect(), split_pheno_ch.flatten()) 
}