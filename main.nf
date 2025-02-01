#!/usr/bin/env nextflow

process SPLIT_PHENO {

    input:
    path matrix

    output:
    path 'Pheno*'

    script: // TODO: Improvement could be done by using the method .baseName on the output 
    """
    split_pheno_matrix.sh $matrix 
    """
}

process FORMAT_PHENO {
    
    input:
    tuple val(sample_id), path(geno_file)
    path pheno_file
    

    output:
    path 'FORMATED_*'

    script:
    """
    matrix_formating.R $pheno_file ${geno_file[2]}
    """
    
}

process FASTLMM {
    container "image_fastlmm_plink:latest"
    maxForks 5
    publishDir params.outdir, mode: 'copy'
    memory = '8GB'

    input:
    tuple val(meta), path(pheno)
    tuple val(sample_id), path(geno_file)
    path kin 

    output:
    tuple val(meta), path("FASTLMM_${meta}.txt")

    script:
        """
        fastlmmc -REML -verboseOut -bfile $sample_id -pheno $pheno -sim $kin -simLearnType Full -out FASTLMM_${meta}.txt -maxThreads 8 -mpheno 1
        """
}

process TO_REMOVE_FOR_MAF {
    input:
    tuple val(meta), path(pheno)

    output:
    tuple val(meta), path("MISSING_DATA.txt")

    shell:
    '''
    #!/bin/bash
    # Save the lines with missing value for phenotype in another file to avoid mistakes when calculating the MAF
    awk '$3 == -9' !{pheno} > MISSING_DATA.txt
    '''

}

process CALCULATE_MAF {
    input:
    tuple val(meta), path("MISSING_DATA.txt")
    tuple val(sample_id), path(geno_file)
    

    output:
    tuple val(meta), path("*.frq")

    script:
    """
    # Calculate the MAF for the markers for a given phenotype while removing individuals with missing measurment
    plink \
        --noweb \
	    --bfile $sample_id \
	    --remove MISSING_DATA.txt \
	    --freq \
	    --out freq_$meta
    """

}

process ADD_INFO_GWAS {
    
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(pheno)
    path bim
    path info
    path position

    output:
    tuple val(meta), path("Formated_${meta}.txt")

    script:
    """
    formating_fastlmm_output.R $meta ${pheno[0]} ${pheno[1]} $bim $info $position 
    """
}

process FILTER_SIGNIF_MARKERS {
    
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(pheno)

    output:
    path "Signif_markers_${meta}.txt"

    script:
    """
    extract_signif_marker.R $meta ${pheno[0]} 
    """
}

workflow {
    Channel // Create channel for the phenotyping matrix
        .fromPath(params.pheno_matrix)
        .set{expr_matrix_ch}
     Channel // Create channel for the bitfiles
    	.fromPath(params.bfile) // Create a channel from the pattern of your file bitfile
        .toSortedList()         // Collect all matching files into a sorted list
        .flatten()              // Flatten the list in case of nested structures
    	.map { it -> [it.name.split('\\.').first(), it] } // Extract the filename prefix (before the first dot) and pair it with the file
    	.groupTuple()           // Group files by their common prefix (e.g., sample1 -> [sample1.bed, sample1.bim, sample1.fam])
        .collect()              // Collect all grouped files into a single list
    	.set{geno_files_ch}
    Channel // Create channel for the kinship matrix
        .fromPath(params.kin)
        .collect()
        .set{kin_ch} 

    formated_pheno_ch = FORMAT_PHENO(geno_files_ch, expr_matrix_ch) 
    split_pheno_ch = SPLIT_PHENO(formated_pheno_ch)

    split_pheno_ch
        .flatten()
        .map{it -> [it.name.split('\\/').last(), it]}
        .groupTuple()
        .set{split_pheno_ch}
    
    // Run the GWAS
    gwas_results_ch = FASTLMM(split_pheno_ch, geno_files_ch, kin_ch)
    // MAF 
    removal_ch = TO_REMOVE_FOR_MAF(split_pheno_ch)
    maf_ch = CALCULATE_MAF(removal_ch, geno_files_ch) 

    // Group the GWAS results with the MAF information per phenotype
    gwas_results_ch
        .concat(maf_ch)
        .groupTuple()
        .set{gwas_results_ch}

    bim = Channel.fromPath(params.bim).collect() // Information about all the markers and their alleles
    info = Channel.fromPath(params.info).collect() // R file containing the allele information for the varieties
    position = Channel.fromPath(params.position).collect() // Convert the position from v2 to v4

    formated_data = ADD_INFO_GWAS(gwas_results_ch, bim, info, position) // Convert marker positions to v4 and add the allele information
    signif_data = FILTER_SIGNIF_MARKERS(formated_data) // Export the significant markers
}