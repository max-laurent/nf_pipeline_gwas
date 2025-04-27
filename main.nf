params.bfile = "./data/geno_data*"
params.pheno_matrix = "./data/expression_matrix.txt"
params.kin = "./data/FASTLMM_kinship_sansChr11.txt"
params.bim = "./data/Matrix_50K_600K_GBS_BGA_SVAmaizing_Dente_NoPrivate_Chr.bim.all"
params.info = "./data/Matrix_ImputedBeagle012_50K_600K_GBS_BGA_SV_Amaizing_Dente_Chr_All_OneMatrix_NotFiltered_NoPrivate.rds"
params.position = "./data/2023-11-17_AmaizingV3_Info_SNPs_InDels_WithConfInt_R2_r2k_SeuilR2_0.1_Model_HillWeir_cor_interval_RefGen_v4.txt" 
params.outdir = "./data/results" 

process SPLIT_PHENO {

    input:
    path matrix

    output:
    path 'Pheno*'

    script: // Improvement could be done by using the method .baseName on the output 
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
    //container "image_fastlmm_plink:latest"
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
    Channel
        .fromPath(params.pheno_matrix)
        .set{expr_matrix_ch}
    Channel
    	.fromPath(params.bfile)
        .toSortedList()
        .flatten()
    	.map { it -> [it.name.split('\\.').first(), it] }
    	.groupTuple()
        .collect()
    	.set{geno_files_ch}
    Channel
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
        
    gwas_results_ch = FASTLMM(split_pheno_ch, geno_files_ch, kin_ch)
    removal_ch = TO_REMOVE_FOR_MAF(split_pheno_ch)
    maf_ch = CALCULATE_MAF(removal_ch, geno_files_ch) 

    gwas_results_ch
        .concat(maf_ch)
        .groupTuple()
        .set{gwas_results_ch}

    bim = Channel.fromPath(params.bim).collect()
    info = Channel.fromPath(params.info).collect()
    position = Channel.fromPath(params.position).collect()

    formated_data = ADD_INFO_GWAS(gwas_results_ch, bim, info, position)
    signif_data = FILTER_SIGNIF_MARKERS(formated_data)
}