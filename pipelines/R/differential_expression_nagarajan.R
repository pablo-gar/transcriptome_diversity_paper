# Performs differential expression from nagarajan et al. accounting for transcriptome diversity
# Wed 15 Jul 2020 04:30:50 PM PDT

library('edgeR')
library('dplyr')
library('readr')
source('../../R/transcriptome_diversity_tools.R')
source('differential_expression_source.R')


main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    cov_file <- cmdArgs[2]
    transcriptome_diversity_file_tpm <- cmdArgs[3]
    outdir <- cmdArgs[4]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan.txt'
    #cov_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan_cov.txt'
    #transcriptome_diversity_file_tpm <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/nagarajan_tpm.txt' 
    #transcriptome_diversity_file_tmm <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/nagarajan_tmm.txt' 
    #outdir <- '.'
    
    ###################################
    # Options (loads options table for this study)
    dge_options <- build_comparison_table(outdir)
    covariate_types <- c('regular', 'transcriptome_diversity_tpm', 'transcriptome_diversity_tpm_scramble')
    #covariate_types <- c('regular', 'transcriptome_diversity_tpm', 'transcriptome_diversity_tpm_scramble')
    
    ###################################
    # Read input files 
    
    exp_mat<- read_tsv(exp_mat_file)
    covariates_orginal <- read_covariates_all(cov_file, transcriptome_diversity_file_tpm)
    
    dge_engine(exp_mat, covariates_orginal, dge_options, covariate_types)
    
    
}

## Builds readme-like table containing all the pairwise comparisons 
build_comparison_table <- function(outdir) {
    
    x <- tibble(ID=c('genotype', 'tamoxifen_wt', 'jq1_wt', 'fulvestrant_wt'), 
                differential_variable=c('genotype', 'treatment', 'treatment', 'treatment'),
                treatment=c(NA, 'tamoxifen', 'jq1', 'fulvestrant'),
                covariate_variables=list(c('clone'), c('clone'), c('clone'), c('clone')),
                ignore_indeces=c(5,99,99,99))
    
    x$path <- file.path(outdir, paste0('edgeR_results_', x$ID, '.Rds'))
    
    return(x)
}


filter_covariates <- function(covariates, opts) {
    
    # Keep only final lineage
    covariates <- covariates %>%
        dplyr::filter(lineage=='final')
    
    if(opts$differential_variable=='genotype') {
        
        # Keep only no treatment
        covariates <- covariates %>%
            dplyr::filter(treatment=='vehicle') %>%
            group_by(genotype) %>%
            mutate(transcriptome_diversity_tpm_scramble=sample(transcriptome_diversity_tpm)) %>%
            ungroup()
        
    } else if(opts$differential_variable=='treatment') {
        
        # Keep only wt
        covariates <- covariates %>%
            dplyr::filter(genotype=='wt') %>%
            dplyr::filter(treatment %in% c(opts$treatment, 'vehicle'))  %>%
            mutate(transcriptome_diversity_tpm_scramble=sample(transcriptome_diversity_tpm)) %>%
            ungroup()
        
    }
    
    
    return(covariates)
}
    

main()
