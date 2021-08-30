# Performs differential expression from nagarajan et al. accounting for transcriptome diversity
# Wed 15 Jul 2020 04:30:50 PM PDT

library('edgeR')
library('dplyr')
library('readr')
source('../../R/transcriptome_diversity_tools.R')
source('./differential_expression_source.R')


main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    cov_file <- cmdArgs[2]
    transcriptome_diversity_file_tpm <- cmdArgs[3]
    outdir <- cmdArgs[4]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/sarantopoulou_truseq.txt'
    #cov_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/sarantopoulou_truseq_cov.txt'
    #transcriptome_diversity_file_tpm <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/sarantopoulou_truseq_tpm.txt' 
    #transcriptome_diversity_file_tmm <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/sarantopoulou_truseq_tmm.txt' 
    #outdir <- '.'
    
    
    ###################################
    # Options (loads options table for this study)
    dge_options <- build_comparison_table(outdir)
    #covariate_types <- c('regular', 'transcriptome_diversity_tpm', 'transcriptome_diversity_tmm', 'transcriptome_diversity_tpm_scramble')
    covariate_types <- c('regular', 'transcriptome_diversity_tpm', 'transcriptome_diversity_tpm_scramble')
    
    ###################################
    # Read input files 
    
    exp_mat<- read_tsv(exp_mat_file)
    covariates_orginal <- read_covariates_all(cov_file, transcriptome_diversity_file_tpm)
    
    dge_engine(exp_mat, covariates_orginal, dge_options, covariate_types)
    
    
}

## Builds readme-like table containing all the pairwise comparisons 
build_comparison_table <- function(outdir) {
    
    x <- tibble(ID=c('ILB'), 
                differential_variable=c('strain'),
                treatment=c(NA),
                covariate_variables=list(c('')),
                ignore_indeces=c(99)
                )
    
    x$path <- file.path(outdir, paste0('edgeR_results_', x$ID, '.Rds'))
    
    return(x)
    
}


filter_covariates <- function(covariates, opts) {
    
    if(opts$differential_variable=='strain') {
        
        keep_going <- T
        
        while(keep_going) {
            # Keep only no treatment
            covariates <- covariates %>%
                group_by(strain) %>%
                mutate(transcriptome_diversity_tpm_scramble=sample(transcriptome_diversity_tpm)) %>%
                ungroup()
            
            if(any(covariates$transcriptome_diversity_tpm == covariates$transcriptome_diversity_tpm_scramble)) {
                keep_going <- T
            } else {
                keep_going <- F
            }
        }
        
        
    }     
    
    return(covariates)
}
    

main()
