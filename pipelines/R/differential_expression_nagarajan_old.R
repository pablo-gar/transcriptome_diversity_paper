# Performs differential expression from nagarajan et al. accounting for transcriptome diversity
# Wed 15 Jul 2020 04:30:50 PM PDT

library('edgeR')
library('dplyr')
library('readr')
source('../../R/transcriptome_diversity_tools.R')


main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    cov_file <- cmdArgs[2]
    transcriptome_diversity_file <- cmdArgs[3]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan.txt'
    #cov_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan_cov.txt'
    #transcriptome_diversity_file_tpm <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/nagarajan_tpm.txt' 
    #transcriptome_diversity_file_tmm <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/nagarajan_tmm.txt' 
    #outdir <- '.'
    
    
    ###################################
    # Options (loads options table for this study)
    dge_options <- build_comparison_table(outdir)
    
    
    ###################################
    # Read input files 
    
    exp_mat<- read_tsv(exp_mat_file)
    
    transcriptome_diversity_tpm <- read_tsv(transcriptome_diversity_file_tpm)
    transcriptome_diversity_tmm <- read_tsv(transcriptome_diversity_file_tmm)
    
    covariates_orginal <- read_covariates_t(cov_file)
    
    # add trasncriptome diversity to covariates
    covariates_orginal <- left_join(covariates_orginal, transcriptome_diversity_tpm) %>%
        rename(transcriptome_diversity_tpm='transcriptome_diversity') %>%
        left_join(transcriptome_diversity_tmm) %>%
        rename(transcriptome_diversity_tmm='transcriptome_diversity') %>%
        rename(Sample='sample_id')
    
    for(i in 1:nrow(dge_options)) {
        
        current_options <- dge_options[i,,drop=F]
        
        # setting options
        
        differential_variable <- current_options$differential_variable
        ignore_indeces <- current_options$ignore_indeces
        outfile <- current_options$path
        covariate_variables <- current_options$covariate_variables[[1]]
        
        ###################################
        # Filtering samples
        covariates <- filter_covariates(covariates_orginal, current_options)
        exp_mat_final <- exp_mat[,c('gene_id', covariates$Sample)]
        
        ###################################
        # EdgeR
        
        # Building edgeR objects
        dge <- DGEList(exp_mat_final[,-1], genes=exp_mat_final[,1, keep=T], samples=covariates)
        
        # Filtering genes
        keep <- filterByExpr(dge, group=differential_variable)
        dge <- dge[keep, , keep.lib.size=F]
        
        ###
        # Exploratory analysis
        #plotMDS(dge)
        
        ###
        # Calculate differential gene expression w and w/o transcriptome diversity as covariates
        
        # Results is a list of lists with the following structure:
        # resluts[[transcriptome diversity included in model?]]
        #           [[design]] contains design matrix
        #           [[dge]] contains edge object after estimating disperssion
        #           [[fit]] contains edge the fit glmQLFit
        
        edgeR_results <- list()
        for(j in c('regular', 'transcriptome_diversity_tpm', 'transcriptome_diversity_tmm', 'transcriptome_diversity_tpm_scramble')) {
            
            if(j == 'regular') {
                to_use <- covariate_variables
            } else {
                to_use <- c(covariate_variables, j)
            }
            design_formula <- as.formula(paste0('~ ', paste(to_use, collapse=' + '), ' + ', differential_variable))
            edgeR_results[[j]] <- fit_dge(dge, design_formula, ignore_indeces=ignore_indeces)
            
        }
        
        saveRDS(edgeR_results, outfile)
    }
    
}

## Builds readme-like table containing all the pairwise comparisons 
build_comparison_table <- function(outdir) {
    
    x <- tibble(ID=c('genotype', 'tamoxifen_wt', 'jq1_wt', 'fulvestrant_wt'), 
                differential_variable=c('genotype', 'treatment', 'treatment', 'treatment'),
                treatment=c(NA, 'tamoxifen', 'jq1', 'fulvestrant'),
                covariate_variables=list(c('clone'), c('clone'), c('clone'), c('clone')),
                ignore_indeces=c(5,99,99,99))
    
    x$path <- file.path(outdir, paste0('edgeR_results_', x$ID, '.RDS'))
    
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
    
# Fits glm given an edgeR object and a design formula
# WARNING ignore_indeces is used to eliminate extra factors from desing matrix, it's hard coded and needs to be checked every time
fit_dge <- function(dge, design_formula, ignore_indeces) {
    
    design <- model.matrix(design_formula, dge$samples)
    message(paste('Eliminating the following from design matrix, make sure this is what you want: \n', paste(colnames(design)[ignore_indeces], collapse='\n')))
    #browser()
    design <- design[,-ignore_indeces, drop=F]
    
    ###
    # estimate dispersion
    dge <- estimateDisp(dge, design, robust=T)
    
    ###
    # Perform fit and find differential expressed genes
    fit <- glmQLFit(dge, design, robust=TRUE)
    differential_genes <- glmQLFTest(fit)
    differential_genes$table$FDR <- p.adjust(differential_genes$table$PValue, method='BH')
    
    return(list(dge=dge, design=design, fit=fit, differential_genes=differential_genes))
}

main()
