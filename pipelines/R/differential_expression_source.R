read_covariates_all <- function(cov_file, transcriptome_diversity_file_tpm, transcriptome_diversity_file_tmm=NULL) {
    
    transcriptome_diversity_tpm <- read_tsv(transcriptome_diversity_file_tpm)
    
    covariates_orginal <- read_covariates_t(cov_file)
    
    # add trasncriptome diversity to covariates
    covariates_orginal <- covariates_orginal %>%
        left_join(transcriptome_diversity_tpm) %>%
        rename(transcriptome_diversity_tpm='transcriptome_diversity')
    
    # Read tmm based covariates if available
    if(!is.null(transcriptome_diversity_file_tmm)) {
        transcriptome_diversity_tmm <- read_tsv(transcriptome_diversity_file_tmm)
        covariates_orginal <- covariates_orginal %>%
            left_join(transcriptome_diversity_tmm) %>%
            rename(transcriptome_diversity_tmm='transcriptome_diversity')
    }

    covariates_orginal <- rename(covariates_orginal, Sample='sample_id')
    
    return(covariates_orginal)
    
}

dge_engine <- function(exp_mat, covariates_orginal, dge_options, covariate_types) {
    
    ## REQUIRE A custom build filter_covariates() function
    
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
        for(j in covariate_types) {
            
            if(j == 'regular') {
                to_use <- c(covariate_variables, differential_variable)
            } else {
                to_use <- c(covariate_variables, j, differential_variable)
            }
            
            to_use <- to_use[to_use!='']
            
            design_formula <- as.formula(paste0('~ ', paste(to_use, collapse=' + ')))
                
            edgeR_results[[j]] <- fit_dge(dge, design_formula, ignore_indeces=ignore_indeces)
            
        }
        
        saveRDS(edgeR_results, outfile)
    }
    
    
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
