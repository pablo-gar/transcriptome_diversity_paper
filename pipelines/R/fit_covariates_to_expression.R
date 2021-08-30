source('../../R/misc.R')
source('../../R/transcriptome_diversity_tools.R')
library('tidyr')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_file <- cmdArgs[1] # A gtex-like covariate file with only peer factors
    cov_file_or_dir <- cmdArgs[2] # A gtex-like covariate file
    out_dir <- cmdArgs[3]
    
    #exp_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_matrices_processed/tpm/gtex_Whole_Blood.txt'
    #cov_file_or_dir <- '/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt'
    
    #exp_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_matrices_processed/tpm/lin.txt'
    #cov_file_or_dir <- '/scratch/users/paedugar/transcriptome_diversity/PEER_covariates/lin.txt'
    
    
    # Get covariate files
    
    if(dir.exists(cov_file_or_dir)) {
        cov_files <- list.files(cov_file_or_dir, full.names=T)
    } else if (file.exists(cov_file_or_dir)) {
        cov_files <- cov_file_or_dir
    } else {
        stop(paste0(cov_file_or_dir, " does not exist"))
    }
    
    # Read expression matrix
    exp_mat <- read_expression(exp_file)
    exp_mat <- filter_rarely_expressed_genes(exp_mat)
    exp_mat <- rankit_normalize_expression(exp_mat)
    
    # correct gtex names if needed
    if(any(grepl('GTEX', colnames(exp_mat)))) {
        colnames(exp_mat) <- gsub('(GTEX\\-.+?)\\-.+', '\\1', colnames(exp_mat))
    }
    
    for (current_file in cov_files) {
        cat('Working with ', current_file, '\n')
        covs <- read.table(current_file, sep='\t', stringsAsFactors=F, header=T, check.names=F)
        covs <- rankit_normalize_expression(covs)
        cov_names <- covs[,1]
        
        # Put data together
        current_exp_mat <- exp_mat
        current_exp_mat[,1] <- gsub('-', '_', current_exp_mat[,1])
        colnames(current_exp_mat)[1] <- colnames(covs)[1]
        current_exp_mat <- current_exp_mat[,colnames(current_exp_mat) %in% colnames(covs)]
        genes <- as.character(current_exp_mat[,1])
        
        covs <- covs[, colnames(current_exp_mat)]
        
        merged <- rbind(covs,current_exp_mat)
        rownames(merged) <- merged[,1]
        merged <- merged[,-1]
        merged <- as.data.frame(t(merged))
        
        results <- rep(0, length(genes))
        for (i in 1:length(genes)) {
        #for (i in 1:200) {
            results[i] <- summary(lm(data=merged, formula=as.formula(paste(genes[i], '~', paste(cov_names, collapse=' + ')))))[['r.squared']]
        }
        
        results <- data.frame(gene=genes, r_squared=results, original_file=basename(current_file))
        write.table(results, file.path(out_dir,  basename(current_file)), sep='\t', quote=F, row.names=F)
    }
}

main()
