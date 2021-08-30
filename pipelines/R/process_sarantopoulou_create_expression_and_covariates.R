library('readr')
library('tidyr') 
library('dplyr') 
library('stringr') 

main <- function(cmdArgs=commandArgs(T)) {
    
    original_file <- cmdArgs[1]
    output_exp_file <- cmdArgs[2]
    output_cov_file <- cmdArgs[3]
    
    #original_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/sarantopoulou/GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Truseq.txt.gz'
    #output_exp_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/sarantopoulou_truseq.txt'
    #output_cov_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/sarantopoulou_truseq_cov.txt'
    
    exp_mat <- read_tsv(original_file)
    
    # Building covariate matrix
    covariates <- colnames(exp_mat)
    # Getting covariate names based on sample ids
    if(any(grep('Illumina', covariates))) {
        covariates <- str_split_fixed(str_remove(covariates, 'Illumina.'),'\\_', 3)
        covariates[,3] <- 'Illumina'
    } else {
        covariates <- str_split_fixed(covariates,'\\.', 3)
    }
    
    # Selerecting only covariates (e.g. igonre gene ids)
    covariates <- covariates[rowSums(covariates=='') == 0, ]
    colnames(covariates) <- c('strain', 'individual', 'library')
    
    # build sample id
    if(any(grep('Illumina', covariates))) {
        covariates_df <- data.frame(sample_id=apply(covariates,1, function(x) paste0(x[3], '.', x[1], '_', x[2])))
    } else {
        covariates_df <- data.frame(sample_id=apply(covariates,1, paste, collapse='.'))
    }
    
    # Getting sequencing depth
    covariates_df$mapped_reads <- colSums(exp_mat[,covariates_df$sample_id])
    
    covariates <- cbind(covariates_df, covariates)
    covariates <- t(covariates)
    colnames(covariates) <- covariates[1,]
    covariates <- as.data.frame(cbind(ID=rownames(covariates)[-1], covariates[-1,]))
    
    
    # Processing matrix
    exp_mat <- select(exp_mat, -c(geneCoordinate, geneSymbol)) %>%
        rename(gene_id='id')
    
    write_tsv(exp_mat, output_exp_file)
    write_tsv(covariates, output_cov_file)
    
    
    
}

main()
