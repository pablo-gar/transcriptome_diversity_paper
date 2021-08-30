# Calculates PEER covariates from a TMM normalized matrix

source('../../R/transcriptome_diversity_tools.R')
library('peer')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    #exp_mat_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_matrices_processed/tmm/lin.txt'
    
    # Reading and filtering expression matrix
    exp_mat <- read_expression(exp_mat_file)
    exp_mat <- filter_rarely_expressed_genes(exp_mat)
    
    rownames(exp_mat) <- exp_mat[,1, drop=T]
    exp_mat <- t(as.matrix(exp_mat[,-1]))
    
    # Calculating the number of hidden factors based on GTEx guidelines
    if(nrow(exp_mat) < 150) {
        hiddenFactors <- 15
    } else if (nrow(exp_mat) < 250) {
        hiddenFactors <- 30
    } else if (nrow(exp_mat) < 350) {
        hiddenFactors <- 45
    } else {
        hiddenFactors <- 60
    }
    
    # Getting PEER covariates
    model = PEER()
    PEER_setPhenoMean(model, exp_mat)
    PEER_setNk(model, hiddenFactors)
    PEER_setNmax_iterations(model, 500) 
    PEER_update(model)

    columnsHidden <- paste0("InferredCov", 1:hiddenFactors)

    factors <- as.data.frame(PEER_getX(model))
    rownames(factors) <- rownames(exp_mat)
    colnames(factors) <- columnsHidden
    factors <- t(factors)
    factors <- cbind(data.frame(ID=rownames(factors)), as.data.frame(factors))

    
    # Write results
    write.table(factors, out_file, sep='\t', row.names=F, col.names=T, quote=F)


    
}

main()
