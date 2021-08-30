source('../../R/misc.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    factors_file <- cmdArgs[1] # A gtex-like covariate file
    peer_file <- cmdArgs[2] # A gtex-like covariate file with only peer factors
    out_prefix <- cmdArgs[3]
    
    #factors_file <- '/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt.only_non_peer'
    #peer_file <- '/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt'
    #out_dir <- '.'
    
    #factors_file <- '/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm_as_covariates/lin.txt'
    #peer_file <- '/scratch/users/paedugar/transcriptome_diversity/PEER_covariates/lin.txt'
    #out_prefix <- '/scratch/users/paedugar/transcriptome_diversity/PEER_covariates/lin.'
    
    factors <- read.table(factors_file, sep='\t', stringsAsFactors=F, header=T, check.names=F)
    peer <- read.table(peer_file, sep='\t', stringsAsFactors=F, header=T, check.names=F)
    
    factors <- factors[,colnames(factors) %in% colnames(peer) , drop=F]
    peer <- peer[,colnames(peer), drop=F]
    
    for(j in 1:nrow(factors)) {
        
        results <- list()
        current_factor <- rankitNormalize_vector(as.numeric(factors[j,-1]))
        current_factor_name <- factors[j,1]
        
        for (i in 1:nrow(peer)) {
            current_peer <- rankitNormalize_vector(as.numeric(peer[i,-1]))
            normalized <- resid(lm(current_peer~current_factor))
            results[[i]] <- normalized
        }
        
        results <- as.data.frame(do.call(rbind, results))
        results <- cbind(peer[,1, drop=F], results)
        colnames(results) <- colnames(peer)
        
        write.table(results, paste0(out_prefix, current_factor_name, '.normalized_peer.txt'), sep='\t', quote=F, row.names=F)
    }
    
}

main()
