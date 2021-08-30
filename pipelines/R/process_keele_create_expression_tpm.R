library('readr')

main <- function(cmdArgs=commandArgs(T)) {
    
    matrix_dir <- cmdArgs[1]
    prefix <- cmdArgs[2]
    out_dir <- cmdArgs[3]
    
    #matrix_dir <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/GEO_matrices'
    #prefix <- 'keele_tpm_'
    #out_dir <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices'
    
    matrix_files <- list.files(matrix_dir, full.names=T)
    info <- strsplit(basename(matrix_files), '_')
    
    samples <- data.frame(strain=sapply(info, function(x) x[2]),
                          tissue=sapply(info, function(x) x[3]),
                          file_c=matrix_files,
                          stringsAsFactors=F
                          )
    
    for(tissue in unique(samples$tissue)) {
        
        current_samples <- samples[samples$tissue==tissue,]
        current_exp_mat <- list()
        
        for(i in 1:nrow(current_samples)) {
            x <- read_tsv(current_samples[i,'file_c'])
            current_exp_mat[[i]] <- x[,c('gene_id', 'TPM')]
            colnames(current_exp_mat[[i]])[2] <- current_samples[i,'strain']
        }
        
        current_exp_mat <- Reduce(function(...) merge(..., by='gene_id', all.x=TRUE), current_exp_mat)
        write_tsv(current_exp_mat, file.path(out_dir, paste0(prefix, tissue, '.txt')))
        
    }
    
    
}

main()
