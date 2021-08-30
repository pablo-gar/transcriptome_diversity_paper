source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    n_reads <- as.numeric(cmdArgs[2])
    out_exp_mat <- cmdArgs[3]
    
    #exp_mat_file <- 'lin_counts_small.txt'
    #n_reads <- 5e6
    
    # Read counts
    exp_mat <- read_expression(exp_mat_file)
    
    # Sample
    exp_mat <- resize_expression(exp_mat, n=n_reads, method='montecarlo')
    
    #Read gene annotation for gene lengths
    write.table(exp_mat, out_exp_mat, sep='\t', quote=F, row.names=F, col.names=T)
    
    
}

main()
