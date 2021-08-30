library("readr")
library("edgeR")
source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    out_exp_mat <- cmdArgs[2]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/lin_counts.txt'
    #out_exp_mat <- './deleteme_package.txt'
   
    #####
    # Read counts
    exp_mat <- read_tsv(exp_mat_file)
    exp_mat <- get_tmm(exp_mat)
    
    write_tsv(exp_mat, out_exp_mat)
    
    
}

main()
