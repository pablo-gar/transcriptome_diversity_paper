# Calculates shannon diversity for all samples in gtex
# Usage
# Rscript get_transcriptome_diversity.R out_table.txt

source("../../R/transcriptome_diversity_tools.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/gtex_Whole_Blood.txt'
    #good_genes_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/Whole_Blood.v8.good_genes.txt'
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices_transcriptomeDiversity_normalized/gtex_Whole_Blood.txt'
    
    exp_mat <- read_expression(exp_mat_file)
    
    #good_genes <- readLines(good_genes_file)
    #exp_mat <- exp_mat[ exp_mat[,1] %in% good_genes,]
    
    exp_mat <- normalize_transcriptome_diversity(exp_mat)
    
    write.table(exp_mat, out_file, sep='\t', quote=F, row.names=F)
    
}

main()
