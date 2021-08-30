# Calculates shannon diversity for all samples in gtex
# Usage
# Rscript get_transcriptome_diversity.R out_table.txt

source("../../R/transcriptome_diversity_tools.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/transcriptome_diversity/exonJunctions/Whole_Blood.txt'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/lin.txt'
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/expression_associations_from_raw/transcriptome_diversity_non_lowly_expressed/lin_expressed_in_all.txt'
    
    exp_mat <- read_expression(exp_mat_file)
    transcriptome_diversity <- get_transcriptome_diversity(exp_mat, count_cutoff=5, expressed_in_all=T)
    
    write.table(transcriptome_diversity, out_file, sep='\t', quote=F, row.names=F)
    
}

main()
