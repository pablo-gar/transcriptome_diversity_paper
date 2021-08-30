# Calculates shannon diversity for all samples in gtex
# Usage
# Rscript get_transcriptome_diversity.R out_table.txt

source("../../R/transcriptome_diversity_tools.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/tcga_expression/TCGA-BRCA.htseq_fpkm-uq.tsv.gz'
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/transcriptome_diversity/exonJunctions/Whole_Blood.txt'
    
    exp_mat <- read_expression(exp_mat_file)
    transcriptome_diversity <- get_transcriptome_diversity(exp_mat)
    
    write.table(transcriptome_diversity, out_file, sep='\t', quote=F, row.names=F)
    
}

main()
