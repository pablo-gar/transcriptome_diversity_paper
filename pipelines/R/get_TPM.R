library("GenomicFeatures")
library("readr")
source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    gtf_file <- cmdArgs[2]
    out_exp_mat <- cmdArgs[3]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/lin_counts.txt'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/lin/GSE60314_dme5_57_ERCC.gtf.gz'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/sarantopoulou_v4.txt'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/gencode.vM9.annotation.gtf.gz'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/Mus_musculus.NCBIM37.67.gtf.gz'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan.txt'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/gtex_Lung_counts.txt'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/gtex/gencode.v26.GRCh38.genes.gtf'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/keele_lung.txt'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/Mus_musculus.NCBIM37.67.gtf.gz'
    #out_exp_mat <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan_tpm.txt'
    #out_exp_mat <- './deleteme_package.txt'
   
    #####
    # Read counts
    exp_mat <- read_tsv(exp_mat_file)
    
    #Read gene annotation for gene lengths
    gtf <- makeTxDbFromGFF(gtf_file, format='gtf')
    gene_lengths <- exonsBy(gtf, 'gene')
    gene_lengths <- sapply(gene_lengths, function(x) sum(width(reduce(x)))) 
    
    if(grepl('ENS', names(gene_lengths[1]))) {
        if(!any(exp_mat[,1, drop=T] %in% names(gene_lengths))) {
            new_names <- gsub('\\..+$', '', names(gene_lengths))
            gene_lengths <- gene_lengths[!duplicated(new_names)]
            names(gene_lengths) <- new_names[!duplicated(new_names)]
        }
    }
        
    
    exp_mat <- exp_mat[exp_mat[,1,drop=T] %in% names(gene_lengths),]
    exp_mat <- get_tpm(exp_mat, gene_lengths)
    
    write_tsv(exp_mat, out_exp_mat)
    
    
}

main()
