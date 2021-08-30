library("readr")
library("dplyr")

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat2_file <- cmdArgs[1]
    cov_file <- cmdArgs[2]
    out_exp_mat <- cmdArgs[3]
    out_cov_mat <- cmdArgs[4]
    
    #exp_mat2_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/lin/GSE60314_5_57_HTSeq_raw_read_counts.txt.gz'
    #cov_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/lin/GSE60314_GEO_run_summary_metadata.txt'
    #gtf_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/lin/GSE60314_dme5_57_ERCC.gtf.gz'
    #out_exp_mat <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/lin_counts.txt'
    #out_cov_mat <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/lin_covariates.txt'
    
    ##########################
    # Reading metadata_lin
    metadata_lin <- read_tsv(cov_file)
    cols <- c("Library_ID", "DGRP_Number", "Environment", "Fly_Number", "Sex", "RNA_Prep_Method", "mapped_reads")
    
    # Getting mapped reads 
    read_count <- read_tsv(cov_file) %>%
        group_by(Library_ID) %>%
        summarise(mapped_reads=sum(Mapped_Reads_5.57)) %>%
        ungroup()
    
    metadata_lin <- left_join(metadata_lin, read_count)
    
    metadata_lin <- unique(as.data.frame(metadata_lin[,cols]))
    rownames(metadata_lin) <- metadata_lin[,1]
    metadata_lin <- t(metadata_lin[,-1])
    metadata_lin <- cbind(data.frame(ID=rownames(metadata_lin), stringsAsFactors=F), metadata_lin)
    
    
    #####
    # Read counts
    exp_mat <- read_tsv(exp_mat2_file)
    exp_mat <- exp_mat[,-2]
    colnames(exp_mat)[1] <- 'gene_id'
    
    # Read gene annotation for gene lengths
    #gtf <- makeTxDbFromGFF(gtf_file, format='gtf', organism="Drosophila melanogaster")
    #gene_lengths <- exonsBy(gtf, 'gene')
    #gene_lengths <- sapply(gene_lengths, function(x) sum(width(reduce(x)))) 
    #gene_lengths <- gene_lengths[exp_mat$X1]
    #
    #TPM <- calculateTPM(as.matrix(exp_mat[,-(1:2)]), gene_lengths)
    
    write_tsv(metadata_lin, out_cov_mat)
    write_tsv(exp_mat, out_exp_mat)
    
    
}

main()
