library('readr')
library("GEOquery")
library("stringr")

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    out_exp_mat <- cmdArgs[2]
    out_cov <- cmdArgs[3]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/nagarajan/arid1a_knockout_COUNTS.tsv'
    #out_exp_mat <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan.txt'
    #out_cov <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/nagarajan_cov.txt'
    
    exp_mat <- read_tsv(exp_mat_file)
    colnames(exp_mat)[1] <- 'gene_id'
    exp_mat$gene_id <- gsub('\\..+$', '', exp_mat$gene_id)
    exp_mat <- exp_mat[!duplicated(exp_mat$gene_id), ]
    
    # Reading geo data to get replicate and clone information
    geo_id  <- "GSE123285"
    gsm <- getGEO(geo_id, GSEMatrix=T, parseCharacteristics=F)
    metadata <- as.character(pData(phenoData(gsm[[1]]))[,'title'])
    
    metadata <- gsub('Clone ', 'Clone-', metadata)
    metadata <- gsub('rep ', 'rep-', metadata)
    metadata <- str_split_fixed(metadata, " ", 3)
    metadata[metadata[,2] == '', 2] <- 'clone-NA'
    rownames(metadata) <- metadata[,1]
    
    # Bulding covariates
    covariates <- data.frame(sample_id=colnames(exp_mat)[-1], genotype='wt', lineage='parental', genotype_lineage='', rep='', treatment='vehicle', stringsAsFactors=F)
    covariates$clone <- metadata[covariates$sample_id, 2]
    covariates$rep <- metadata[covariates$sample_id, 3]
    covariates$genotype[grepl('arid', covariates$sample_id)] <- 'arid1a'
    covariates$lineage[!grepl('parental', covariates$sample_id)] <- 'final'
    covariates$genotype_lineage <- paste0(covariates$genotype, '_', covariates$lineage)
    covariates$treatment[grepl('jq1', covariates$sample_id)] <- 'jq1'
    covariates$treatment[grepl('tamoxifen', covariates$sample_id)] <- 'tamoxifen'
    covariates$treatment[grepl('fulvestrant', covariates$sample_id)] <- 'fulvestrant'
    
    covariates <- t(covariates)
    colnames(covariates) <- covariates[1,]
    covariates <- covariates[-1,]
    
    covariates <- cbind(data.frame(ID=rownames(covariates)), covariates)
    
    
    # Outputing tables
    write_tsv(covariates, out_cov)
    write_tsv(exp_mat, out_exp_mat)
    

    
}

main()
