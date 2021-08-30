# Takes the big expression matrix from gtex v8 and its annotation file and creates individual expression matrices

library('readr')
library('dplyr')

main <- function(cmdArgs = commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    sample_annotation_file <- cmdArgs[2]
    out_sample_table_file <- cmdArgs[3]
    root_dir <- cmdArgs[4]
    relative_dir <- cmdArgs[5]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
    #sample_annotation_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/GTEx_Analysis_sample_annotation_file_v8.txt'
    
    # Read sample annotation and select RNA-seq
    sample_annotation <- read_tsv(sample_annotation_file)
    sample_annotation <- sample_annotation %>%
        filter(SMAFRZE == 'RNASEQ') %>%
        select(SAMPID, SMTSD) %>%
        mutate(SMTSD = code_friendly_tissues(SMTSD))
    
    
    # Get sample table ready
    sample_table <- data.frame(id=paste0('gtex_', unique(sample_annotation$SMTSD)))
    sample_table$name <- sample_table$id
    sample_table$path <- file.path(relative_dir, paste0(sample_table$id, '.txt'))
    sample_table$group <- 'gtex'
    sample_table$subgroup <- 'gtex'
    sample_table$method <- 'gtex_original'
    sample_table$count_type <- 'TPM'
    sample_table$tissue <- unique(sample_annotation$SMTSD)
    
    write_tsv(sample_table, out_sample_table_file)
    
    # Extracting expression values
    cat('Reading expression matrix \n')
    exp_mat <- read_tsv(exp_mat_file, skip=2)
    
    for(tissue in unique(sample_annotation$SMTSD)) {
       
        cat('Working with ', tissue, '\n')
        current_samples <- filter(sample_annotation, SMTSD == tissue)
        current_samples_index <- colnames(exp_mat) %in% current_samples$SAMPID
        current_samples_index[1] <- T
        current_exp_mat <- exp_mat[,current_samples_index]
        
        write_tsv(current_exp_mat, file.path(root_dir, relative_dir, paste0('gtex_', tissue, '.txt')))
        
    }
    
}

code_friendly_tissues <- function(x) {
    
    x <- gsub('-', '', x)
    x <- gsub('\\(', '', x)
    x <- gsub('\\)', '', x)
    x <- gsub('\\s+', '_', x)
    
    x[x=='Brain_Spinal_cord_cervical_c1'] <- 'Brain_Spinal_cord_cervical_c-1'
    x[x=='Cells_EBVtransformed_lymphocytes'] <- 'Cells_EBV-transformed_lymphocytes'
    
    return(x)
}

main()
