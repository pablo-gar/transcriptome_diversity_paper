library('tidyr')
library('readr')
library('dplyr')
source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    genome_file <- cmdArgs[2]
    converstion_table <- cmdArgs[3]
    out_exp_file <- cmdArgs[4]
    out_cov_file <- cmdArgs[5]
    out_bed_file <- cmdArgs[6]
    
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/kidney_expression.csv.zip'
    #genome_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/refseq_mm9_tss.txt.zip'
    #out_exp_file <- 'test_exp.txt'
    #out_bed_file <- 'test_exp_bed.txt.gz'
    #out_cov_file <- 'test_cov.txt'
    
    
    exp_mat <- read_csv(exp_mat_file)
    genome <- read_tsv(genome_file, col_names=F)
    genome <- select(genome, -X5)
    colnames(genome) <- c('#chr', 'start', 'end', 'gene_id')
    
    
    # Create covariate files
    cov_mat <- exp_mat %>%
        select(BATCH, SUBJECT.NAME) %>%
        unique() %>%
        mutate(BATCH=convert_batch_to_string(BATCH)) %>%
        pivot_wider(names_from=SUBJECT.NAME, values_from=BATCH) %>%
        mutate(ID='BATCH')
    
    cov_mat <- cov_mat[,c(ncol(cov_mat), 1:(ncol(cov_mat)-1))]
    
    # Create expression mat file
    exp_mat <- exp_mat %>%
        select(SUBJECT.NAME, starts_with('gene')) %>%
        t()
    
    exp_mat_subject <- exp_mat[1,]
    exp_mat_gene <- gsub('gene_', '', rownames(exp_mat)[-1])
    exp_mat <- as.data.frame(apply(as.matrix(exp_mat[-1,]), 2, as.numeric))
    colnames(exp_mat) <- exp_mat_subject
    exp_mat <- bind_cols(data.frame(gene_id=exp_mat_gene, stringsAsFactors=F), exp_mat)
    
    # Convert ids to ensembl
    conversion_vector <- read.table(converstion_table, sep='\t', stringsAsFactors=F, header=T)
    conversion_vector <- conversion_vector[!duplicated(conversion_vector[,2]),]
    conversion_vector <- setNames(conversion_vector[,1], conversion_vector[,2])
    
    exp_mat[,1] <- conversion_vector[exp_mat[,1]]
    exp_mat <- exp_mat[!is.na(exp_mat[,1]),]
    
    genome$gene_id <- conversion_vector[genome$gene_id]
    genome <- genome[!is.na(genome$gene_id),]
    
    # Create expression mat if bed format
    for_bed <- filter_rarely_expressed_genes(exp_mat)
    for_bed <- rankit_normalize_expression(for_bed)
    exp_mat_bed <- inner_join(genome, for_bed)
    
    write_tsv(exp_mat, out_exp_file)
    write_tsv(exp_mat_bed, out_bed_file)
    write_tsv(cov_mat, out_cov_file)
    
    
}

convert_batch_to_string <- function(x) {
    
    un <- unique(x)
    un <- setNames(paste0('batch', 1:length(un)), un)
    
    x <- un[as.character(x)]
    
    
    return(x)
    
}

main()
