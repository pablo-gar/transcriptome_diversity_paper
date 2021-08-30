# Calculate how much of the variance from gene expression can be explained by transcriptome diversity

library('readr')
library('dplyr')
library('tidyr')
source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    exp_mat_file <- cmdArgs[2]
    out_pca_matrix <- cmdArgs[3]
    out_pca_stats <- cmdArgs[4]
    out_pca_cors <- cmdArgs[5]
    
    #transcriptome_diversity_file <- '/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm/gtex_Whole_Blood.txt'
    #exp_mat_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_matrices/gtex_Whole_Blood.txt'
    #out_pca_matrix <-'' 
    #out_pca_stats <- '' 
    #out_pca_cors <- '' 
    
    # Read transcriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file)
    #transcriptome_diversity$transcriptome_diversity <- rankitNormalize_vector(transcriptome_diversity$transcriptome_diversity)
    
    # Read expression
    exp_mat <- read_expression(exp_mat_file)
    
    # Filter lowly expressed genes
    exp_mat <- filter_rarely_expressed_genes(exp_mat, percentage=0.2, gene_names_as_rownames=F, count_cutoff=1)
    
    # Performs PCA
    pca_results <- pca_expression(exp_mat)
    pca_stats <- as.matrix(summary(pca_results)$importance)
    pc_matrix <- as.data.frame(pca_results$x)
    
    # Correlating PCs with transcriptome_diversity
    pc_matrix$sample_id <- rownames(pc_matrix)
    pc_matrix <- pivot_longer(pc_matrix, -sample_id, names_to='PC', values_to='PC_value')
    
    # Merge
    pc_matrix <- inner_join(transcriptome_diversity, pc_matrix)

    cors <- pc_matrix %>%
        group_by(PC) %>%
        mutate(transcriptome_diversity = rankitNormalize_vector(transcriptome_diversity), 
               PC_value = rankitNormalize_vector(PC_value)) %>%
        summarise(pearson_cor = cor_test(transcriptome_diversity, PC_value, val='estimate'),
                  pvalue = cor_test(transcriptome_diversity, PC_value, val='p.value'),
                  signed_pvalue = sign(pearson_cor) * -log10(pvalue)) %>%
        ungroup() %>%
        filter(!is.na(pearson_cor)) %>%
        mutate(pvalue_bonf = p.adjust(pvalue)) %>%
        arrange(desc(pearson_cor))
        
    # Write results
    pca_stats <- as.data.frame(pca_stats)
    pca_stats$stat <- rownames(pca_stats)
    write_tsv(pc_matrix, out_pca_matrix)
    write_tsv(cors, out_pca_cors)
    write_tsv(pca_stats, out_pca_stats)
}

main()
