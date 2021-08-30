# Calculates correlation between the transcriptome diversity and Principal Components of expression

source("../../R/transcriptome_diversity_tools.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
library('tidyr')
library('dplyr')
library('ggplot2')
library('readr')

main <- function(cmdArgs = commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    transcriptome_diversity_file <- cmdArgs[2]
    out_pca_matrix <- cmdArgs[3]
    out_pca_stats <- cmdArgs[4]
    out_pca_cors <- cmdArgs[5]
    
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/expression_associations_from_raw/tmm/lin.txt'
    #transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/expression_associations_from_raw/transcriptome_diversity/lin.txt'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/arora_TCGA_recount2_TPM_TCGATHCA.txt'
    #exp_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/arora_TCGA_xena_RPKM_BREAST_INVASIVE_CARCINOMA.txt'
    
    # Read and get transcriptome diversity
    exp_mat <- read_expression(exp_mat_file)
    transcriptome_diversity <- read.table(transcriptome_diversity_file, sep='\t', stringsAsFactors=F, header=T)
    transcriptome_diversity$transcriptome_diversity <- rankitNormalize_vector(transcriptome_diversity$transcriptome_diversity)
    
    # Filter exp for PCA
    exp_mat2 <- filter_rarely_expressed_genes(exp_mat)
    exp_mat2 <- rankit_normalize_expression(exp_mat2)
    
    pca_results <- pca_expression(exp_mat2)
    
    pca_stats <- as.matrix(summary(pca_results)$importance)
    pca_stats <- as.data.frame(pca_stats)
    pca_stats$stat <- rownames(pca_stats)
    
    pca_loadings <- as.data.frame(pca_results$x)
    pca_loadings$sample_id <- rownames(pca_loadings)
    
    # Get correlations between transcriptome_diversity and pc loadings
    cors <- cor_pca_transcriptome_diverisity(transcriptome_diversity, pca_loadings)
    
    # append transcriptome_diversity to pca loading
    pca_loadings <- left_join(pca_loadings, transcriptome_diversity)
    
    # Plot 
    toPlot <- pca_loadings
    toPlot <- toPlot %>%
        select(c(sample_id,PC1, PC2)) %>%
        inner_join(transcriptome_diversity)
    
    p <- ggplot(toPlot, aes(x=PC1, y=PC2)) +
            geom_point(aes(colour=transcriptome_diversity), size=1.8) +
            scale_colour_gradient(low = "yellow", high = "red") +
            labs(subtitle=paste0('PC1 pearson correlation with transcriptome diversity:\nr=', signif(cors[cors$feature=='PC1', 'pearson_cor'], 2))) +
            xlab(paste0('PC1 - var explained=', pca_stats['Proportion of Variance',1])) +
            ylab(paste0('PC2 - var explained=', pca_stats['Proportion of Variance',2])) +
            theme_sleek() +
            theme(legend.position='top')
    
    # Write results
    write_tsv(pca_loadings, out_pca_matrix)
    write_tsv(cors, out_pca_cors)
    write_tsv(pca_stats, out_pca_stats)
    ggsave(file.path(dirname(out_pca_cors), 'PC1_PC2_transcriptome_diversity.pdf'), p)
}

cor_pca_transcriptome_diverisity <- function(transcriptome_diversity, pca_loadings) {
    
    pca_loadings$sample_id <- rownames(pca_loadings)
    pca_loadings <- pivot_longer(pca_loadings, -sample_id, names_to='feature', values_to='feature_value')
    
    # Merge
    merged <- inner_join(transcriptome_diversity, pca_loadings)

    cors <- merged %>%
        group_by(feature) %>%
        mutate(transcriptome_diversity = rankitNormalize_vector(transcriptome_diversity), 
               feature_value = rankitNormalize_vector(feature_value)) %>%
        summarise(pearson_cor = cor_test(transcriptome_diversity, feature_value, val='estimate'),
                  pvalue = cor_test(transcriptome_diversity, feature_value, val='p.value'),
                  signed_pvalue = sign(pearson_cor) * -log10(pvalue)) %>%
        ungroup() %>%
        filter(!is.na(pearson_cor)) %>%
        mutate(pvalue_bonf = p.adjust(pvalue)) %>%
        arrange(desc(pearson_cor))
    
    return(cors)
    
}

main()
