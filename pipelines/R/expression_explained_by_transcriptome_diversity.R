# Calculate how much of the variance from gene expression can be explained by transcriptome diversity

library('readr')
library('dplyr')
library('tidyr')
library('ggplot2')
source('../../R/gtex.R', chdir=T)
source('../../R/misc.R', chdir=T)
source('../../R/ggthemes.R', chdir=T)

main <- function(cmdArgs=commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    gtex_version <- cmdArgs[2]
    out_pca_matrix <- cmdArgs[3]
    out_pca_stats <- cmdArgs[4]
    out_pca_cors <- cmdArgs[5]
    
    if(!gtex_version %in% c('v7', 'v8'))
        stop('Gtex version has to be v7 or v8 ')
    
    #transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/transcriptome_diversity/Whole_Blood.txt'
    #tissue <- 'Whole_Blood'
    #gtex_version <- 'v8'
    #out_pca_matrix <-'/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/PCA_var_explained_by_transcriptome_diversity/Whole_Blood/pca_matrix.txt' 
    #out_pca_stats <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/PCA_var_explained_by_transcriptome_diversity/Whole_Blood/pca_stats.txt' 
    #out_pca_cors <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/PCA_var_explained_by_transcriptome_diversity/Whole_Blood/cors.txt' 
    
    # Read transcriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file)
    transcriptome_diversity$transcriptome_diversity <- rankitNormalize_vector(transcriptome_diversity$transcriptome_diversity)
    transcriptome_diversity$ind <- gtexLongToShort(transcriptome_diversity$gtexId)
    
    # Read expression
    if(gtex_version=='v8') {
        exp_mat <- readAllGtexExpression(transcriptome_diversity$gtexId,  GTEX_CON$expresionAllTissuesV8)
    } else {
        exp_mat <- readAllGtexExpression(transcriptome_diversity$gtexId)
    }
    
    gene_ids <- gsub('\\..+', '', exp_mat[,1])
    exp_mat <- exp_mat[!duplicated(gene_ids), ]
    rownames(exp_mat) <- gsub('\\..+', '', exp_mat[,1])
    exp_mat <- exp_mat[,-1]
    
    # Select good expression
    exp_mat <- exp_mat[rowSums(exp_mat > 1) >= (ncol(exp_mat) * 0.2),]
    exp_mat <- rankitNormalize(as.matrix(exp_mat))
    exp_mat <- t(exp_mat)
    
    
    # Performs PCA
    pca_results <- prcomp(exp_mat)
    pca_stats <- as.matrix(summary(pca_results)$importance)
    peer <- as.data.frame(pca_results$x)
    
    # Correlating PCs with transcriptome_diversity
    peer$gtexId <- rownames(peer)
    peer <- pivot_longer(peer, -gtexId, names_to='feature', values_to='feature_value')
    
    # Merge
    merged <- inner_join(transcriptome_diversity, peer)

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
    
    
    # Plot 
    toPlot <- as.data.frame(pca_results$x)
    toPlot$gtexId <- rownames(toPlot)
    toPlot <- toPlot %>%
        select(c(gtexId,PC1, PC2)) %>%
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
    pca_stats <- as.data.frame(pca_stats)
    pca_stats$stat <- rownames(pca_stats)
    write_tsv(peer, out_pca_matrix)
    write_tsv(cors, out_pca_cors)
    write_tsv(pca_stats, out_pca_stats)
    ggsave(file.path(dirname(out_pca_cors), 'PC1_PC2_transcriptome_diversity.pdf'), p)

    
}

main()
