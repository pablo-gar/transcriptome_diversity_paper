# Calculates correlation between the transcriptome diversity and Principal Components of expression

source("../../R/transcriptome_diversity_tools.R", chdir = T)
source("../../R/ggthemes.R", chdir = T)
library('tidyr')
library('dplyr')
library('ggplot2')
library('readr')

main <- function(cmdArgs = commandArgs(T)) {
    
    pca_loadings_file <- cmdArgs[1]
    pca_stats_file <- cmdArgs[2]
    covariate_file <- cmdArgs[3]
    out_file <- cmdArgs[4]
    
    #pca_stats_file <- '/scratch/users/paedugar/transcriptome_diversity/var_explained_by_transcriptome_diversity/pca_method/from_raw/tpm/lin/PC_stats.txt'
    #pca_loadings_file <- '/scratch/users/paedugar/transcriptome_diversity/var_explained_by_transcriptome_diversity/pca_method/from_raw/tpm/lin/PC_matrix.txt'
    #covariate_file <- '/scratch/users/paedugar/transcriptome_diversity/PEER_covariates/lin.txt'
    
    # Read and get transcriptome diversity
    pca_loadings <- read_tsv(pca_loadings_file)
    pca_stats <- read_tsv(pca_stats_file)
    
    pca_loadings <- pivot_longer(pca_loadings, contains("PC"), names_to='feature', values_to='feature_value')
    
    pca_stats <- pivot_longer(pca_stats, contains("PC"), names_to='feature', values_to='feature_value') %>%
        dplyr::filter(stat=='Proportion of Variance') %>%
        select(c(feature, feature_value))
    colnames(pca_stats)[2] <- 'pc_r_squared'
    
    # Read covariate file
    covariates <- read_tsv(covariate_file)
    covariates <- pivot_longer(covariates, -ID, names_to='sample_id', values_to='covariate_value')
    
    # Get correlations between transcriptome_diversity and pc loadings
    cors <- cor_pca_transcriptome_diverisity(covariates, pca_loadings)
    
    # Join both results
    cors_merged <- left_join(cors, pca_stats) %>%
        mutate(var_explained_per_pc=r_squared*pc_r_squared) %>%
        group_by(ID) %>%
        summarise(var_explained=sum(var_explained_per_pc)) %>%
        ungroup()
    
    
    # Write results
    write_tsv(cors_merged, out_file)
}

cor_pca_transcriptome_diverisity <- function(covariates, pca_loadings) {
    
    
    # Merge
    merged <- inner_join(covariates, pca_loadings)

    cors <- merged %>%
        group_by(feature,ID) %>%
        mutate(covariate_value = rankitNormalize_vector(covariate_value), 
               feature_value = rankitNormalize_vector(feature_value)) %>%
        summarise(pearson_cor = cor_test(covariate_value, feature_value, val='estimate'),
                  pvalue = cor_test(covariate_value, feature_value, val='p.value'),
                  signed_pvalue = sign(pearson_cor) * -log10(pvalue)) %>%
        ungroup() %>%
        #filter(!is.na(pearson_cor)) %>%
        mutate(pvalue_bonf = p.adjust(pvalue), r_squared = pearson_cor^2)
        #arrange(desc(pearson_cor))
    
    return(cors)
    
}

main()
