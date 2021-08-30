# Create a table of variance explained by transcriptome diversity of the gene expression matrix using the PCA method

library('tidyr')
library('dplyr')
source('../../R/misc.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    
    path <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    path <- "/scratch/users/paedugar/transcriptome_diversity/var_explained_by_transcriptome_diversity/pca_method/from_raw"
    out_file <- "../figures/data/supp_tables/variance_explained_by_trans_diversity.tsv"
    
    types <- c('tpm', 'tmm')
    
    results <- list()
    for(type in types) {
        files <- list.files(file.path(path, type), recursive=T, full.names=T, pattern='var_explained.txt')
        results[[type]] <- concatenate_table_files2(files, id_names = rep(type, length(files)))
        results[[type]]$dataset <- basename(dirname(files))
    }
    
    results <- do.call(rbind, results) %>%
        pivot_wider(id_cols=dataset, names_from=id_names, values_from=var_explained) %>%
        filter(!dataset %in% c('keele_kidney', 'keele_liver', 'keele_lung', 'nagarajan', 'sarantopoulou_pico', 'sarantopoulou_truseq', 'sarantopoulou_v4'))
    
    write_tsv(results, out_file)
    
}
