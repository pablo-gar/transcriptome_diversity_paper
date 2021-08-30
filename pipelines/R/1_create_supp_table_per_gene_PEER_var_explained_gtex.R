
library('tidyr')
library('dplyr')
source('../../R/misc.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    
    path <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    path <- "/scratch/users/paedugar/transcriptome_diversity/PEER_analyses/var_explained_per_gene"
    out_file <- "../figures/data/supp_tables/per_gene_PEER_var_explained.tsv"
    
    # Get lin data
    files <- list.files(file.path(path, 'inhouse_PEER'), recursive=T, full.names=T)
    lin <- concatenate_table_files2(files, id_names=files)
    lin$type <- basename(dirname(dirname(lin$id_names)))
    lin$PEER_controlled <- ifelse(grepl('lin.transcriptome_diversity.normalized_peer.txt', lin$id_names), 'transcriptome_diversity', 'intact')
    lin <- lin %>%
        group_by(type, PEER_controlled) %>%
        summarise(median_r_sqr = median(r_squared))
    lin$dataset <- 'lin'
    
    
    # get gtex data
    files <- list.files(file.path(path, 'original_PEER'), recursive=T, full.names=T)
    tissues <- unique(gsub('.*(gtex_.+?)\\/.*', '\\1', files))
    
    results <- list()
    for(tissue in tissues) {
        
        current_files <- files[grepl(tissue, files)]
        intact_files <- current_files[grep("v8.covariates.txt.only_peer", current_files)]
        other_files <- current_files[!current_files %in% intact_files]
        
        intact <- concatenate_table_files2(intact_files, id_names=intact_files)
        current <- concatenate_table_files2(other_files, id_names=other_files)
        
        intact$PEER_controlled <- 'intact'
        intact$type <- basename(dirname(dirname(intact$id_names)))
        
        current$PEER_controlled <- gsub('(.+?)\\..+', '\\1', basename(current$id_names))
        current$type <- basename(dirname(dirname(dirname(current$id_names))))
        
        current <- bind_rows(intact, current)
        
        current <- current %>%
            group_by(type, PEER_controlled) %>%
            summarise(median_r_sqr = median(r_squared))
        
        current$dataset <- tissue
        results[[tissue]] <- current
        
    }
    
    results <- do.call(bind_rows, results)
    
    
    # Combine data and rearrange
    merged <- bind_rows(lin, results)
    #merged <- results %>%
    merged <- merged %>%
        pivot_wider(id_cols=c(type, dataset), values_from=median_r_sqr, names_from=PEER_controlled) %>%
        mutate(transcriptome_diversity_of_PEER=(intact-transcriptome_diversity)/intact)
    
    write_tsv(merged, out_file)
}

main()
