#Finds PEER covariates associated with transcriptome diversity

library('readr')
library('dplyr')
library('tidyr')
source('../../R/transcriptome_diversity_tools.R')

main <- function(cmdArgs=commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    covariate_file <- cmdArgs[2]
    cor_pvalue_cutoff <- as.numeric(cmdArgs[3])
    out_cors <- cmdArgs[4]
    out_stats <- cmdArgs[5]
    out_trans_PEER <- cmdArgs[6]
    
    
    # Read trasncriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file)
    transcriptome_diversity$transcriptome_diversity_rankit <- rankitNormalize_vector(transcriptome_diversity$transcriptome_diversity)
    transcriptome_diversity$ind <- gtexLongToShort(transcriptome_diversity$sample_id)
    
    # Read PEER
    peer <-  as.data.frame(t(read.table(covariate_file, sep='\t', stringsAsFactors=F, header=T, row.names=1, check.names=F)))
    
    peer <- as.data.frame(rankitNormalize(as.matrix(peer), 2))
    peer <- select(peer, starts_with('InferredCov'))
    peer$ind <- rownames(peer)
    peer <- pivot_longer(peer, -ind, names_to='feature', values_to='feature_value')
    
    # Merge
    merged <- inner_join(transcriptome_diversity, peer)

    cors <- merged %>%
        group_by(feature) %>%
        summarise(pearson_cor = cor_test(transcriptome_diversity_rankit, feature_value, method='spearman', val='estimate'),
                  pvalue = cor_test(transcriptome_diversity_rankit, feature_value, method='spearman', val='p.value'),
                  signed_pvalue = sign(pearson_cor) * -log10(pvalue)) %>%
        ungroup() %>%
        filter(!is.na(pearson_cor)) %>%
        mutate(pvalue_bonf = p.adjust(pvalue), fdr = p.adjust(pvalue, method='BH'))
    
    stats <- data.frame(total_peer=length(unique(peer$feature)), cor_with_transcriptome_diversity_peer=nrow(cors), cor_pvalue_cutoff=cor_pvalue_cutoff)
    
    # Save results
    write_tsv(cors, out_cors)
    write_tsv(stats, out_stats)
    write_tsv(merged, out_trans_PEER)
    
    
}

# convert long gtex ids to short ones 
gtexLongToShort <- function(x) {
     gsub('(GTEX-\\w+?)-.+', '\\1', x)
}

main()
