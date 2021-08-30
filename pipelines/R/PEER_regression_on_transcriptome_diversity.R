# Calculates how much variance of transcriptome diversity can PEER covariates explain

library('readr')
library('dplyr')
library('tidyr')
library('broom')
source('../../R/transcriptome_diversity_tools.R', chdir=T)

main <- function(cmdArgs=commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    covariate_file <- cmdArgs[2]
    out_file <- cmdArgs[3]
    
    transcriptome_diversity_file <- '/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm/lin.txt'
    covariate_file <- '/scratch/users/paedugar/transcriptome_diversity/PEER_covariates/lin.txt'
    
    # Read trasncriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file)
    transcriptome_diversity <- transcriptome_diversity %>%
        mutate(transcriptome_diversity=rankitNormalize_vector(transcriptome_diversity),
               ind=sample_id) %>%
        select(-sample_id)
    
    # Read PEER
    peer <-  as.data.frame(t(read.table(covariate_file, sep='\t', stringsAsFactors=F, header=T, row.names=1, check.names=F)))
    
    peer <- as.data.frame(rankitNormalize(as.matrix(peer), 2))
    peer <- select(peer, starts_with('InferredCov'))
    peer$ind <- rownames(peer)
    
    # Correct for gtex singularity
    if(any(grepl('GTEX', transcriptome_diversity$ind))) {
        
        if(all(!transcriptome_diversity$ind %in% peer$ind))
            transcriptome_diversity$ind <- gtexLongToShort(ind)
        
    }
    
    # Merge
    merged <- inner_join(transcriptome_diversity, peer)
    
    # Do linear regressions
    lm_result <- lm(transcriptome_diversity ~ . , data=select(merged, -ind))
    lm_result_tidy <- tidy(lm_result)
    lm_result_tidy$model_r_squared <- summary(lm_result)[['r.squared']]
    
    # Write results
    write_tsv(lm_result_tidy, out_file)

    
}

# convert long gtex ids to short ones 
gtexLongToShort <- function(x) {
     gsub('(GTEX-\\w+?)-.+', '\\1', x)
}


main()
