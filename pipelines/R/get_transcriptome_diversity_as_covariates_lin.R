# re formats transcriptome diversity file to square wide format
# Usage

library('readr')
library('tidyr')
library('dplyr')
source('../../R/misc.R', chdir=T)


main <- function(cmdArgs = commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    out_file <- cmdArgs[2]
    
    #transcriptome_diversity_file <- '/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm/lin.txt'
    #out_file <- '/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm_as_covariates/lin.txt'
    
    # Read ordering of samples
    
    # Read transcriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file) %>%
        filter(!duplicated(sample_id)) %>%
        mutate(transcriptome_diversity=rankitNormalize_vector(transcriptome_diversity)) %>%
        pivot_wider(names_from=sample_id, values_from=transcriptome_diversity)
    
    transcriptome_diversity$ID <- 'transcriptome_diversity'
    transcriptome_diversity <- transcriptome_diversity[,c(ncol(transcriptome_diversity), 1:(ncol(transcriptome_diversity)-1))]
    
    write_tsv(transcriptome_diversity, out_file)
}

main()
