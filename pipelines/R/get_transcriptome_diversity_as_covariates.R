# re formats transcriptome diversity file to square wide format
# Usage

library('readr')
library('tidyr')
library('dplyr')
source('../../R/transcriptome_diversity_tools.R', chdir=T)


main <- function(cmdArgs = commandArgs(T)) {
    
    transcriptome_diversity_file <- cmdArgs[1]
    bed_phenotype_file <- cmdArgs[2]
    out_file <- cmdArgs[3]
    
    #transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/transcriptome_diversity_all/Whole_Blood.txt'
    #bed_phenotype_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz'
    
    # Read ordering of samples
    sample_order <- unlist(strsplit(readLines(bed_phenotype_file, n=1), '\t'))
    sample_order <- sample_order[grep('GTEX', sample_order)]
    
    # Read transcriptome diversity
    transcriptome_diversity <- read_tsv(transcriptome_diversity_file) %>%
        mutate(sample_id=gtexLongToShort(sample_id)) %>%
        filter(!duplicated(sample_id)) %>%
        mutate(transcriptome_diversity=rankitNormalize_vector(transcriptome_diversity)) %>%
        pivot_wider(names_from=sample_id, values_from=transcriptome_diversity)
    
    transcriptome_diversity <- transcriptome_diversity[,sample_order]
    
    transcriptome_diversity$ID <- 'transcriptome_diversity'
    transcriptome_diversity <- transcriptome_diversity[,c(ncol(transcriptome_diversity), 1:(ncol(transcriptome_diversity)-1))]
    
    write_tsv(transcriptome_diversity, out_file)
}

gtexLongToShort <- function(x) {
     gsub('(GTEX-\\w+?)-.+', '\\1', x)
}

main()
