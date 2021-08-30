source('../../R/misc.R')
library('readr')
library('stringr')

main <- function(cmdArgs=commandArgs(T)) {
    
    eqtl_trans_folder <- cmdArgs[1] # Covariates: transcriptome diversity
    eqtl_peer_folder <- cmdArgs[2] # Covariates: PEER covariates
    eqtl_no_peer_folder <- cmdArgs[3] # Covariates: non-peer covariates
    eqtl_peer_plus_folder <- cmdArgs[4] # Covariates: PEER covariates + non-peer covariates
    eqtl_no_cov_folder <- cmdArgs[5] # No covariates
    eqtl_trans_no_peer_folder <- cmdArgs[6] # Covariates: transcriptome diversity + non-peer covariates
    out_file <- cmdArgs[7]  
    
    eqtl_trans_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_trans_as_covariate"
    eqtl_peer_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_PEER_covariate"
    eqtl_no_peer_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_nonPEER_covariate"
    eqtl_peer_plus_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_full_covariate"
    eqtl_no_cov_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_no_covariate"
    eqtl_trans_no_peer_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_trans_as_covariate_plus_non_peer"
    eqtl_peer_regressed_out_folder <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/eqtls_peer_covariate_trasn_div_regressed_out"
    out_file <- "/scratch/users/paedugar/transcriptome_diversity/eqtl_results/gtex/all_results/gtex_eqtl_results_new.txt"
    
    
    # Read files
    eqtl_trans <- read_eqtl(eqtl_trans_folder, 'transcriptome_diversity_as_covariates')
    eqtl_trans_no_peer <- read_eqtl(eqtl_trans_no_peer_folder, 'transcriptome_diversity_as_covariates_plus_nonPeer')
    eqtl_peer <- read_eqtl(eqtl_peer_folder, 'peer')
    eqtl_no_peer <- read_eqtl(eqtl_no_peer_folder, 'no_peer')
    eqtl_no_trans <- read_eqtl(eqtl_no_cov_folder, 'no_covariates')
    eqtl_peer_plus <- read_eqtl(eqtl_peer_plus_folder, 'peer_plus')
    eqtl_peer_regressed <- read_eqtl(eqtl_peer_regressed_out_folder, 'peer_regressed')

    # Bind results
    eqtl_results <- bind_rows(eqtl_trans, eqtl_trans_no_peer, eqtl_no_trans, eqtl_peer, eqtl_no_peer, eqtl_peer_plus, eqtl_peer_regressed)
    eqtl_results$group <- factor(eqtl_results$group, levels=c('no_covariates', 'transcriptome_diversity_as_covariates', 'no_peer', 'transcriptome_diversity_as_covariates_plus_nonPeer', 'peer', 'peer_plus', 'peer_regressed'))
    colnames(eqtl_results)[1:10] <- c('gene_id', 'variants_tested', 'MLE_1_beta_shape', 'MLE_2_beta_shape', 'dummy', 'variant_id', 'ditance_to_tss', 'nominal_pvalue', 'permutation1_pvalue', 'permutation_pvalue')

    # Calculate FDR 
    eqtl_results <- eqtl_results %>%
        group_by(group) %>%
        mutate(fdr=p.adjust(permutation_pvalue, method='BH')) %>%
        ungroup()
    
    write_tsv(eqtl_results, out_file)
    
}

read_eqtl <- function(eqtl_dir, group) {
    
    files <- list.files(eqtl_dir, full.names=T, recursive=T, pattern='txt.gz$')
    ids <- str_split_fixed(basename(files), '\\.', n=2)[,1]
    cors <- concatenate_table_files2(files, id_names=ids, delim=' ', header=F, read_function=read_delim, progress=F, col=cols(), skip_err=T)
    cors$group <- group
    
    return(cors)
}

main()


