proj_dir='/scratch/users/paedugar/transcriptome_diversity/'

Rscript ./get_transcriptome_diversity_as_covariates_lin.R $proj_dir/transcriptome_diversity/tpm/lin.txt $proj_dir/transcriptome_diversity/tpm_as_covariates/lin.txt

Rscript ./regress_out_from_peer.R $proj_dir/transcriptome_diversity/tpm_as_covariates/lin.txt $proj_dir/PEER_covariates/lin.txt $proj_dir/PEER_covariates/lin.

# TPM
Rscript ./fit_covariates_to_expression.R $proj_dir/expression_matrices_processed/tpm/lin.txt $proj_dir/PEER_covariates/lin.txt $proj_dir/PEER_analyses/var_explained_per_gene/inhouse_PEER/tpm/lin/
Rscript ./fit_covariates_to_expression.R $proj_dir/expression_matrices_processed/tpm/lin.txt $proj_dir/PEER_covariates/lin.transcriptome_diversity.normalized_peer.txt $proj_dir/PEER_analyses/var_explained_per_gene/inhouse_PEER/tpm/lin

# TMM
Rscript ./fit_covariates_to_expression.R $proj_dir/expression_matrices_processed/tmm/lin.txt $proj_dir/PEER_covariates/lin.txt $proj_dir/PEER_analyses/var_explained_per_gene/inhouse_PEER/tmm/lin/
Rscript ./fit_covariates_to_expression.R $proj_dir/expression_matrices_processed/tmm/lin.txt $proj_dir/PEER_covariates/lin.transcriptome_diversity.normalized_peer.txt $proj_dir/PEER_analyses/var_explained_per_gene/inhouse_PEER/tmm/lin

