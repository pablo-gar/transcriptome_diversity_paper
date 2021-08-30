# To make expression tables first run:
# cd R; Rscript create_sample_table_DNA_samples.R # HAS NOT BEEN UPDATED 

import pandas as pd
import os
import re

configfile: '../config.json'
localrules: all

keele_tissues=['kidney', 'liver', 'lung']
#keele_tissues=['kidney']

#---------------------------------#
# Load sample table raw counts 
sample_table = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_for_expression_associations']))
sample_table = sample_table.set_index('id')

all_datasets = set(sample_table.index.get_level_values('id'))



#---------------------------------#
# Load sample table original matrices
sample_table_original = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files']))
sample_table_original = sample_table_original.set_index('id')

all_datasets_original = set(sample_table_original.index.get_level_values('id'))


# getting arora
arora_datasets = [i for i in list(all_datasets_original) if 'arora_' in i]

#---------------------------------#
# Load smaple table for gtex eqtl
sample_table_gtex_bed = pd.read_table(os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files_gtex_bed']))
sample_table_gtex_bed = sample_table_gtex_bed.set_index('group')
all_tissues_gtex = set(sample_table_gtex_bed.index.get_level_values('group'))
tissues_eqtl_gtex = set(['Whole_Blood', 'Muscle_Skeletal', 'Skin_Sun_Exposed_Lower_leg'])


all_datasets_subset_gtex = [i for i in all_datasets if 'gtex' in i]
all_datasets_subset_gtex_noPrefix = [re.sub('gtex_', '', i) for i in all_datasets_subset_gtex]
all_datasets_subset_gtex_noPrefix_temp = list(set(all_datasets_subset_gtex_noPrefix).intersection(set(all_tissues_gtex)))
all_datasets_subset_gtex_noPrefix = all_datasets_subset_gtex_noPrefix_temp + ['Cells_EBVtransformed_lymphocytes', 'Brain_Spinal_cord_cervical_c1']
all_datasets_subset_gtex = ['gtex_' + i for i in all_datasets_subset_gtex_noPrefix]



#---------------------------------#
# Pipeline
# 

rule all:
    input: 
        ## Transcriptome diversity from down-sampled reads
        expand(os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'sampled_reads', '{n_reads}', '{dataset}.txt'), n_reads=['2.5e5', '5e5', '1e6', '2.5e6', '5e6', '10e6'], dataset='lin'),
        #### Transcriptome diversity for Arora et al. data
        expand(os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_sources', '{dataset}.txt'),  dataset=arora_datasets),
        ### Transcriptome diversitye non_lowly expressed
        ### Get PEER covariates
         expand(os.path.join(config['projectDir'], config['PEER_dir'], '{dataset}.txt'),  dataset=list(all_datasets)),
        ### Get expression associations
        expand(os.path.join(config['projectDir'], config['expression_associations_dir'], 'expression_associations', '{normalization}', '{dataset}.txt'),  dataset=list(all_datasets), normalization=['tmm', 'tpm']),
        ###------------------------------------------------------------------
        ### Do PCA analysis on multiple processing pipelines, Arora et al
        expand(os.path.join(config['projectDir'], 'PCA_arora', '{project}_{count_type}.Rds'), project=['gtex', 'TCGA'], count_type=['TPM', 'RPKM']),
        expand(os.path.join(config['projectDir'], 'PCA_arora', '{project}_{count_type}.pca.txt'), project=['gtex', 'TCGA'], count_type=['TPM', 'RPKM']),
        ###------------------------------------------------------------------
        ## Look for variance explained of gene expression by transcreiptome diversity (PCA method)
        #
        expand(os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'PC_matrix.txt'), normalization=['tmm', 'tpm'], dataset=list(all_datasets)),
        expand(os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'PC_stats.txt'), normalization=['tmm', 'tpm'], dataset=list(all_datasets)),
        expand(os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'cors.txt'), normalization=['tmm', 'tpm'], dataset=list(all_datasets)),
        expand(os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'var_explained.txt'), normalization=['tmm', 'tpm'], dataset=list(all_datasets))
        ###------------------------------------------------------------------
        ## Do correlations with PEER in gtex
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', '{group}', 'cors.txt'), group=list(all_tissues_gtex)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', '{group}', 'stats.txt'), group=list(all_tissues_gtex)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', '{group}', 'PEER_transcriptome_diversity.txt'), group=list(all_tissues_gtex)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'var_explained_by_transcriptome_diversity', '{group}', 'lm_results.txt'), group=list(all_tissues_gtex)),
        ####------------------------------------------------------------------
        ### Do correlations with PEER in all datasets
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', '{dataset}', 'cors.txt'), dataset=list(all_datasets)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', '{dataset}', 'stats.txt'), dataset=list(all_datasets)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', '{dataset}', 'PEER_transcriptome_diversity.txt'), dataset=list(all_datasets)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'var_explained_by_transcriptome_diversity', '{dataset}', 'lm_results.txt'), dataset=list(all_datasets)),
        ###--------------------
        ## Do total variance explained of gene expression by PEER
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_by_PEER', '{dataset}', '{group}.txt'), zip, dataset=all_datasets_subset_gtex, group=all_datasets_subset_gtex_noPrefix),
        #expand(os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_by_PEER_nonGtex', '{dataset}.txt'), dataset=['lin']),
        ####------------------------------------------------------------------
        #### Do differential_gene_expression
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_genotype.Rds'), 
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_fulvestrant_wt.Rds'), 
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_tamoxifen_wt.Rds'), 
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_jq1_wt.Rds'),
        expand(os.path.join(config['projectDir'], config['differential_expression_dir'], 'sarantopoulou_{platform}', 'edgeR_results_ILB.Rds'),  platform=['truseq', 'v4', 'pico']),
        ###------------------------------------------------------------------
        ### Get how much var explained by PEER for each gene
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tpm', '{dataset}', '{group}.v8.covariates.txt.only_peer'), zip, group=list(all_datasets_subset_gtex_noPrefix), dataset=list(all_datasets_subset_gtex)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tpm', '{dataset}', '{group}', 'transcriptome_diversity.normalized_peer.txt'), zip, group=list(all_datasets_subset_gtex_noPrefix), dataset=list(all_datasets_subset_gtex)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tmm', '{dataset}', '{group}.v8.covariates.txt.only_peer'), zip, group=list(all_datasets_subset_gtex_noPrefix), dataset=list(all_datasets_subset_gtex)),
        expand(os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tmm', '{dataset}', '{group}', 'transcriptome_diversity.normalized_peer.txt'), zip, group=list(all_datasets_subset_gtex_noPrefix), dataset=list(all_datasets_subset_gtex))
        
        
        
        
#---------------------------------#
# Variance explained by all PEER covariates for each individual gene
# 

# Get peer normalized by trans diversity and other metadata
rule get_PEER_normalized:
    input:
        trans_diversity_as_cov = os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm_as_covariates', '{group}.txt'),
        cov = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
    params:
        out_prefix = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}')
    output:
        peer_out = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}', 'transcriptome_diversity.normalized_peer.txt')
    shell:
        '''
        cd R
        
        # Getting non-peer + trans diversity
        Rscript rbind.R {input.cov}.temp {input}
        grep -v InferredCov {input.cov}.temp > {input.cov}.non_PEER_transcriptome_diversity
        rm {input.cov}.temp
        
        # getting only peer
        if [ ! -f "{input.cov}.only_peer" ]
        then 
            grep -e ID -e InferredCov {input.cov} > {input.cov}.only_peer
        fi
        Rscript regress_out_from_peer.R {input.cov}.non_PEER_transcriptome_diversity {input.cov}.only_peer {params.out_prefix}/
        '''

rule PEER_normalized_gene_var_explained:
    input:
        exp_mat = lambda wc: os.path.join(config['projectDir'], sample_table_original.loc[wc.dataset, 'path']),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}', 'transcriptome_diversity.normalized_peer.txt')
    params: 
        in_dir = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}') 
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tpm', '{dataset}', '{group}', 'transcriptome_diversity.normalized_peer.txt')
    shell:
        '''
        cd R
        Rscript fit_covariates_to_expression.R {input.exp_mat} {params.in_dir} $(dirname {output})
        '''
        
rule PEER_gene_var_explained:
    input:
        exp_mat = lambda wc: os.path.join(config['projectDir'], sample_table_original.loc[wc.dataset, 'path']),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tpm', '{dataset}', '{group}.v8.covariates.txt.only_peer')
    shell:
        '''
        # Select only peer covariates
        if [ ! -f "{input.peer}.only_peer" ]
        then 
            grep -e ID -e InferredCov {input.peer} > {input.peer}.only_peer
        fi
        
        cd R
        
        Rscript fit_covariates_to_expression.R {input.exp_mat} {input.peer}.only_peer $(dirname {output})
        '''
        
rule PEER_normalized_gene_var_explained_tmm:
    input:
        exp_mat = os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tmm', '{dataset}.txt'),
        #exp_mat = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path'])
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}', 'transcriptome_diversity.normalized_peer.txt')
    params: 
        in_dir = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}') 
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tmm', '{dataset}', '{group}', 'transcriptome_diversity.normalized_peer.txt')
    shell:
        '''
        cd R
        Rscript fit_covariates_to_expression.R {input.exp_mat} {params.in_dir} $(dirname {output})
        '''
        
rule PEER_gene_var_explained_tmm:
    input:
        exp_mat = os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tmm', '{dataset}.txt'),
        #exp_mat = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path'])
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_per_gene', 'original_PEER', 'tmm', '{dataset}', '{group}.v8.covariates.txt.only_peer')
    shell:
        '''
        # Select only peer covariates
        if [ ! -f "{input.peer}.only_peer" ]
        then 
            grep -e ID -e InferredCov {input.peer} > {input.peer}.only_peer
        fi
        
        cd R
        
        Rscript fit_covariates_to_expression.R {input.exp_mat} {input.peer}.only_peer $(dirname {output})
        '''
        
        
#---------------------------------#
# Down-sample reads 
 
rule sampled_count_matrices:
    input:
        lambda wc: os.path.join(config['projectDir'], sample_table.loc[wc.dataset, 'path'])
    params:
        n_reads='{n_reads}'
    output:
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'],  'sampled_reads', '{n_reads}', '{dataset}.txt')
    shell:
        '''
        cd R;  Rscript get_sampled_reads.R {input} {params} {output}
        '''
        
        
#---------------------------------#
# Get expression estimates

rule get_tmm:
    input:
        lambda wc: os.path.join(config['projectDir'], sample_table.loc[wc.dataset, 'path'])
    output:
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tmm', '{dataset}.txt')
    shell:
        '''
        cd R; Rscript get_tmm.R {input} {output}
        '''
        
rule get_tpm:
    input: 
        lambda wc: os.path.join(config['projectDir'], sample_table.loc[wc.dataset, 'path']),
        lambda wc: os.path.join(config['projectDir'], sample_table.loc[wc.dataset, 'gtf_file'])
    output:
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tpm', '{dataset}.txt')
    shell:
        '''
        cd R; Rscript get_TPM.R {input} {output}
        '''
       
#---------------------------------#
# Get PEER covariates

rule get_PEER:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tmm', '{dataset}.txt')
    output:
        os.path.join(config['projectDir'], config['PEER_dir'], '{dataset}.txt')
    shell:
        '''
        cd R; Rscript PEER.R {input} {output}
        '''
       
#---------------------------------#
# get transcriptome diversity
 
rule get_transcriptome_diversity:
    input: 
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tpm', '{dataset}.txt')
    output: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', '{dataset}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general.R {input} {output}'
        
rule get_transcriptome_diversity_original_sources:
    input: 
        lambda wc: os.path.join(config['projectDir'], sample_table_original.loc[wc.dataset, 'path'])
    output: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_sources', '{dataset}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general.R {input} {output}'
        
# Downsampling reads 
rule get_transcriptome_diversity_from_sampled:
    input: 
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'sampled_reads', '{n_reads}', '{dataset}.txt')
    output: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'sampled_reads', '{n_reads}', '{dataset}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general.R {input} {output}'
        
# Getting genes that are expressed across all samples of a datasert
rule get_transcriptome_diversity_non_lowly_expressed:
    input: 
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], 'tpm', '{dataset}.txt')
    output: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'non_lowly_expressed_tpm', '{dataset}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general_remove_lowlyExpressed.R {input} {output}'
        
        
#---------------------------------#
# Do associations between transcreiptome diversity and expression on standardized matrices
        
rule expression_association_transcriptome_diversity:
    input: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', '{dataset}.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], '{normalization}', '{dataset}.txt')
    output:
        os.path.join(config['projectDir'], config['expression_associations_dir'], 'expression_associations', '{normalization}', '{dataset}.txt')
    shell:
        '''
        cd R; Rscript expression_vs_transcriptome_diversity_general.R {input} {output}
        '''
        
#---------------------------------#
# Do PCA on standardized matrices
        
rule get_cors_with_PCs:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_processed_dir'], '{normalization}', '{dataset}.txt'),
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', '{dataset}.txt')
    output:
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'PC_stats.txt'),
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'cors.txt')
    shell:
        'cd R/; Rscript get_expression_explained_transcriptome_diversity_general.R {input} {output}'
        
rule get_var_explained_by_transcriptome_diversity:
    input:
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'PC_stats.txt')
    output:
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', '{normalization}', '{dataset}', 'var_explained.txt')
    shell:
        'cd R/; Rscript get_expression_explained_transcriptome_diversity_general_PCA.R {input} {output}'
        
#---------------------------------#
# Analyzing multiple processing pipelines Arora
 
rule perform_PCA_arora:
    input:
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['expression_matrix_files'])
    params:
        root_dir=config['projectDir'],
        project='{project}',
        count_type='{count_type}'
    output:
        os.path.join(config['projectDir'], 'PCA_arora', '{project}_{count_type}.Rds'),
        os.path.join(config['projectDir'], 'PCA_arora', '{project}_{count_type}.pca.txt')
    shell:
        '''
        cd R; Rscript PCA_arora.R {input} {params} {output} 
        '''
    
        
#---------------------------------#
# eqtl analysis keele with as done in their paper
# 
 
rule get_transcriptome_diversity_as_covariates:
    input: 
        transcriptome_diversity=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'keele_{keele_tissue}.txt'),
        expression_bed=os.path.join(config['projectDir'], config['expression_matrices_dir'], 'keele_{keele_tissue}.bed.gz')
    output:
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'as_covariates', 'keele_{keele_tissue}.txt')
    shell:
        'cd R/; Rscript get_transcriptome_diversity_as_covariates_general.R {input} {output}'

rule perform_eqtl_as_paper_keele:
    input: 
        cov=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'as_covariates', 'keele_{keele_tissue}.txt'),
        expression_original=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', '{keele_tissue}_expression.csv.zip'),
        gene_info=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'refseq_mm9_tss.txt.zip'),
        genome_cache_dir=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'cc_genome_cache_full_l2_0.1/'),
        var_dir=os.path.join(config['projectDir'], config['expression_datasets_dir'], 'keele', 'data_all', 'isvdb_var_db/')
    params:
        covariate_mode='{covariate_mode}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'keele', '{covariate_mode}', '{keele_tissue}.txt.gz'),
        os.path.join(config['projectDir'], config['eqtl_dir'], 'keele', '{covariate_mode}', '{keele_tissue}.top_pvalues_per_gene.txt')
    shell:
        '''
        cd R
        Rscript eqtl_keele_method.R {input} {params} {output}
        '''

        
#---------------------------------#
# Differential gene expression analysis

rule dge_satantopoulou:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_{platform}.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'sarantopoulou_{platform}_cov.txt'),
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'sarantopoulou_{platform}.txt')
    params:
        os.path.join(config['projectDir'],  config['differential_expression_dir'], 'sarantopoulou_{platform}')
    output:
        os.path.join(config['projectDir'],  config['differential_expression_dir'], 'sarantopoulou_{platform}', 'edgeR_results_ILB.Rds')
    shell:
        '''
        cd R
        Rscript differential_expression_sarantopoulou.R {input} {params}
        '''
        
rule dge_nagarajan:
    input:
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan.txt'),
        os.path.join(config['projectDir'], config['expression_matrices_dir'], 'nagarajan_cov.txt'),
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', 'nagarajan.txt')
    params:
        os.path.join(config['projectDir'],  config['differential_expression_dir'], 'nagarajan')
    output:
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_genotype.Rds'),
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_fulvestrant_wt.Rds'),
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_tamoxifen_wt.Rds'),
        os.path.join(config['projectDir'], config['differential_expression_dir'], 'nagarajan', 'edgeR_results_jq1_wt.Rds')
    shell:
        'cd R; Rscript differential_expression_nagarajan.R {input} {params}'
        
#---------------------------------#
# eQTL search on GTEx data
# 
rule get_transcriptome_diversity_gtex_original_tpm:
    input: 
        os.path.join(config['projectDir'], 'expression_matrices', 'gtex_{group}.txt')
    output: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm', '{group}.txt')
    shell: 
        'cd R; Rscript get_transcriptome_diversity_general.R {input} {output}'
 
 
rule get_transcriptome_diversity_as_covariates_gtex:
    input: 
        trans_diversity=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm', '{group}.txt'),
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path'])
    output:
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm_as_covariates', '{group}.txt')
    shell:
        'cd R/; Rscript get_transcriptome_diversity_as_covariates.R {input} {output}'

rule merge_non_peer_covs_with_transcriptome_diversity:
    input: 
        trans_diversity=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm_as_covariates', '{group}.txt'),
        peer=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
    output:
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm_as_covariates_plus_non_peer', '{group}.txt')
    shell:
        '''
        Rscript R/rbind.R {output}.temp {input}
        grep -v InferredCov {output}.temp > {output}
        rm {output}.temp
        '''
        

# Performs eqtl search. Covariates: transcriptome_diversity
rule perform_eqtl_with_covariates:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        covariate_file=os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm_as_covariates', '{group}.txt'),
        vcf=config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_trans_as_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        '''
        
# Performs eqtl search. Covariates: transcriptome_diversity and non-peer covariates
rule perform_eqtl_with_covariates_plus:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        covariate_file = os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_gtex_tpm_as_covariates_plus_non_peer', '{group}.txt'),
        vcf=config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_trans_as_covariate_plus_non_peer', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        #../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --region chr11:5225000-5239000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.covariate_file} --permute 1000 10000 --out {output}
        '''
        
# Performs eqtl search. Covariates: no coovariates
rule perform_eqtl_with_no_covariates:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        vcf=config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_no_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --permute 1000 10000 --out {output}
        #../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --region chr11:5225000-5239000 --vcf {input.vcf} --bed {input.expression_bed} --permute 1000 10000 --out {output}
        '''
        
# Performs eqtl search. Covariates: original gtex covariates
rule perform_eqtl_with_full_covariates:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
        vcf = config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_full_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer} --permute 1000 10000 --out {output}
        '''
        
# Performs eqtl search. Covariates: PEER covariates
rule perform_eqtl_with_PEER_covariates:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
        vcf = config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_PEER_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        # Select only peer covariates
        if [ ! -f "{input.peer}.only_peer" ]
        then 
            grep -e ID -e InferredCov {input.peer} > {input.peer}.only_peer
        fi
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer}.only_peer --permute 1000 10000 --out {output}
        '''
        
# Performs eqtl search. Covariates: non-peer covariates
rule perform_eqtl_with_nonPEER_covariates:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt'),
        vcf = config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_nonPEER_covariate', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        # Select only peer covariates
        if [ ! -f "{input.peer}.only_non_peer" ]
        then 
            grep -v InferredCov {input.peer} > {input.peer}.only_non_peer
        fi
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer}.only_non_peer --permute 1000 10000 --out {output}
        '''
        
# Performs eqtl search. Covariates: PEER covariates for which transcriptome diversity has been regressed out 
rule perform_eqtl_with_PEER_covariates_regressed:
    input: 
        expression_bed = lambda wc: os.path.join(config['projectDir'], sample_table_gtex_bed.loc[wc.group, 'path']),
        peer = os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], 'normalized_PEER', '{group}', 'transcriptome_diversity.normalized_peer.txt'),
        vcf = config['vcfFileGTEXCommonV8']
    params:
        chunk = '{chunk}'
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_peer_covariate_trasn_div_regressed_out', '{group}.{chunk}.txt.gz')
    shell:
        '''
        module load gsl/1.16
        ../fastqtl/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --chunk {params.chunk} 2000 --vcf {input.vcf} --bed {input.expression_bed} --cov {input.peer} --permute 1000 10000 --out {output}
        '''

rule join_gtex_eqtl_results:
    input:
        #expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_trans_as_covariate', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_PEER_covariate', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_nonPEER_covariate', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_full_covariate', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_no_covariate', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001)),
        #expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_trans_as_covariate_plus_non_peer', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001)),
        expand(os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_peer_covariate_trasn_div_regressed_out', '{group}.{chunk}.txt.gz'), group=list(tissues_eqtl_gtex), chunk=range(1,2001))
    params:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_trans_as_covariate'),
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_peer_covariate'),
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_nonpeer_covariate'),
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_full_covariate'),
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_no_covariate'),
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'eqtls_trans_as_covariate_plus_non_peer')
    output:
        os.path.join(config['projectDir'], config['eqtl_dir'], 'gtex', 'all_results', 'gtex_eqtl_results.txt')
    shell:
        'cd R; #Rscript eqtl_join_gtex_results.r {params} {output}' 
        
        
#---------------------------------#
# PCA (var explained by transcriptome_diversity)
 
        
rule get_cors_with_PCs_original_sources:
    input: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_sources', '{dataset}.txt'),
        lambda wc: os.path.join(config['projectDir'], sample_table_original.loc[wc.dataset, 'path'])
    output:
        os.path.join(config['projectDir'], 'var_explained_by_transcriptome_diversity', 'pca_method', 'original_sources', '{dataset}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'], 'var_explained_by_transcriptome_diversity', 'pca_method', 'original_sources', '{dataset}', 'PC_stats.txt'),
        os.path.join(config['projectDir'], 'var_explained_by_transcriptome_diversity', 'pca_method', 'original_sources', '{dataset}', 'cors.txt')
    shell:
        'cd R/; Rscript expression_explained_by_transcriptome_diversity_general.R {input} {output}'
        

## GTEx PEER 
        
rule get_cors_with_PEER:
    input: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_sources', 'gtex_{group}.txt'),
        peer=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
    params:
        pvalue='0.05'
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', '{group}', 'cors.txt'),
        os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', '{group}', 'stats.txt'),
        os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'cors', '{group}', 'PEER_transcriptome_diversity.txt')
    shell:
        'cd R/; Rscript PEER_vs_transcriptome_diversity.R {input} {params} {output}'
        
rule get_var_explained_transcriptome_diversity_by_PEER:
    input: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'original_sources', 'gtex_{group}.txt'),
        peer=os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'gtex_original', 'var_explained_by_transcriptome_diversity', '{group}', 'lm_results.txt'),
    shell:
        'cd R/; Rscript PEER_regression_on_transcriptome_diversity.R {input} {output}'
        
## all PEER 
        
rule get_cors_with_PEER_all:
    input: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', '{dataset}.txt'),
        os.path.join(config['projectDir'], config['PEER_dir'], '{dataset}.txt')
    params:
        pvalue='0.05'
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', '{dataset}', 'cors.txt'),
        os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', '{dataset}', 'stats.txt'),
        os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'cors', '{dataset}', 'PEER_transcriptome_diversity.txt')
    shell:
        'cd R/; Rscript PEER_vs_transcriptome_diversity.R {input} {params} {output}'
        
rule get_var_explained_transcriptome_diversity_by_PEER_all:
    input: 
        os.path.join(config['projectDir'], config['transcriptome_diversity_dir'], 'tpm', '{dataset}.txt'),
        os.path.join(config['projectDir'], config['PEER_dir'], '{dataset}.txt')
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'tpm', 'var_explained_by_transcriptome_diversity', '{dataset}', 'lm_results.txt'),
    shell:
        'cd R/; Rscript PEER_regression_on_transcriptome_diversity.R {input} {output}'
        
        
rule get_var_explained_by_PEER_gtex_original:
    input: 
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', 'tpm', '{dataset}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', 'tpm', '{dataset}', 'PC_stats.txt'),
        os.path.join(config['projectDir'], config['auxiliary_Files']['dir'], config['auxiliary_Files']['gtexPeerV8Dir'], '{group}.v8.covariates.txt')
        #os.path.join(config['projectDir'], config['PEER_dir'], '{dataset}.txt')
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_by_PEER', '{dataset}', '{group}.txt'),
    shell:
        'cd R/; Rscript get_expression_explained_covariates_gtex.R {input} {output}'
        
rule get_var_explained_by_PEER_all:
    input: 
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', 'tpm', '{dataset}', 'PC_matrix.txt'),
        os.path.join(config['projectDir'],  'var_explained_by_transcriptome_diversity', 'pca_method', 'from_raw', 'tpm', '{dataset}', 'PC_stats.txt'),
        os.path.join(config['projectDir'], config['PEER_dir'], '{dataset}.txt')
    output:
        os.path.join(config['projectDir'], 'PEER_analyses', 'var_explained_by_PEER_nonGtex', '{dataset}.txt'),
    shell:
        'cd R/; Rscript get_expression_explained_covariates_gtex.R {input} {output}'
        

