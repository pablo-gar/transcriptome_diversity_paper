library('SummarizedExperiment')

main <- function(cmdArgs=commandArgs(T)) {

    base_dir <- cmdArgs[1]
    out_dir_base <- cmdArgs[2]
    out_dir_final <- cmdArgs[3]
    out_file_table <- cmdArgs[4]
    
    #base_dir <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/arora/OriginalTCGAGTExData/SE_objects/'
    #out_dir <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/'
    #out_file_table <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/expression_datasets_arora_2.txt'
    
    datasets <- get_dataset_info(base_dir)
    gtex_tissues <- read_gtex_tissues()
    
    all_datasets_table <- list()
    
    for(i in 1:nrow(datasets)) {
        
        current_name <- as.character(datasets$name[i])
        current_file <- datasets$filename[i]
        is_log2 <- datasets$is_log2[i]
        
        cat("Working with: ", current_name, "\n")
        
        if(!file.exists(current_file))
            stop('\n',current_file, 'does not exist. Are you sure all Arora et al. data has been dowloaded?')
        
        current_data <- get(load(current_file))
        
        gene_mat <- as.matrix(assay(current_data))
        
        if(is_log2) {
            gene_mat <- 2^gene_mat
            gene_mat[is.infinite(gene_mat)] <- 0
            gene_mat <- as.data.frame(gene_mat)
        }
        
        sample_attrs <- as.data.frame(colData(current_data))
        gene_attrs <- as.data.frame(rowData(current_data))
        
        if(grepl('gtex_v6', current_name)) {
            group_index <- 7 
            gene_index <- 5
        } else if (grepl('gtex_xena', current_name)){
            sample_attrs <- gtex_tissues
            sample_attrs <- sample_attrs[colnames(gene_mat)[-1],]
            group_index <- 2
        } else if (grepl('TCGA_recount2', current_name)){
            group_index <- 3
        } else if (grepl('TCGA_piccolo', current_name)){
            gene_index <- 5
            group_index <- 2
        } else if (grepl('TCGA_gdc', current_name)){
            gene_index <- 3
            group_index <- 2
        } else if (grepl('gtex_recount2', current_name)){
            gene_index <- 1
            group_index <- 3
        } else {
            group_index <- 2
            gene_index <- 1
        }
        
        if(grepl('xena', current_name)) {
            gene_index <- 5
        }
        
        
        # Appeding gene names
        gene_id <-data.frame(name = fix_gene_names(gene_attrs[,gene_index]))
        gene_mat <- cbind(gene_id, gene_mat)
        groups <- unique(sample_attrs[,group_index])
        
        for(group in groups) {
            
            
            # Merging colon samples
            if(grepl('olon', group)) {
                group <- groups[grep('olon', groups)]
            } else if(grepl('ervix', group)) {
                group <- groups[grep('ervix', groups)]
            } else if(grepl('alivar', group)) {
                group <- groups[grep('alivar', groups)]
            }
            
            cat("    Group: ", group, "\n")
            
            current_samples <- rownames(sample_attrs[sample_attrs[,group_index] %in% group, ,drop=F])
            current_gene_mat <- gene_mat[,c('name',current_samples)]
            
            group <- group[1]
            
            group <- code_friendly_tissues(group)
            group <- standardized_tissues(group)
            
            if(grepl('TCGA', current_name)) {
                group <- toupper(group)
            }
            
            id_file <- paste0('arora', '_', current_name, '_', group)
            out_file <- file.path(out_dir_base, out_dir_final, paste0(id_file, '.txt'))
            relative_path <- file.path(out_dir_final, paste0(id_file, '.txt'))
            
            info <- unlist(strsplit(current_name, '_'))
            all_datasets_table <- c(all_datasets_table, list(data.frame(id=id_file, name=id_file, path=relative_path, group='arora', subgroup=info[1], method=info[2], count_type=info[3], tissue=group)))
            
            write.table(current_gene_mat, out_file, sep='\t', quote=F, row.names=F)
        }
        
    }
    all_datasets_table <- do.call(rbind, all_datasets_table)
    all_datasets_table <- unique(all_datasets_table)
    
    write.table(all_datasets_table, out_file_table, sep='\t', col.names=T, row.names=F, quote=F)
}


fix_gene_names <- function(x) {
    
    gsub('\\..+$', '', x)
    
}

get_dataset_info <- function(base_dir) {
    
    
    results <- rbind(data.frame(name='gtex_mskccBatch_TPM', is_log2=T, filename='gtex_mskcc_batch_log2_TPM.RData'),
                          data.frame(name='gtex_mskcc_TPM', is_log2=T, filename='gtex_mskcc_norm_log2_TPM.RData'),
                          data.frame(name='gtex_mskccBatch_RPKM', is_log2=F, filename='RPKM_gtex_mskcc_batch.RData'),
                          data.frame(name='gtex_mskcc_RPKM', is_log2=F, filename='RPKM_gtex_mskcc_norm.RData'),
                          
                          data.frame(name='gtex_recount2_TPM', is_log2=T, filename='gtex_recount2_log2_TPM.RData'),
                          data.frame(name='gtex_recount2_RPKM', is_log2=F, filename='RPKM_gtex_recount2.RData'),
                          
                          data.frame(name='gtex_v6_TPM', is_log2=T, filename='gtex_v6_log2_TPM.RData'),
                          data.frame(name='gtex_v6_RPKM', is_log2=F, filename='RPKM_gtex_original.RData'),
                          
                          data.frame(name='gtex_xena_TPM', is_log2=T, filename='gtex_xena_log2_TPM.RData'),
                          data.frame(name='gtex_xena_RPKM', is_log2=T, filename='RPKM_gtex_xena.RData'),
                          
                          #data.frame(name='gtex_xena_TPM', is_log2=T, filename='gtex_xena_TPM.RData'),
                          #data.frame(name='gtex_xena_RPKM', is_log2=F, filename='RPKM_xena.RData'),
                          
                          
                          data.frame(name='TCGA_mskccBatch_RPKM', is_log2=F, filename='RPKM_mskcc_batch.RData'), 
                          data.frame(name='TCGA_mskcc_RPKM', is_log2=F, filename='RPKM_mskcc_norm.RData'), 
                          data.frame(name='TCGA_mskccBatch_TPM', is_log2=T, filename='tcga_mskcc_batch_log2_TPM.RData'),
                          data.frame(name='TCGA_mskcc_TPM', is_log2=T, filename='tcga_mskcc_norm_log2_TPM.RData'),
                          
                          data.frame(name='TCGA_recount2_RPKM', is_log2=F, filename='RPKM_rse_tcga_recount2.RData'), 
                          data.frame(name='TCGA_recount2_TPM', is_log2=T, filename='tcga_recount2_log2_TPM.RData'),
                          
                          
                          data.frame(name='TCGA_gdc_TPM', is_log2=T, filename='tcga_gdc_log2_TPM.RData'),
                          data.frame(name='TCGA_gdc_RPKM', is_log2=F, filename='RPKM_gdc.RData'),
                          
                          data.frame(name='TCGA_piccolo_RPKM', is_log2=F, filename='RPKM_TCGA_gse62944_tumor.RData'), 
                          data.frame(name='TCGA_piccolo_TPM', is_log2=T, filename='tcga_piccolo_log2_TPM.RData'),
                          
                          data.frame(name='TCGA_xena_RPKM', is_log2=T, filename='RPKM_xena.RData'), 
                          data.frame(name='TCGA_xena_TPM', is_log2=T, filename='tcga_xena_log2_TPM.RData')
                     )
    results$filename <- file.path(base_dir, results$filename)
    
    return(results)
}

code_friendly_tissues <- function(x) {
    
    x <- gsub('-', '_', x)
    x <- gsub('\\(', '', x)
    x <- gsub('\\)', '', x)
    x <- gsub('\\s+', '_', x)
    x <- gsub('_+', '_', x)
    x <- gsub('&', '', x)
    
    x <- tolower(x)
    
    return(x)
}

standardized_tissues <- function(x) {
    
    if(x=='esophagus_gastroesophageal_junction')
        return('esophagus_gas')
    
    if(x=='esophagus_mucosa')
        return('esophagus_muc')
    
    if(x=='esophagus_muscularis')
        return('esophagus_mus')
    
    if(x=='kidney_cortex')
        return('kidney')
    
    if(grepl('olon', x))
        return('colon')
    
    if(grepl('alivar', x))
        return('salivary')
    
    if(grepl('ervix', x))
        return('cervix')
    
    if(x=='bladder_urothelial_carcinoma')
        return('BLCA')
    
    if(x=='breast_invasive_carcinoma')
        return('BRCA')
    
    if(x=='colon_adenocarcinoma')
        return('COAD')
    
    if(x=='head__neck_squamous_cell_carcinoma')
        return('HNSC')
    
    if(x=='kidney_chromophobe')
        return('KICH')
    
    if(x=='kidney_clear_cell_carcinoma')
        return('KIRC')
    
    if(x=='kidney_papillary_cell_carcinoma')
        return('KIRP')
    
    if(x=='liver_hepatocellular_carcinoma')
        return('LIHC')
    
    if(x=='lung_adenocarcinoma')
        return('LUAD')
    
    if(x=='lung_squamous_cell_carcinoma')
        return('LUSC')
    
    if(x=='prostate_adenocarcinoma')
        return('PRAD')
    
    if(x=='rectum_adenocarcinoma')
        return('READ')
    
    if(x=='stomach_adenocarcinoma')
        return('STAD')
    
    if(x=='thyroid_carcinoma')
        return('THCA')
    
    if(x=='uterine_corpus_endometrioid_carcinoma')
        return('UCEC')
    
    if(grepl('reast', x))
        return('breast')
    
    
    if(grepl('tcga', x)) {
        
        x <- unlist(strsplit(x, '_'))
        
        if(length(x) == 2) {
            return(x[2])
        } else {
            return(x[1])
        }
    }
    
        
    
    return(x)
}

read_gtex_tissues <- function() {
    
    x <- read.table('./process_arora_data_gtex_table.txt', header=F, sep='\t', stringsAsFactors=F)
    rownames(x) <- x[,1]
    return(x)
    
}

main()
