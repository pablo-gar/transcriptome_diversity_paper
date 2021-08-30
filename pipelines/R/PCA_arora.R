source('../../R/ggthemes.R', chdir=T)
source('../../R/transcriptome_diversity_tools.R', chdir=T)
library('readr')
library('dplyr')
library('tidyr')
library('ggplot2')
library('stringr')


main <- function(cmdArgs = commandArgs(T)) {
    
    file_descriptions_file <- cmdArgs[1]
    root_dir <- cmdArgs[2]
    project <- cmdArgs[3]
    count_type <- cmdArgs[4]
    out_file <- cmdArgs[5]
    out_file_pca <- cmdArgs[6]
    
    #file_descriptions_file <- '/scratch/users/paedugar/transcriptome_diversity/auxiliary_files/expression_datasets.txt'
    #root_dir <- '/scratch/users/paedugar/transcriptome_diversity'
    #project <- 'gtex'
    #count_type <- 'TPM'
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/PCA_arora/gtex_TPM.Rds'
    #out_file_pca <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/PCA_arora/gtex_TPM_pca.txt'
    
    # Read data
    cat("Reading data \n")
    descriptions <- read_tsv(file_descriptions_file)
    descriptions <- descriptions[descriptions$group=='arora',]
    descriptions$path <- file.path(root_dir, descriptions$path)

    current <- descriptions[descriptions$subgroup==project,]
    current <- current[current$count_type==count_type,]

    all_exp <- read_expression_files(current)

    # Filter lowly expressed genes
    #all_exp$exp_mat <- filter_rarely_expressed_genes(all_exp$exp_mat, percentage=0.3, gene_names_as_rownames=T)
    #all_exp$exp_mat_controlled_transcriptome_diversity <- all_exp$exp_mat_controlled_transcriptome_diversity[rownames(all_exp$exp_mat),]

    # Do pca
    cat("Doing PCA \n")
    all_exp$pca <- prcomp(t(log2(as.matrix(all_exp$exp_mat))), center=F, scale.=F)
    all_exp$pca_controlled_transcriptome_diversity <- prcomp(t(log2(as.matrix(all_exp$exp_mat_controlled_transcriptome_diversity))), center=F, scale.=F)
    
    # join pca results with sample info
    x <- join_PCA_info(all_exp$pca$x, all_exp$info)
    y <- join_PCA_info(all_exp$pca_controlled_transcriptome_diversity$x, all_exp$info)
    
    x$pca_type <- 'original'
    y$pca_type <- 'transcriptome_diversity_controlled'
    pca_results <- bind_rows(x,y)

    # Save results
    saveRDS(all_exp, out_file)
    write_tsv(pca_results, out_file_pca)

}

read_expression_files <- function(x) {
    
    x$n_samples <- 0

    # initialize matrices with fir expression profile
    exp_mat <- read.table(x$path[1], sep='\t', stringsAsFactors=F, header=T, row.names=1)
    x$n_samples[1] <- ncol(exp_mat)
    rows <- rownames(exp_mat)
    colnames(exp_mat) <- paste0(x$id[1], ',', colnames(exp_mat))

    # Normlaized by trans diversity
    exp_mat_controlled_transcriptome_diversity <- get_normalized(exp_mat)
    for(i in 2:nrow(x)) {
    #for(i in 2:4) {
        
        cat(x$path[i], "\n")
        c_mat <- read.table(x$path[i], sep='\t', stringsAsFactors=F, header=T, row.names=1) 
        
        if(nrow(c_mat) != nrow(exp_mat))
            next
        
        colnames(c_mat) <- paste0(x$id[i], ',', colnames(c_mat))
        
        # get normalized by trans diversity
        c_mat_normalized <- get_normalized(c_mat)
        
        # Appending
        c_mat <- c_mat[rows,]
        c_mat_normalized <- c_mat_normalized[rows,]
        
        exp_mat <- cbind(exp_mat, c_mat)
        exp_mat_controlled_transcriptome_diversity <- cbind(exp_mat_controlled_transcriptome_diversity, c_mat_normalized)
        
        # Counting samples in current expression
        x$n_samples[i] <- ncol(c_mat)
        
    }
    
    x <- x[x$n_samples>0,]
    
    exp_mat <- exp_mat + 0.001
    exp_mat_controlled_transcriptome_diversity <- exp_mat_controlled_transcriptome_diversity + 0.001
    
    exp_mat[is.na(exp_mat)] <- min(as.matrix(exp_mat), na.rm=T)
    exp_mat_controlled_transcriptome_diversity[is.na(exp_mat_controlled_transcriptome_diversity)] <- min(as.matrix(exp_mat_controlled_transcriptome_diversity), na.rm=T)
    return(list(exp_mat=exp_mat, exp_mat_controlled_transcriptome_diversity=exp_mat_controlled_transcriptome_diversity, info=x))
    
}

get_normalized <- function(x) {
    
    # Filter lowly expressed genes
    c_mat_filtered <- filter_rarely_expressed_genes(x, gene_names_as_rownames=T)
    
    # Get transcriptome diversity
    trans <- get_transcriptome_diversity(c_mat_filtered, gene_names_as_rownames=T, do_rankitNormalization=F, percentage=F)
    
    # Get normalized matrix
    c_mat_normalized <- x / matrix(rep(trans[colnames(x), 'transcriptome_diversity'], each=nrow(x)), byrow=F, ncol=ncol(x))
   
    return(c_mat_normalized)
}

join_PCA_info <- function(x, info) {
    
    
    x <- as.data.frame(x[,1:6])
    x$id<- str_split_fixed(rownames(x), ',', 2)[,1]
    
    x <- left_join(x, info)
    
    return(x)
    
}

main()
