main <- function(cmdArgs=commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    anno_file <- cmdArgs[2]
    out_dir <- cmdArgs[3]
    
    #exp_mat_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_datasets/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
    #anno_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_datasets/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
    #out_dir <- '/scratch/users/paedugar/transcriptome_diversity/expression_matrices/'
    
    cat("Read expression matrix\n")
    exp_mat <- read.table(gzfile(exp_mat_file), skip=2, header=T, sep='\t', stringsAsFactors=F, check.names=F)
    colnames(exp_mat)[1] <- 'gene_id'
    
    cat("Read annotation file\n")
    anno <- read.delim(anno_file, sep='\t', stringsAsFactors=F, header=T)[,c(1,7)]
    anno[,2] <- code_friendly_tissues(anno[,2])
    
    for(tissue in unique(anno[,2])) {
        cat(tissue, "\n")
        samples <- anno[anno[,2] == tissue, 1]
        samples <- samples[samples %in% colnames(exp_mat)]
        current <- exp_mat[,c('gene_id', samples), drop=F]
        
        if(ncol(current) > 1) {
            write.table(current, file.path(out_dir, paste0('gtex_', tissue, '_counts.txt')), sep='\t', col.names=T, row.names=F, quote=F)
            cat('Done\n')
        }
    }
    
    
}

code_friendly_tissues <- function(x) {
    
    x <- gsub('-', '', x)
    x <- gsub('\\(', '', x)
    x <- gsub('\\)', '', x)
    x <- gsub('\\s+', '_', x)
    
    return(x)
}

main()
