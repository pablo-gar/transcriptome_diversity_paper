library('MatrixEQTL')

if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

#source("../../R/GWAS_methods.R", chdir = T)
source("../../R/transcriptome_diversity_tools.R", chdir = T)

main <- function(cmdArgs=commandArgs(T)) {

    transcriptome_diversity_file <- cmdArgs[1]
    expression_mat_file <- cmdArgs[2]
    out_file <- cmdArgs[3]

    #transcriptome_diversity_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity/keele_kidney.txt'
    #expression_mat_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/keele_kidney.txt'
    #transcriptome_diversity_file <- '/scratch/users/paedugar/transcriptome_diversity/transcriptome_diversity/tpm/sarantopoulou_pico.txt'
    #expression_mat_file <- '/scratch/users/paedugar/transcriptome_diversity/expression_matrices_processed/tpm/sarantopoulou_pico.txt'
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/expression_associations/keele_kidney.txt'


    #---------------------------
    # MAIN

    # Reading mutation data and metadata
    transcriptome_diversity <- read.table(transcriptome_diversity_file, header=T, sep="\t", stringsAsFactors=F, check.names=F)
    transcriptome_diversity <- matrix(transcriptome_diversity[,2], nrow=1, dimnames=list('transcriptome_diversity', transcriptome_diversity[,1]))

    # Reading expression mat
    expression_mat <- read_expression(expression_mat_file)
    rownames(expression_mat) <- expression_mat[,1]
    expression_mat <- as.matrix(expression_mat[,colnames(transcriptome_diversity)])
    
    # Filter genes that are lowly expressed in 20% of samples of more 
    #expression_mat <- expression_mat[rowSums(expression_mat>1) > (0.2 *(ncol(expression_mat))), ,drop=F]
    expression_mat <- filter_rarely_expressed_genes(expression_mat, percentage=0.2, gene_names_as_rownames=T, count_cutoff=1)
    
    transcriptome_diversity <- rankitNormalize(transcriptome_diversity)
    expression_mat <- rankitNormalize(expression_mat)


    x <- lmMatFunction(transcriptome_diversity, expression_mat)
    dfFull <- x$param$dfFull
    
    x <- x$all$eqtls
    
    # Getting r2
    tstat <- x$statistic;
    r <- tstat / sqrt( dfFull + tstat^2 )
    x$r2 <- r^2
    
    # Ordering
    x$pvalue_bonf <- p.adjust(x$pvalue)
    x <- x[order(x$pvalue),]
        
    write.table(x, out_file, sep = "\t", col.names = T, row.names = F, quote = F)
}

lmMatFunction <- function(x, y, useModel = modelLINEAR, cvrt = NULL, pvalCutoff = 1, errorCovariance = numeric(), outFile = tempfile(), min.pv.by.genesnp = F, noFDRsaveMemory = F){
    
    xMat <- SlicedData$new()
    yMat <- SlicedData$new()
    cvrtMat <- SlicedData$new()
    
    xMat$CreateFromMatrix(x)
    yMat$CreateFromMatrix(y)
    
    if(!is.null(cvrt)) 
        cvrtMat$CreateFromMatrix(cvrt)
    
    
    results <- Matrix_eQTL_engine(snps = xMat, gene = yMat, cvrt = cvrtMat, output_file_name = outFile,
                      pvOutputThreshold = pvalCutoff, useModel = useModel, min.pv.by.genesnp = min.pv.by.genesnp, noFDRsaveMemory = noFDRsaveMemory)
    
    return(results)
}

main()
