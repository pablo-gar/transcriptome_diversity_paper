# Calculates shannon diversity for all samples in gtex
# Usage
# Rscript get_transcriptome_diversity.R out_table.txt

source("../../R/gtex.R", chdir = T)
source("../../R/misc.R", chdir = T)

main <- function(cmdArgs = commandArgs(T)) {
    
    exp_mat_file <- cmdArgs[1]
    samples <- parseArg(cmdArgs[2], sep=',')
    do_splicing <- as.logical(cmdArgs[3])
    out_file <- cmdArgs[4]
    
    #do_splicing <- T
    #exp_mat_file <- GTEX_CON$exonJunctionCountsV8
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/transcriptome_diversity/exonJunctions/Whole_Blood.txt'
    #tissue <- 'Whole_Blood'
    #do_splicing <- F
    #exp_mat_file <- GTEX_CON$readsGenesAllTissues
    #samples <- tissueToGtexLong('Whole_Blood')[,2]
    #out_file <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/transcriptome_diversity/expression/Whole_Blood.txt'
    #aaa <- read.table('/scratch/users/paedugar/cnv_gtex_project/simple_cnv_caller/chromosme_arm_counts_from_matrix/Whole_Blood.txt', header=T, check.names=F)
    #samples <- colnames(aaa)[grep('GTEX', colnames(aaa))]
    
    samples <- unique(samples)
    if(do_splicing) {
        exp_mat <- readExonJunctionCounts(samples=samples, expressionFile=exp_mat_file)
    } else {
        exp_mat <- readAllGtexExpression(samples=samples, expressionFile=exp_mat_file)
    }
    
    out_table <- get_shannon(exp_mat)
    write.table(out_table, out_file, sep = "\t", quote = F, row.names = F, col.names = T)
}


# calculates the shanon entropy of a vector
shannonEntropy <- function(x, percentage = T) {
    
    x <- x[ x != 0 & !is.na(x) ]
    
    shan <- -sum( (x/sum(x)) * log2(x/sum(x)) )
    
    if(percentage)
        shan <- shan / log2(length(x))
    
    return(shan)
}

# Calcautes the shanon entropy across columns in a matrix in an efficent way
# IT DOES NOT ELIMINATE 0s
shannonEntropyMatrix <- function(x) {
    
    small <- 0.000001
    
    # Doing frequencies by column
    x <- x / matrix(colSums(x), nrow = nrow(x), ncol = ncol(x), byrow=T)
    
    shan <- -colSums((x + small) * log2(x + small))
    
    return(shan)
}

get_shannon <- function(geneExp, ignoreZero = T) {

    
    if(ncol(geneExp) < 2 | nrow(geneExp) < 2)
        return(NULL)
    
    rownames(geneExp)  <- geneExp[,1]
    geneExp <- geneExp[,-1]
    
    if(ignoreZero) {
        shan <- apply(as.matrix(geneExp), 2, shannonEntropy)
    } else {
        shan <- shannonEntropyMatrix(as.matrix(geneExp))
    }
    

    shan <- data.frame(gtexId = names(shan), transcriptome_diversity = shan)
    
    return(shan)
    
}

main()
