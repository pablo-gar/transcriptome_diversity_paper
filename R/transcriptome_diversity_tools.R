#' Read expression matrix
#' @param x path to file
#' @return data.frame 

read_expression <- function(x) {
    
    result <- read.table(x, sep='\t', header=T, check.names=F)
    
    colnames(result)[1] <- 'gene_id'
    return(result)
    
}

#' Reads a covariate matrix that's ready for fastQTL and transposes for better use in R
#' @param x path to file
#' @return data.frame 

read_covariates_t <- function(x) {
    
    covariates <- t(read_tsv(x))
    colnames(covariates) <- covariates[1,]
    covariates <- as.data.frame(covariates[-1,], stringsAsFactors=F)
    covariates$sample_id <- rownames(covariates)
    
    return(covariates)
    
}

#' Filters out genes that are expressed under a percentage of samples
#' @param x expression matrix (first col is gene id, others are samples)
#' @return data.frame with genes ketp

filter_rarely_expressed_genes <- function(x, percentage=0.2, gene_names_as_rownames=F, count_cutoff=1) {
    if(gene_names_as_rownames) {
        return(x[c(rowSums(as.matrix(x) > count_cutoff) >= (ncol(x) * percentage)),])
    } else {
        return(x[c(rowSums(as.matrix(x[,-1]) > count_cutoff) >= (ncol(x[,-1]) * percentage)),])
    }
}


#' Perform PCA on expression matrix
#' @param x expression matrix (first col is gene id, others are samples)
#' @param ... to be passed to prcomp
#' @return prcomp object

pca_expression <- function(x, ...) {
    rownames(x) <- x[,1]
    prcomp(t(x[,-1]), ...)
}

#' Rankit normalize expression across genes
#' @param x expression matrix (first col is gene id, others are samples)
#' @return data.frame same dimension as x
rankit_normalize_expression <- function(x, gene_names_as_rownames=F) {
    
    if(gene_names_as_rownames) {
        return(rankitNormalize(as.matrix(x)))
    } else {
        return(cbind(x[,1, drop=F], as.data.frame(rankitNormalize(as.matrix(x[,-1])))))
    }
}

rankitNormalize <- function(x, IND = 1) {

    # Normalizes rows (IND = 1) or columns (IND = 2)
    # to a quantile standard normalization (i.e. rankit)

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(x))
    
    rowNames <- rownames(x)
    colNames <- colnames(x)

    x <- apply(x, IND, rankitNormalize_vector)
    if(IND == 1)
        x <- t(x)
    
    rownames(x) <- rowNames
    colnames(x) <- colNames

    return(x)

}

rankitNormalize_vector <- function(x) {

    stopifnot(is.numeric(x))
    noNa <- !is.na(x)
    x[noNa] <- qnorm((rank(x[noNa]) - 0.5) / sum(noNa))
    return(x)

}
#' Calculates transcriptome diverisity per sample
get_transcriptome_diversity <- function(geneExp, percentage=T, gene_names_as_rownames=F, do_rankitNormalization=F, count_cutoff=0, expressed_in_all=F) {
    
    if(ncol(geneExp) < 2 | nrow(geneExp) < 2)
        return(NULL)
    
    if(!gene_names_as_rownames) {
        rownames(geneExp)  <- geneExp[,1]
        geneExp <- geneExp[,-1]
    }
    
    if(expressed_in_all)
        geneExp <- filter_rarely_expressed_genes(geneExp, percentage=1, count_cutoff=count_cutoff)
    
    shan <- apply(as.matrix(geneExp), 2, shannonEntropy, percentage=percentage, count_cutoff=count_cutoff)
    
    if(do_rankitNormalization)
        shan <- rankitNormalize_vector(shan)

    shan <- data.frame(sample_id = names(shan), transcriptome_diversity = shan)
    
    return(shan)
    
}

#' Calculates the shanon entropy of a vector
#'
#' @param count_cutoff eliminate elements having counts less or equal than this value
#' @return data.frame 
shannonEntropy <- function(x, percentage = T, count_cutoff = 0) {
    
    #x <- x[ x != 0 & !is.na(x) ]
    x <- x[!is.na(x)]
    x <- x[ x > count_cutoff ]
    
    shan <- -sum( (x/sum(x)) * log2(x/sum(x)) )
    
    if(percentage)
        shan <- shan / log2(length(x))
    
    return(shan)
}

#' Normalize transcriptome diversity in expression matrix

normalize_transcriptome_diversity <- function(geneExp, percentage=T) {
    
    transcriptome_diversity <- get_transcriptome_diversity(geneExp, percentage=percentage)
    transcriptome_diversity <- transcriptome_diversity[colnames(geneExp)[-1], 2]
    
    transcriptome_diversity <- rankitNormalize_vector(transcriptome_diversity)
    geneExp <- rankit_normalize_expression(geneExp)
    
    normalized <- t(apply(geneExp[,-1, drop=F], 1, function(x) resid(lm(x ~ transcriptome_diversity))))
    
    geneExp[,-1] <- normalized[,colnames(geneExp)[-1]]
    
    return(geneExp)
    
}

#' Prabilistic resizing of count vector
#' 
#' @param x orignal vector
#' @param n total final counts
resize_count_vector <- function(x, n, method='montecarlo') {
    
    if(method=='montecarlo') {
        
        probs <- x/sum(x)
        y <- sample(1:length(x), size=n, replace=T, prob=probs)
        y <- table(y)
        
        x <- setNames(rep(0, length(x)), names(x))
        x[as.numeric(names(y))] <- as.numeric(y)
        
        
    } else if(method=='arithmetic') { 
        names_x <- names(x)
        x <- arithmetic_resize(x=x, n=n)
        names(x) <- names_x
        
    } else {
        stop ('method has to be one of the following: montecarlo, arithmetic')
    }
    return(x)
}

arithmetic_resize <- function(x, n) {
    
    if (sum(x) <= n | sum(x) <=  0)
        return(x)
    
    total <- sum(x)
    down <- round((total-n) / sum(x>0))
    down <- ifelse(down>min(x[x>0]), min(x[x>0]), down)
    #cat('down: ', down, "total: ", total, "n: ", n, "\n")
    x[x>0] <- x[x>0] - down
    
    arithmetic_resize(x, n)
    
}

#' Prabilistic resizing of count matrix
#' 
#' @param x orignal matrix, where the columns will be idepently resized
#' @param n total final counts, either a vector to be applied to each column or a single number applied to all columns

resize_count_matrix <- function(x, n, method='montecarlo') {
    
    if(!(is.matrix(x) & is.numeric(x)))
        stop("x has to be a numeric matrix")
    
    if(!(length(n)!=1 | length(n)!=ncol(x)))
        stop("n has to be either a sinlge value or a vector with size equal to the number of columns in x ")
    
    for(i in 1:ncol(x)) {
        
        if(length(n) == 1) {
            current_n <- n 
        } else {
            current_n <- n[i]
        }
        
        x[,i] <- resize_count_vector(x[,i], current_n, method=method)
    }
    
    return(x)
    
}

#' Prabilistic resizing of count matrix
#' 
#' @param x orignal matrix, where the columns will be idepently resized
#' @param n total final counts, either a vector to be applied to each column or a single number applied to all columns

resize_expression <- function(x, n, method='montecarlo', gene_names_as_rownames=F) {
    
    if(gene_names_as_rownames) {
        r_names <- rownames(x)
        x <-  resize_count_matrix(as.matrix(x), n, method=method)
        rownames(x) <- r_names
    } else {
        x <- as.data.frame(x)
        x[,-1] <-  resize_count_matrix(as.matrix(x[,-1]), n, method=method)
    }
    
    return(x)
    
}

#' Performs a correlation test and returns the specified 'val' from the output of cor.test
#' It returns NA in case there's an error

cor_test <- function(..., val) {
    
    tryCatch(cor.test(...)[[val]], error = function(e) NA)
    
}

get_tpm <- function(x, gene_lengths, reorder_gene_lengths=T) {
    
    x <- as.data.frame(x)
    
    if(! all(x[,1,drop=T] %in% names(gene_lengths)))
        stop('Lengths for all genes in the matrix were not found in the gene_lengths vector')
    
    if(reorder_gene_lengths)
        gene_lengths <- gene_lengths[x[,1,drop=T]]
    
    x[,-1] <- calculateTPM(as.matrix(x[,-1]), gene_lengths)
    
    return(x)
}

get_tmm <- function(x) {
    
    library('edgeR')
    
    x <- as.data.frame(x)
    
    dge <- DGEList(x[,-1])
    dge <- calcNormFactors(dge, method='TMM')
    x[,-1] <- cpm(dge)
    
    return(x)
}


########################################
# The functions below were taken from the scuttle package


#' Calculate TPMs
#'
#' Calculate transcripts-per-million (TPM) values for expression from feature-level counts.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size.factors A numeric vector containing size factors to adjust the library sizes.
#' If \code{NULL}, the library sizes are used directly. 
#' @param lengths Numeric vector providing the effective length for each feature in \code{x}.
#' Alternatively \code{NULL}, see Details.
#' @param assay.type A string specifying the assay of \code{x} containing the count matrix.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the ANY method, further arguments to pass to \code{\link{calculateCPM}}.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' This is because the number of UMIs is a direct (albeit biased) estimate of the number of transcripts.
#'
#' @return A numeric matrix of TPM values with the same dimensions as \code{x} (unless \code{subset.row} is specified).
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#' @seealso
#' \code{\link{calculateCPM}}, on which this function is based.
#'
#' @examples
#' example_sce <- mockSCE()
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' tout <- calculateTPM(example_sce, lengths = eff_len)
#' str(tout)
#'
#' @name calculateTPM
NULL

.calculate_tpm <- function(x, lengths=NULL, ...) {
    if (!is.null(lengths)) {
        x <- x/lengths
    }
    .calculate_cpm(x, ...)
}

#' @export
#' @rdname calculateTPM
setGeneric("calculateTPM", function(x, ...) standardGeneric("calculateTPM"))

#' @export
#' @rdname calculateTPM
setMethod("calculateTPM", "ANY", .calculate_tpm)


#' Calculate CPMs
#'
#' Calculate counts-per-million (CPM) values from the count data.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size.factors A numeric vector containing size factors to adjust the library sizes.
#' If \code{NULL}, the library sizes are used directly. 
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the count matrix.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param size_factors,subset_row,exprs_values Soft-deprecated counterparts to the arguments above.
#'
#' @details 
#' If \code{size.factors} are provided or available in \code{x}, they are used to define the effective library sizes. 
#' This is done by scaling all size factors such that the mean factor is equal to the mean sum of counts across all features. 
#' The effective library sizes are then used as the denominator of the CPM calculation.
#'
#' @return A numeric matrix of CPM values with the same dimensions as \code{x} (unless \code{subset.row} is specified).
#'
#' @name calculateCPM
#' @author Aaron Lun
#' @seealso 
#' \code{\link{normalizeCounts}}, on which this function is based.
#'
#' @examples
#' example_sce <- mockSCE()
#' cpm(example_sce) <- calculateCPM(example_sce)
#' str(cpm(example_sce))
NULL

#' @importFrom Matrix colSums
.calculate_cpm <- function(x, size.factors=NULL, subset.row=NULL, size_factors=NULL, subset_row=NULL) {
    size.factors <- .replace(size.factors, size_factors)
    subset.row <- .replace(subset.row, subset_row)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    lib.sizes <- colSums(x) / 1e6
    if (!is.null(size.factors)) {
        lib.sizes <- size.factors / mean(size.factors) * mean(lib.sizes)
    }

    normalizeCounts(x, size.factors=lib.sizes, log=FALSE, center.size.factors=FALSE)
}

#' @export
#' @rdname calculateCPM
setGeneric("calculateCPM", function(x, ...) standardGeneric("calculateCPM"))

#' @export
#' @rdname calculateCPM
setMethod("calculateCPM", "ANY", .calculate_cpm)

#' Compute normalized expression values
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A numeric matrix-like object containing counts for cells in the columns and features in the rows.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the count matrix.
#' @param size.factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add when log-transforming expression values.
#' @param center.size.factors Logical scalar indicating whether size factors should be centered at unity before being used.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param downsample Logical scalar indicating whether downsampling should be performed prior to scaling and log-transformation.
#' @param down.target Numeric scalar specifying the downsampling target when \code{downsample=TRUE}.
#' If \code{NULL}, this is defined by \code{down.prop} and a warning is emitted.
#' @param down.prop Numeric scalar between 0 and 1 indicating the quantile to use to define the downsampling target.
#' Only used when \code{downsample=TRUE}.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY or \linkS4class{DelayedMatrix} methods.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how library size factor calculations should be parallelized.
#' Only used if \code{size.factors} is not specified.
#' @param exprs_values,size_factors,pseudo_count,center_size_factors,subset_row,down_target,down_prop
#' Soft-deprecated equivalents to the arguments described previously.
#'
#' @details 
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This removes cell-specific scaling biases due to differences in sequencing coverage, capture efficiency or total RNA content.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{pseudo.count} to the normalized count and performing a log2-transformation.
#'
#' If no size factors are supplied, they are determined automatically from \code{x}:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} are supplied, they will override any size factors present in \code{x}.
#'
#' @section Centering the size factors:
#' If \code{center.size.factors=TRUE}, size factors are centred at unity prior to calculation of normalized expression values.
#' This ensures that the computed expression values can be interpreted as being on the same scale as original counts.
#' We can then compare abundances between features normalized with different sets of size factors; the most common use of this fact is in the comparison between spike-in and endogenous abundances when modelling technical noise (see \code{\link[scran]{modelGeneVarWithSpikes}} package for an example).
#'
#' More generally, when \code{log=TRUE}, centering of the size factors ensures that the value of \code{pseudo.count} can be interpreted as being on the same scale as the counts, i.e., the pseudo-count can actually be thought of as a \emph{count}.
#' This is important as it implies that the pseudo-count's impact will diminish as sequencing coverage improves.
#' Thus, if the size factors are centered, differences between log-normalized expression values will more closely approximate the true log-fold change with increasing coverage, whereas this would not be true of other metrics like log-CPMs with a fixed offset.
#'
#' The disadvantage of using centered size factors is that the expression values are not directly comparable across different calls to \code{\link{normalizeCounts}}, typically for multiple batches.
#' In theory, this is not a problem for metrics like the CPM, but in practice, we have to apply batch correction methods anyway to perform any joint analysis - see \code{\link[batchelor]{multiBatchNorm}} for more details. 
#'
#' @section Downsampling instead of scaling:
#' If \code{downsample=TRUE}, counts for each cell are randomly downsampled instead of being scaled.
#' This is occasionally useful for avoiding artifacts caused by scaling count data with a strong mean-variance relationship.
#' Each cell is downsampled according to the ratio between \code{down.target} and that cell's size factor.
#' (Cells with size factors below the target are not downsampled and are directly scaled by this ratio.)
#' If \code{log=TRUE}, a log-transformation is also performed after adding \code{pseudo.count} to the downsampled counts.
#'
#' We automatically set \code{down.target} to the 1st percentile of size factors across all cells involved in the analysis,
#' but this is only appropriate if the resulting expression values are not compared across different \code{normalizeCounts} calls.
#' To obtain expression values that are comparable across different \code{normalizeCounts} calls
#' (e.g., in \code{\link[scran]{modelGeneVarWithSpikes}} or \code{\link[batchelor]{multiBatchNorm}}),
#' \code{down_target} should be manually set to a constant target value that can be considered a low size factor in every call.
#'  
#' @return A numeric matrix-like object containing (log-)normalized expression values.
#' This has the same dimensions as \code{x} (unless \code{subset.row} is specified)
#' and is of the same class as the original count matrix.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{logNormCounts}}, which wraps this function for convenient use with SingleCellExperiment instances.
#'
#' \code{\link{librarySizeFactors}}, to compute the default size factors.
#'
#' \code{\link{downsampleMatrix}}, to perform the downsampling.
#' @examples
#' example_sce <- mockSCE()
#'
#' # Standard scaling + log-transformation:
#' normed <- normalizeCounts(example_sce)
#' normed[1:5,1:5]
#'
#' # Scaling without transformation:
#' normed <- normalizeCounts(example_sce, log=FALSE)
#' normed[1:5,1:5]
#'
#' # Downscaling with transformation:
#' normed <- normalizeCounts(example_sce, downsample=TRUE)
#' normed[1:5,1:5]
#'
#' # Using custom size factors:
#' with.meds <- computeMedianFactors(example_sce)
#' normed <- normalizeCounts(with.meds)
#' normed[1:5,1:5]
#' 
#' @name normalizeCounts
NULL

#' @export
#' @rdname normalizeCounts
setGeneric("normalizeCounts", function(x, ...) standardGeneric("normalizeCounts"))

#' @export
#' @rdname normalizeCounts
#' @importFrom BiocParallel SerialParam
setMethod("normalizeCounts", "ANY", function(x, size.factors=NULL, 
    log=TRUE, pseudo.count=1, center.size.factors=TRUE, subset.row=NULL,
    downsample=FALSE, down.target=NULL, down.prop=0.01, BPPARAM=SerialParam(),
    size_factors=NULL, pseudo_count=NULL, center_size_factors=NULL,
    subset_row=NULL, down_target=NULL, down_prop=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    size.factors <- .replace(size.factors, size_factors)
    center.size.factors <- .replace(center.size.factors, center_size_factors)
    down.target <- .replace(down.target, down_target)
    down.prop <- .replace(down.prop, down_prop)
    pseudo.count <- .replace(pseudo.count, pseudo_count)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (nrow(x)==0L) {
        return(x + 0) # coerce to numeric.
    }

    size.factors <- .get_default_sizes(x, size.factors, center.size.factors, BPPARAM=BPPARAM)
    if (length(size.factors)!=ncol(x)) {
        stop("number of size factors does not equal 'ncol(x)'")
    }
    if (!all(is.finite(size.factors) & size.factors > 0)) {
        stop("size factors should be positive")
    }

    if (downsample) {
        down.out <- .downsample_counts(x, size.factors, down.prop=down.prop, down.target=down.target)
        x <- down.out$x
        size.factors <- down.out$size.factors
    }

    .internal_transformer(x, size.factors, log, pseudo.count) 
})

.get_default_sizes <- function(x, size.factors, center.size.factors, ...) {
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, ...)
    }
    .center.size.factors(size.factors, center.size.factors)
}

.center.size.factors <- function(size.factors, center.size.factors) {
    if (center.size.factors) {
        size.factors <- size.factors/mean(size.factors)
    }
    size.factors
}

#' @importFrom stats quantile
.downsample_counts <- function(x, size.factors, down.prop, down.target) {
    if (is.null(down.target)) {
        down.target <- quantile(size.factors, probs=down.prop)
        warning("'down.target' defined as the 1st percentile of size factors")
    }
    down_rate <- pmin(1, down.target/size.factors)
    x <- downsampleMatrix(x, down_rate, bycol=TRUE)
    size.factors <- size.factors * down_rate/down.target
    list(x=x, size.factors=size.factors)
}



###########################################

setGeneric(".internal_transformer", function(x, ...) standardGeneric(".internal_transformer"))

#' @importFrom Matrix t
setMethod(".internal_transformer", "ANY", function(x, size.factors, log, pseudo.count) {
    norm_exprs <- t(t(x) / size.factors)
    if (log) {
        norm_exprs <- log2(norm_exprs + pseudo.count)
    }
    norm_exprs
})


.transform_sparse <- function(x, expanded_sf, log, pseudo.count) {
    x@x <- x@x/expanded_sf
    if (log) {
        x@x <- log2(x@x + pseudo.count)
    }
    x
}

.replace <- function(new, old) {
    if (!is.null(old)) old else new
}



