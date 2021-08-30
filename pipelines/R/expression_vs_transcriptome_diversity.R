do_nomarlized_expression Performs asociations between individual gene expression and average expression (in z-score) of a the entire genome (as averaged for all chromosome z-scores)
#' 
#' @param numeric - bonferroni pvalue cutoff to significant associations
#' @param character - path to file containing median expression per chromosome arm
#' @param character - path to output file containing signficant associations
#' @param character - path to output file containing non-signficant associations
#' @param character - path to output file having GO analyses of the significant associations
#'
##' @usage Rscript expression_vs_genome_instability_all_genome.R 

options(error=function()traceback(2))


if(!dir.exists("../../R"))
    stop("Can't find the shared R path, make sure to run script from directory where it lives")

source("../../R/gtex.R", chdir = T)
source("../../R/plots.R", chdir = T)
source("../../R/GO_analysis.R", chdir = T)
source("../../R/GWAS_methods.R", chdir = T)
source("../../R/misc.R", chdir = T)
source("../../R/mutationGeneAnnotation.R", chdir = T)

cmdArgs <- commandArgs(T)
pvalueCutoff <- as.numeric(cmdArgs[1])
do_exon_junctions <- as.logical(cmdArgs[2])
do_nomarlized_expression <- as.logical(cmdArgs[3])
do_covariates <- as.logical(cmdArgs[4])
tissue <- cmdArgs[5] # only important if do_nomarlized_expression
transcriptome_diversity_file <- cmdArgs[6]
outputSignificant <- cmdArgs[7]
outputNon <- cmdArgs[8]
outputGO <- cmdArgs[9]

if(do_exon_junctions & do_nomarlized_expression)
    stop('Specify just one exon junction mode')


# Covariates
vars_meta <- c('AGE', 'GENDER', 'RACE', 'BMI')
vars_sample_anno <- c('SMTSISCH', 'SMRIN')
sample_metadata_file <- file.path(GTEX_CON$root, GTEX_CON$sampleMetadataGTEXV8)

#transcriptome_diversity_file <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/transcriptome_diversity/Whole_Blood.txt"
#tissue <- 'Whole_Blood'
#do_exon_junctions <- F
#do_nomarlized_expression <- T
#pvalueCutoff <- 0.05
#outputNon <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/expression_associations/Whole_Blood/nonSignificant.txt.gz"
#outputSignificant <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/expression_associations/Whole_Blood/significant.txt.gz"
#outputGO <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/expression_associations/Whole_Blood/GO_results.txt"
#outputNon <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/expression_normalized_associations/Whole_Blood/nonSignificant.txt.gz"
#outputSignificant <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/expression_normalized_associations/Whole_Blood/significant.txt.gz"
#outputGO <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity/expression_normalized_associations/Whole_Blood/GO_results.txt"


#---------------------------
# MAIN

# Reading mutation data and metadata
median_exp_chr <- read.table(transcriptome_diversity_file, header=T, sep="\t", stringsAsFactors=F, check.names=F)
rownames(median_exp_chr) <- median_exp_chr[,1]
median_exp_chr <- t(median_exp_chr[,-1,drop=F])

# Read metadata
meta <- readMetadataGTEX(vars_meta)
meta <- meta[meta$RACE ==3,]
meta <- meta[,colnames(meta) != 'RACE']
sample_anno <- readSampleAnnotationGTEX(vars_sample_anno)
PCs <- readGenotypePCA()[1:3,]

# Select shared sampes
shared_samples <- colnames(median_exp_chr)[colnames(median_exp_chr) %in% rownames(sample_anno)]

median_exp_chr <- median_exp_chr[ ,shared_samples,drop=F]
sample_anno <- t(sample_anno[shared_samples, ])

shared_samples_long <- colnames(median_exp_chr)
colnames(median_exp_chr) <- gtexLongToShort(colnames(median_exp_chr))
colnames(sample_anno) <- gtexLongToShort(colnames(sample_anno))

# Reads expression file
if(do_exon_junctions) {
    expressionMat <- readExonJunctionCounts(shared_samples_long)
    rownames(expressionMat) <- expressionMat[,1]
} else if(do_nomarlized_expression) {
    expressionMat <- readNormalizedExpression(tissue)
    
    genes_ids <- gsub('\\..+', '', rownames(expressionMat))
    expressionMat <- expressionMat[!duplicated(genes_ids),]
    genes_ids <- gsub('\\..+', '', rownames(expressionMat))
    rownames(expressionMat) <- genes_ids
} else {
    expressionMat <- readAllGtexExpression(shared_samples_long, expressionFile=GTEX_CON$expresionAllTissuesV8)
    gene_ids <- gsub('\\..+', '', expressionMat[,1])
    expressionMat <- expressionMat[!duplicated(gene_ids), ]
    rownames(expressionMat) <- gsub('\\..+', '', expressionMat[,1])
}


# Gather shared samples
shared_samples_long <- shared_samples_long[colnames(median_exp_chr) %in% rownames(meta)]
shared_samples <- colnames(median_exp_chr)[colnames(median_exp_chr) %in% rownames(meta)]

shared_samples_long <- shared_samples_long[shared_samples %in% colnames(PCs)]
shared_samples <- shared_samples[shared_samples %in% colnames(PCs)]

if(!do_nomarlized_expression) {
    shared_samples <- shared_samples[shared_samples_long %in% colnames(expressionMat)]
    shared_samples_long <- shared_samples_long[shared_samples_long %in% colnames(expressionMat)]
    
    expressionMat <- expressionMat[, shared_samples_long]
} else {
    shared_samples_long <- shared_samples_long[shared_samples %in% colnames(expressionMat)]
    shared_samples <- shared_samples[shared_samples %in% colnames(expressionMat)]
    
    expressionMat <- expressionMat[, shared_samples]
    
}

median_exp_chr <- median_exp_chr[ ,shared_samples,drop=F]
sample_anno <- sample_anno[, shared_samples]
meta <- t(meta[shared_samples,])
PCs <- PCs[, shared_samples]


# Merge covariates
covariates <- rbind(sample_anno, meta, PCs)

# Select good expression
if(!do_nomarlized_expression) {
    expressionMat <- expressionMat[rowSums(expressionMat > 1) >= (ncol(expressionMat) * 0.2),]
    expressionMat <- as.data.frame(rankitNormalize(as.matrix(expressionMat)))

    inds_exp <- gtexLongToShort(colnames(expressionMat))
    expressionMat <- expressionMat[,!duplicated(inds_exp)]
    inds_exp <- inds_exp[!duplicated(inds_exp)]
    colnames(expressionMat) <- inds_exp
} else {
    expressionMat <- as.data.frame(rankitNormalize(as.matrix(expressionMat)))
}

################
# Goes over all mutation signatures

allSignatureGWAS <- list()

covariates <- as.matrix(covariates)
median_exp_chr <- as.matrix(median_exp_chr)

covariates <- covariates[apply(covariates, 1, function(x) length(unique(x)) > 1),]
covariates <- rankitNormalize(covariates, IND=1)

median_exp_chr <- rankitNormalize(median_exp_chr, IND=1)

# Normalizing by expression if looking at exon junctions
if(do_exon_junctions) {
    library('purrr')
    
    # Read original expression matrix
    expressionMat_for_normalization <- readAllGtexExpression(shared_samples_long)
    rownames(expressionMat_for_normalization) <- gsub('\\..+', '', expressionMat_for_normalization[,1])
    
    # Read sequencing depth
    sample_metadata <- readSampleAnnotationGTEX(columns='SMMPPDUN', location=sample_metadata_file)
    inds <- rownames(sample_metadata)[rownames(sample_metadata) %in% colnames(expressionMat_for_normalization)]
    
    sample_metadata <- sample_metadata[inds,,drop=F]
    expressionMat_for_normalization <- expressionMat_for_normalization[,inds]
    
    # Gathering individuals and genes of interest
    # exp mat
    inds_exp <- gtexLongToShort(colnames(expressionMat_for_normalization))
    expressionMat_for_normalization <- expressionMat_for_normalization[,!duplicated(inds_exp)]
    inds_exp <- inds_exp[!duplicated(inds_exp)]
    colnames(expressionMat_for_normalization) <- inds_exp
    
    # metadata
    inds_exp <- gtexLongToShort(rownames(sample_metadata))
    sample_metadata <- sample_metadata[!duplicated(inds_exp),,drop=F]
    inds_exp <- inds_exp[!duplicated(inds_exp)]
    rownames(sample_metadata) <- inds_exp
    
    # Correcting gene names
    expressionMat_for_normalization <- expressionMat_for_normalization[gsub('\\..+', '', rownames(expressionMat)), colnames(expressionMat)]
    
    # Do rankit normalization
    expressionMat_for_normalization <- as.data.frame(rankitNormalize(as.matrix(expressionMat_for_normalization)))
    sample_metadata <- as.data.frame(rankitNormalize(as.matrix(sample_metadata), 2))
    
    # Normalizing
    expressionMat_for_normalization <- t(expressionMat_for_normalization)
    expressionMat <- t(expressionMat)
    expressionMat <- expressionMat[rownames(expressionMat_for_normalization),]
    
    # Eliminate genes with NA expression
    expressionMat <- expressionMat[,colSums(is.na(expressionMat_for_normalization)) != nrow(expressionMat_for_normalization)]
    expressionMat_for_normalization <- expressionMat_for_normalization[,colSums(is.na(expressionMat_for_normalization)) != nrow(expressionMat_for_normalization)]
    
    norma <- map2_dfc(.x=as.data.frame(expressionMat), 
                      .y=as.data.frame(expressionMat_for_normalization), 
                      .sample_metadata=sample_metadata,
                      .f=function(x,y, .sample_metadata) resid(lm(x ~ ., data=cbind(data.frame(x=x, y=y), .sample_metadata), na.action=na.exclude))
                      )
    norma <- as.matrix(norma)
    
    dimnames(norma) <- dimnames(expressionMat)
    expressionMat <- t(norma)
}

# Exit if less than 50 individuals
if(ncol(median_exp_chr) < 50) {
    
    warning('Not enough individuals to do GWAS')
    all_GO_genes <- data.frame()
    all_GO <- data.frame()
    allHits <- data.frame()
    
    write.table(all_GO_genes, paste0(outputGO, "_genes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(all_GO, outputGO, sep = "\t", col.names = T, row.names = F, quote = F) 
    write.table(allHits, gzfile(outputSignificant), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(allHits, gzfile(outputNon), sep = "\t", col.names = T, row.names = F, quote = F)
    
    quit(save='no')
}


all_GO <- list()
all_GO_genes <- list()
for (i in rownames(median_exp_chr)) {
    
    
    phenoVector <- median_exp_chr[i,,drop=F]
    inds <- colnames(phenoVector)

    # Performing tests
    if(do_covariates) {
        x <- lmMatFunction(phenoVector, as.matrix(expressionMat[,inds]), cvrt = as.matrix(covariates[,inds]))
    } else {
        x <- lmMatFunction(phenoVector, as.matrix(expressionMat[,inds]), cvrt = as.matrix(covariates[,inds]))
    }
    
    x <- x$all$eqtls
    x$pvalue_bonf <- p.adjust(x$pvalue)
    x <- x[order(x$pvalue),]
    
    if(do_exon_junctions) {
        x$exon <- x$gene
        x$gene <- gsub("\\..+", "", x$gene)
    }     
    
    x$alias <- queryRefCDS(x$gene, reference = "gene_id", query = "gene_name")
    
    
    # Keeping results
    if(!exists("allHits")) {
        allHits <- x
    } else {
        allHits <- rbind(allHits, x)
    }
    
    # Skip if no significant results
    if (all(x$pvalue_bonf >= pvalueCutoff))
        next
    
    ###
    # Performing GO analyses
    if(!do_exon_junctions) {
        x_GO <- x
    } else {
        x_GO <- x[!duplicated(x$gene),]
    }
    
    GO_result <- performGO(setNames(x_GO$pvalue_bonf, x_GO$gene), gene_selection_function = function(p) p < pvalueCutoff)
    top_GO_categories <- get_GO_results(GO_result, fdr_signif=F, bonf_signif=F, n = 100)
        
    if(nrow(top_GO_categories) > 0) {
        
        toPrintGO_n <- ifelse(nrow(top_GO_categories) > 15, 15, nrow(top_GO_categories))
        
        top_GO_categories$instability_type <- i
        
        # Get genes in to GO category
        genes_in_GO <- get_genes_GO(x_GO$gene[x_GO$pvalue_bonf < pvalueCutoff], top_GO_categories, GO_result, toPrintGO_n)
        genes_in_GO$alias <- queryRefCDS(genes_in_GO$gene, reference = "gene_id", query = "gene_name")
        genes_in_GO$instability_type <-  i
        
        all_GO[[i]] <- top_GO_categories
        all_GO_genes[[i]] <- genes_in_GO
    }
    
    
    ####
    # Creating plots association plots
    phenoVectNorm <- resid(lm(t(phenoVector) ~ t(covariates[,inds]), na.action=na.exclude))
    counter <- 1
    significant <-  x[x$pvalue_bonf < pvalueCutoff,]
    for(k in 1:nrow(significant)) {
        if(do_exon_junctions) {
            j <- as.character(significant[k,"exon"])
        } else {
            j <- as.character(significant[k,"gene"])
        }
        alias_gene <-  significant[k,"alias"]
        if (counter > 10)
            break
        plot_title <- paste(c("instability_type", i, counter, "Gene", j, alias_gene), collapse = "_")
        p <- scatter(data.frame(y = phenoVectNorm, x = as.numeric(expressionMat[j,inds])), x = "x", y = "y", regression =T, method_cor='spearman')
        p <- p + ggtitle(plot_title) + xlab("gene expression") + ylab("transcriptome diversity")
        ggsave(file.path(dirname(outputSignificant), paste0("plot_", plot_title, ".pdf")), p)
        counter <- counter + 1 
    }

}


all_GO <- do.call(rbind, all_GO)
all_GO_genes <- do.call(rbind, all_GO_genes)

if(is.null(all_GO)) {
    all_GO <- data.frame()
    all_GO_genes <- data.frame()
} else {
    all_GO <- all_GO[order(-all_GO$fdr),]
}
    
# Saving instability files and covariates as well
median_exp_chr <- as.data.frame(t(median_exp_chr))
median_exp_chr$ind <- rownames(median_exp_chr)

covariates <- as.data.frame(t(covariates))
covariates$ind <- rownames(covariates)

expressionMat <- as.data.frame(t(expressionMat[unique(allHits[allHits$pvalue_bonf < pvalueCutoff, 'exon']),]))
expressionMat$ind <- rownames(expressionMat)

write.table(median_exp_chr, paste0(outputSignificant, "_instability.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(covariates, paste0(outputSignificant, "_covariates.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(expressionMat, paste0(outputSignificant, "_expression.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

write.table(all_GO_genes, paste0(outputGO, "_genes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(all_GO, outputGO, sep = "\t", col.names = T, row.names = F, quote = F) 
write.table(allHits[allHits$pvalue_bonf < pvalueCutoff,], gzfile(outputSignificant), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(allHits[allHits$pvalue_bonf >= pvalueCutoff,], gzfile(outputNon), sep = "\t", col.names = T, row.names = F, quote = F)
