##################################################################################
##################################################################################
###
###     File S14: Example code for using miQTL as employed
###               in Keele & Quach et al
###     -- Multi-tissue analysis of gene expression and chromatin accessibility
###        in 47 Collaborative Cross strains
###
###
#keele_{keele_tissue}##     Tasks: 1. QTL scans
###            2. Variant association
###            3. QTL scan permutations
###            4. Conditional QTL scans
###            5. Mediation
###            6. Mediation permutations
###            7. Example of large scale analysis
###     
###     Author: Greg Keele
###
##################################################################################
##################################################################################
#### Necessary R packages
library("tidyverse")
library("readr")
library("miqtl")

main <- function(cmdArgs=commandArgs(T)) {
    
    trans_file <- cmdArgs[1]
    expression_file <-  cmdArgs[2]
    gene_file <- cmdArgs[3]
    genome_cache_dir <- cmdArgs[4]
    var_dir <- cmdArgs[5]
    covariate_mode <- cmdArgs[6] #"transcriptome_diversity|batch|transcriptome_diversity+batch|none"
    out_file <- cmdArgs[7]
    out_file_top <- cmdArgs[8]
    
    #trans_file <- "/scratch/users/paedugar/cnv_gtex_project/transcriptome_diversity_other_datasets/transcriptome_diversity_as_covariates/keele_liver.txt"
    #expression_file <-  "/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/liver_expression.csv.zip"
    #gene_file <- "/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/refseq_mm9_tss.txt.zip"
    #genome_cache_dir <- "/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/cc_genome_cache_full_l2_0.1/"
    #var_dir <- "/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/isvdb_var_db/"
    #covariate_mode <- "transcriptome_diversity+batch"
    #out_file <- "./delemete_pvalues.txt"
    #out_file_top <- "./delemete_pvalues_top.txt"
    
    
    # Loading data 
    gene_dat <- readr::read_tsv(gene_file, col_names = c("gene_chr", "gene_start", "gene_end", "gene", "alt_name")) %>%
        mutate(gene_start = gene_start/1e6, gene = gsub(x = gene, pattern = "-", replacement = ".")) %>%
        as.data.frame()
    
    expression_dat <- read_csv(expression_file) %>%
      mutate(BATCH = factor(BATCH)) %>%
      as.data.frame
  
    trans <- read_tsv(trans_file) %>%
        pivot_longer(-ID, names_to='SUBJECT.NAME', values_to='transcriptome_diversity') %>%
        select(-ID)

    expression_dat <- left_join(expression_dat, trans)
    
    # Perform one example of haplotype eqtl to obtain number of tests
    cox7c_scan <- scan.h2lmm(genomecache = genome_cache_dir,
                                   data = expression_dat, 
                                   model = "additive",
                                   formula = rint(gene_Cox7c) ~ 1 + BATCH, 
                                   p.value.method = "ANOVA",
                                   use.multi.impute = FALSE)
    
    # Perform QR "transcriptome_diversity|batch|transcriptome_diversity+batch|none"
    
    if(covariate_mode=='batch') {
        current_qr <- extract.qr(genomecache = genome_cache_dir, 
                               data = expression_dat, 
                               formula = ~ 1 + BATCH)
    } else if(covariate_mode=='transcriptome_diversity') {
        current_qr <- extract.qr(genomecache = genome_cache_dir, 
                               data = expression_dat, 
                               formula = ~ 1 + transcriptome_diversity)
    } else if(covariate_mode=='transcriptome_diversity+batch') {
        current_qr <- extract.qr(genomecache = genome_cache_dir, 
                               data = expression_dat, 
                               formula = ~ 1 + BATCH + transcriptome_diversity)
    } else if(covariate_mode=='none') {
        current_qr <- extract.qr(genomecache = genome_cache_dir, 
                               data = expression_dat, 
                               formula = ~ 1)
    }
    
    ###########################
    # Perform genome-wide scan
    genes <- grep(x = names(expression_dat), pattern = "gene", value = TRUE)

    ## Set the number of genes to 100 for this example
    m <- length(genes)
    #m <- 3

    ## Initializing p-value matrix
    pval_mat <- matrix(NA, nrow = length(cox7c_scan$p.value), ncol = m)
    ## Setting up eQTL scan loop
    for (i in 1:m) {
      this_scan <- scan.qr(qr.object = current_qr, 
                           data = expression_dat, 
                           phenotype = paste0("rint(", genes[i], ")"))
      pval_mat[,i] <- this_scan$p.value
    }
    rownames(pval_mat) <- this_scan$loci
    colnames(pval_mat) <- genes[1:m]
    
    # Getting top pvalues per gene
    pval_mat_top <- data.frame(variant_id=rep('', ncol(pval_mat)), gene_id='', nominal_pvalue=0, stringsAsFactors=F)
    for(i in 1:ncol(pval_mat)) {
        min_i <- which.min(pval_mat[,i,drop=T])
        pval_mat_top[i,] <- c(rownames(pval_mat)[i], colnames(pval_mat)[i], pval_mat[min_i,i])
    }
    
    
    pval_mat <- cbind(data.frame(variant_id=rownames(pval_mat), stringsAsFactors=F), pval_mat)
    ###
    # Saving results
    write_tsv(as.data.frame(pval_mat), out_file)
    write_tsv(as.data.frame(pval_mat_top), out_file_top)
    
}

##################################################################################
##################################################################################
###
###     File S13: Functions used in Keele & Quach et al to produce 
###               supplemental tables
###     -- Multi-tissue analysis of gene expression and chromatin accessibility
###        in 47 Collaborative Cross strains
###     
###     Author: Greg Keele
###
##################################################################################
##################################################################################
### Local eQTL
pull_local_eqtl <- function(method1_dat, method2_dat, method3_dat, gene_dat) {
  method1_table <- method1_dat %>% 
    dplyr::filter(eqtl_qval < 0.2) %>%
    dplyr::filter(eqtl_status == "local") %>%
    dplyr::select(gene, gene_start, gene_chr, eqtl_chr, eqtl_pos, eqtl_locus, fixef_eqtl_effect_size, ranef_eqtl_effect_size, eqtl_qval, eqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(eqtl_qval_genomewide = eqtl_qval,
                  eqtl_pval_genomewide = eqtl_pval) %>%
    dplyr::mutate(gene_chr = as.character(gene_chr),
                  eqtl_chr = as.character(eqtl_chr))
  method2_table <- method2_dat %>% 
    dplyr::filter(eqtl_qval < 0.2) %>%
    dplyr::filter(eqtl_status == "local") %>%
    dplyr::select(gene, gene_start, gene_chr, fixef_eqtl_effect_size, ranef_eqtl_effect_size, eqtl_pos, eqtl_locus, eqtl_qval, eqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(eqtl_qval_chrwide = eqtl_qval,
                  eqtl_pval_chrwide = eqtl_pval) %>%
    dplyr::mutate(gene_chr = as.character(gene_chr),
                  eqtl_chr = gene_chr)
  method3_table <- method3_dat %>%
    dplyr::filter(eqtl_pval_chrwide < 0.05) %>%
    dplyr::select(gene, gene_start, gene_chr, eqtl_pos, eqtl_locus, fixef_eqtl_effect_size, ranef_eqtl_effect_size, eqtl_pval_genomewide, eqtl_pval_chrwide, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(eqtl_pval_method3genomewide = eqtl_pval_genomewide,
                  eqtl_pval_method3chrwide = eqtl_pval_chrwide) %>%
    dplyr::mutate(gene_chr = as.character(gene_chr),
                  eqtl_chr = gene_chr)
  local_eqtl_table <- dplyr::full_join(dplyr::full_join(method1_table,
                                                        method2_table),
                                       method3_table) %>%
    dplyr::mutate(gene = gsub(x=gene, pattern="gene_", replacement=""))# %>%

  ## Adding gene data
  local_eqtl_table <- gene_dat %>% 
    dplyr::select(gene, gene_end) %>%
    dplyr::mutate(gene = gsub(x = gene, pattern="_|-", replacement="."), 
                  gene_end = gene_end/1e6) %>%
    dplyr::right_join(local_eqtl_table) %>%
    dplyr::select(gene, everything())
  
  ## Detection method
  local_eqtl_table$detection_method <- NA
  mod_index <- !is.na(local_eqtl_table$eqtl_qval_genomewide) & local_eqtl_table$eqtl_qval_genomewide <= 0.1 & is.na(local_eqtl_table$detection_method)
  local_eqtl_table$detection_method[mod_index] <- "method 1"
  mod_index <- !is.na(local_eqtl_table$eqtl_qval_chrwide) & local_eqtl_table$eqtl_qval_chrwide <= 0.1 & is.na(local_eqtl_table$detection_method)
  local_eqtl_table$detection_method[mod_index] <- "method 2"
  mod_index <- !is.na(local_eqtl_table$eqtl_pval_method3genomewide) & local_eqtl_table$eqtl_pval_method3genomewide <= 0.05 & is.na(local_eqtl_table$detection_method)
  local_eqtl_table$detection_method[mod_index] <- "method 3 genomewide"
  mod_index <- !is.na(local_eqtl_table$eqtl_pval_method3chrwide) & local_eqtl_table$eqtl_pval_method3chrwide <= 0.05 & is.na(local_eqtl_table$detection_method)
  local_eqtl_table$detection_method[mod_index] <- "method 3 chrwide"
  mod_index <- !is.na(local_eqtl_table$eqtl_qval_genomewide) & local_eqtl_table$eqtl_qval_genomewide <= 0.2 & is.na(local_eqtl_table$detection_method)
  local_eqtl_table$detection_method[mod_index] <- "method 1 q<0.2"
  mod_index <- !is.na(local_eqtl_table$eqtl_qval_chrwide) & local_eqtl_table$eqtl_qval_chrwide <= 0.2 & is.na(local_eqtl_table$detection_method)
  local_eqtl_table$detection_method[mod_index] <- "method 2 q<0.2"
  
  local_eqtl_table
}


### Distal eQTL
pull_distal_eqtl <- function(method1_dat, method2_dat, gene_dat) {
  method1_table <- method1_dat %>% 
    dplyr::filter(eqtl_qval < 0.2) %>%
    dplyr::filter(eqtl_status == "distal") %>%
    dplyr::select(gene, gene_start, gene_chr, eqtl_chr, eqtl_pos, eqtl_locus, 
                  fixef_eqtl_effect_size, ranef_eqtl_effect_size,
                  eqtl_qval, eqtl_pval, paste("ranef", LETTERS[1:8], sep = "_"),
                  exmediation_pval, exmediator_pos, fixef_eqtl_effect_size_mediator) %>%
    dplyr::rename(eqtl_qval_genomewide = eqtl_qval,
                  eqtl_pval_genomewide = eqtl_pval) %>%
    dplyr::mutate(gene_chr = as.character(gene_chr),
                  eqtl_chr = as.character(eqtl_chr))
  method2_table <- method2_dat %>% 
    dplyr::filter(eqtl_qval < 0.2) %>%
    dplyr::filter(eqtl_status == "distal") %>%
    dplyr::select(gene, gene_start, gene_chr, fixef_eqtl_effect_size, ranef_eqtl_effect_size, eqtl_pos, eqtl_locus, eqtl_qval, eqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(eqtl_qval_chrwide = eqtl_qval,
                  eqtl_pval_chrwide = eqtl_pval) %>%
    dplyr::mutate(gene_chr = as.character(gene_chr),
                  eqtl_chr = gene_chr)
  distal_eqtl_table <- dplyr::full_join(method1_table,
                                        method2_table) %>%
    dplyr::mutate(gene = gsub(x=gene, pattern="gene_", replacement=""))# %>%
  
  ## Adding gene data
  distal_eqtl_table <- gene_dat %>% 
    dplyr::select(gene, gene_end) %>%
    dplyr::mutate(gene = gsub(x = gene, pattern="_|-", replacement="."), 
                  gene_end = gene_end/1e6) %>%
    dplyr::right_join(distal_eqtl_table) %>%
    dplyr::select(gene, everything())
  
  ## Detection method
  distal_eqtl_table$detection_method <- NA
  mod_index <- !is.na(distal_eqtl_table$eqtl_qval_genomewide) & distal_eqtl_table$eqtl_qval_genomewide <= 0.1 & is.na(distal_eqtl_table$detection_method)
  distal_eqtl_table$detection_method[mod_index] <- "method 1"
  mod_index <- !is.na(distal_eqtl_table$eqtl_qval_chrwide) & distal_eqtl_table$eqtl_qval_chrwide <= 0.1 & is.na(distal_eqtl_table$detection_method)
  distal_eqtl_table$detection_method[mod_index] <- "method 2"
  mod_index <- !is.na(distal_eqtl_table$eqtl_qval_genomewide) & distal_eqtl_table$eqtl_qval_genomewide <= 0.2 & is.na(distal_eqtl_table$detection_method)
  distal_eqtl_table$detection_method[mod_index] <- "method 1 q<0.2"
  mod_index <- !is.na(distal_eqtl_table$eqtl_qval_chrwide) & distal_eqtl_table$eqtl_qval_chrwide <= 0.2 & is.na(distal_eqtl_table$detection_method)
  distal_eqtl_table$detection_method[mod_index] <- "method 2 q<0.2"
  
  distal_eqtl_table
}

### Local cQTL
pull_local_cqtl <- function(method1_dat, method2_dat, method3_dat) {
  method1_table <- method1_dat %>% 
    dplyr::filter(cqtl_qval < 0.2) %>%
    dplyr::filter(cqtl_status == "local") %>%
    dplyr::select(chromatin, chromatin_chr, chromatin_pos, cqtl_chr, cqtl_pos, cqtl_locus, fixef_cqtl_effect_size, ranef_cqtl_effect_size, cqtl_qval, cqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(cqtl_qval_genomewide = cqtl_qval,
           cqtl_pval_genomewide = cqtl_pval) %>%
    dplyr::mutate(chromatin_chr = as.character(chromatin_chr),
           cqtl_chr = as.character(cqtl_chr))
  method2_table <- method2_dat %>% 
    dplyr::filter(cqtl_qval < 0.2) %>%
    dplyr::filter(cqtl_status == "local") %>%
    dplyr::select(chromatin, chromatin_chr, chromatin_pos, fixef_cqtl_effect_size, ranef_cqtl_effect_size, cqtl_pos, cqtl_locus, cqtl_qval, cqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(cqtl_qval_chrwide = cqtl_qval,
           cqtl_pval_chrwide = cqtl_pval) %>%
    dplyr::mutate(chromatin_chr = as.character(chromatin_chr),
                  cqtl_chr = chromatin_chr)
  method3_table <- method3_dat %>%
    dplyr::filter(cqtl_pval_chrwide < 0.05) %>%
    dplyr::select(chromatin, chromatin_chr, chromatin_pos, cqtl_pos, cqtl_locus, fixef_cqtl_effect_size, ranef_cqtl_effect_size, cqtl_pval_genomewide, cqtl_pval_chrwide, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(cqtl_pval_method3genomewide = cqtl_pval_genomewide,
                  cqtl_pval_method3chrwide = cqtl_pval_chrwide) %>%
    dplyr::mutate(chromatin_chr = as.character(chromatin_chr),
                  cqtl_chr = chromatin_chr)
  local_cqtl_table <- dplyr::full_join(dplyr::full_join(method1_table,
                                                        method2_table),
                                       method3_table)
  
  ## Detection method
  local_cqtl_table$detection_method <- NA
  mod_index <- !is.na(local_cqtl_table$cqtl_qval_genomewide) & local_cqtl_table$cqtl_qval_genomewide <= 0.1 & is.na(local_cqtl_table$detection_method)
  local_cqtl_table$detection_method[mod_index] <- "method 1"
  mod_index <- !is.na(local_cqtl_table$cqtl_qval_chrwide) & local_cqtl_table$cqtl_qval_chrwide <= 0.1 & is.na(local_cqtl_table$detection_method)
  local_cqtl_table$detection_method[mod_index] <- "method 2"
  mod_index <- !is.na(local_cqtl_table$cqtl_pval_method3genomewide) & local_cqtl_table$cqtl_pval_method3genomewide <= 0.05 & is.na(local_cqtl_table$detection_method)
  local_cqtl_table$detection_method[mod_index] <- "method 3 genomewide"
  mod_index <- !is.na(local_cqtl_table$cqtl_pval_method3chrwide) & local_cqtl_table$cqtl_pval_method3chrwide <= 0.05 & is.na(local_cqtl_table$detection_method)
  local_cqtl_table$detection_method[mod_index] <- "method 3 chrwide"
  mod_index <- !is.na(local_cqtl_table$cqtl_qval_genomewide) & local_cqtl_table$cqtl_qval_genomewide <= 0.2 & is.na(local_cqtl_table$detection_method)
  local_cqtl_table$detection_method[mod_index] <- "method 1 q<0.2"
  mod_index <- !is.na(local_cqtl_table$cqtl_qval_chrwide) & local_cqtl_table$cqtl_qval_chrwide <= 0.2 & is.na(local_cqtl_table$detection_method)
  local_cqtl_table$detection_method[mod_index] <- "method 2 q<0.2"
  
  local_cqtl_table
}

### Distal cQTL
pull_distal_cqtl <- function(method1_dat, method2_dat) {
  method1_table <- method1_dat %>% 
    dplyr::filter(cqtl_qval < 0.2) %>%
    dplyr::filter(cqtl_status == "distal") %>%
    dplyr::select(chromatin, chromatin_chr, chromatin_pos, cqtl_chr, cqtl_pos, cqtl_locus, fixef_cqtl_effect_size, ranef_cqtl_effect_size, cqtl_qval, cqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(cqtl_qval_genomewide = cqtl_qval,
                  cqtl_pval_genomewide = cqtl_pval) %>%
    dplyr::mutate(chromatin_chr = as.character(chromatin_chr),
                  cqtl_chr = as.character(cqtl_chr))
  method2_table <- method2_dat %>% 
    dplyr::filter(cqtl_qval < 0.2) %>%
    dplyr::filter(cqtl_status == "distal") %>%
    dplyr::select(chromatin, chromatin_chr, chromatin_pos, fixef_cqtl_effect_size, ranef_cqtl_effect_size, cqtl_pos, cqtl_locus, cqtl_qval, cqtl_pval, paste("ranef", LETTERS[1:8], sep = "_")) %>%
    dplyr::rename(cqtl_qval_chrwide = cqtl_qval,
                  cqtl_pval_chrwide = cqtl_pval) %>%
    dplyr::mutate(chromatin_chr = as.character(chromatin_chr),
                  cqtl_chr = chromatin_chr)
  distal_cqtl_table <- dplyr::full_join(method1_table,
                                        method2_table)
  
  ## Detection method
  distal_cqtl_table$detection_method <- NA
  mod_index <- !is.na(distal_cqtl_table$cqtl_qval_genomewide) & distal_cqtl_table$cqtl_qval_genomewide <= 0.1 & is.na(distal_cqtl_table$detection_method)
  distal_cqtl_table$detection_method[mod_index] <- "method 1"
  mod_index <- !is.na(distal_cqtl_table$cqtl_qval_chrwide) & distal_cqtl_table$cqtl_qval_chrwide <= 0.1 & is.na(distal_cqtl_table$detection_method)
  distal_cqtl_table$detection_method[mod_index] <- "method 2"
  mod_index <- !is.na(distal_cqtl_table$cqtl_qval_genomewide) & distal_cqtl_table$cqtl_qval_genomewide <= 0.2 & is.na(distal_cqtl_table$detection_method)
  distal_cqtl_table$detection_method[mod_index] <- "method 1 q<0.2"
  mod_index <- !is.na(distal_cqtl_table$cqtl_qval_chrwide) & distal_cqtl_table$cqtl_qval_chrwide <= 0.2 & is.na(distal_cqtl_table$detection_method)
  distal_cqtl_table$detection_method[mod_index] <- "method 2 q<0.2"
  
  distal_cqtl_table
}

## Local eQTL with chromatin mediators
pull_chromatin_mediation <- function(method3_dat, gene_dat) {
  mediation_table <- method3_dat %>%
    dplyr::filter(grepl(x=pass_status, pattern="e:l")) %>%
    dplyr::filter(grepl(x=pass_status, pattern="c:l")) %>%
    dplyr::filter(grepl(x=pass_status, pattern="m:l")) %>%
    dplyr::filter(mediator_status == "local") %>%
    dplyr::filter(fixef_eqtl_effect_size <= fixef_cqtl_effect_size) %>%
    dplyr::select(gene, gene_start, gene_chr, 
           eqtl_pos, eqtl_locus, fixef_eqtl_effect_size, ranef_eqtl_effect_size, eqtl_pval_genomewide, eqtl_pval_chrwide, paste("ranef", LETTERS[1:8], sep = "_"),  
           cqtl_pos, cqtl_locus, cqtl_locus, fixef_cqtl_effect_size, ranef_cqtl_effect_size, cqtl_pval_genomewide, cqtl_pval_chrwide,
           mediator_outcome, mediator_pos, mediation_pval_genomewide, mediation_pval_chrwide) %>%
    dplyr::rename(chromatin_mediator = mediator_outcome) %>%
    dplyr::mutate(gene = gsub(x=gene, pattern="gene_", replacement=""))
  ## Adding gene data
  mediation_table <- gene_dat %>% 
    dplyr::select(V3, V4) %>%
    dplyr::rename(gene_end = V3,
           gene = V4) %>%
    dplyr::mutate(gene = gsub(x = gene, pattern="_|-", replacement="."), 
           gene_end = gene_end/1e6) %>%
    dplyr::right_join(mediation_table) %>%
    dplyr::select(gene, everything())
  ## Adding chromatin information
  mediation_table$chromatin_start <- sapply(1:length(mediation_table$chromatin_mediator),
                                            function(i) as.numeric(unlist(strsplit(x=as.character(mediation_table$chromatin_mediator[i]), split=".", fixed=TRUE))[2])/1e6)
  mediation_table$chromatin_end <- sapply(1:length(mediation_table$chromatin_mediator),
                                          function(i) as.numeric(unlist(strsplit(x=as.character(mediation_table$chromatin_mediator[i]), split=".", fixed=TRUE))[3])/1e6)
  mediation_table
}

### eQTL
summarize_all_eqtl_for_table <- function(qtl_dat, expression_dat) {
  n <- length(unique(qtl_dat$gene))
  prop <- n/sum(grepl(x = names(expression_dat), pattern = "gene"))
  cat(round(n, 1), paste0("(", round(prop*100, 1), ")"))
}
summarize_local_eqtl_for_table <- function(qtl_dat) {
  n <- qtl_dat %>%
    dplyr::filter(eqtl_status == "local") %>%
    dplyr::select(gene) %>%
    unique %>% 
    nrow
  prop <- n/length(unique(qtl_dat$gene))
  cat(round(n, 1), paste0("(", round(prop*100, 1), ")"))
}
summarize_distal_eqtl_for_table <- function(qtl_dat) {
  n <- qtl_dat %>%
    dplyr::filter(eqtl_status == "distal") %>%
    dplyr::select(gene) %>%
    unique %>% 
    nrow
  prop <- n/length(unique(qtl_dat$gene))
  cat(round(n, 1), paste0("(", round(prop*100, 1), ")"))
}

### cQTL
summarize_all_cqtl_for_table <- function(qtl_dat, chromatin_dat) {
  n <- length(unique(qtl_dat$chromatin))
  prop <- n/sum(grepl(x = names(chromatin_dat), pattern = "chr"))
  cat(round(n, 1), paste0("(", round(prop*100, 1), ")"))
}
summarize_local_cqtl_for_table <- function(qtl_dat) {
  n <- qtl_dat %>%
    filter(cqtl_status == "local") %>%
    select(chromatin) %>%
    unique %>% 
    nrow
  prop <- n/length(unique(qtl_dat$chromatin))
  cat(round(n, 1), paste0("(", round(prop*100, 1), ")"))
}
summarize_distal_cqtl_for_table <- function(qtl_dat) {
  n <- qtl_dat %>%
    filter(cqtl_status == "distal") %>%
    select(chromatin) %>%
    unique %>% 
    nrow
  prop <- n/length(unique(qtl_dat$chromatin))
  cat(round(n, 1), paste0("(", round(prop*100, 1), ")"))
}

##### Chromatin mediation
summarize_mediation_for_table <- function(qtl_dat) {
  n <- length(unique(qtl_dat$gene))
  n
}

main()
