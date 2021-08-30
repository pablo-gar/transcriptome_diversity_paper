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
