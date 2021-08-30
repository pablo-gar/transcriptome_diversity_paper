# Convert isvdb variant file from keele to vcf files

library('dplyr')
library('tidyr')
library('purrr')


main <- function(cmdArgs=commandArgs(T)) {
    
    
    input_rds <- cmdArgs[1]
    output_vcf <- cmdArgs[2]
    
    #input_rds <- '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/keele/data_all/isvdb_var_db/19.rds'
    #output_vcf <- '19.vcf'
    
    chrom <- paste0('chr', gsub('\\..+', '', basename(input_rds)))
    
    
    # Read in rds file
    rds <- readRDS(input_rds)
    
    # Create vcf header
    vcf_header <- create_header(input_rds, chrom)
    
    # Create vcf matrix
    vcf_matrix <- create_matrix(rds, chrom)
    
    #Write file
    writeLines(vcf_header, output_vcf)
    write.table(vcf_matrix, output_vcf, sep='\t', quote=F, row.names=F, append=T)
    
    
}

create_header <- function(x, chrom) {
    
    header <- c('##fileformat=VCFv4.2',
                 '##fileDate=2020',
                 paste0('##contig=<ID=', chrom, '>'),
                 '##source=keele_et_al',
                 paste0('##reference=', x),
                 #'##INFO=<ID=GENE,Number=1,Type=String,Description="ENSEMBLE gene">',
                 #'##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="ENSEMBLE transcript">',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                 '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="conditional genotype quality, encoded as a phred quality−10log10 p(genotype call is wrong, conditioned on the site’s being variant)">'
    )

return(header)
    
}

create_matrix <- function(x, chrom) {
    
    for (i in 1:length(x)) 
        x[[i]]$sample_id <- names(x)[i]
    
    x <- do.call(rbind, x)
    
    x <- unique(select(x, -c(is_max, consequence_1, consequence_2, gene_name, transcript_name)))
    x$chrom <- chrom
    
    # Creating allele strings
    x <- x %>%
        group_by(variant_id,pos) %>%
        mutate(alleles_n = length(unique(c(allele_1, allele_2)))) %>%
        mutate(alleles = paste(names(sort(table(c(allele_1, allele_2)), decreasing=T)), collapse=',')) %>%
        ungroup() 
               
    # Createing genotypes
    alleles <- strsplit(x$alleles, ',')
    
    x$genotype <- paste0(map2_dbl(as.list(x$allele_1), alleles, function(x,y) which(x == y) - 1), 
                         '/',
                         map2_dbl(as.list(x$allele_2), alleles, function(x,y) which(x == y) - 1)
                         )
    x$ref <- map_chr(alleles, function(x) x[1])
    x$alt <- map_chr(alleles, function(x) paste(x[-1], collapse=','))
    
    # Creating left side matrix
    x_left <- x %>%
        select(chrom, pos, variant_id, ref, alt) %>%
        mutate(qual='30', filter='.', info = '.', format='GT:GQ') %>%
        unique() 
    
    colnames(x_left) <- c('#CHROM', 'POS', 'ID', 'REF', 'ATL', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
    
    # Creating right side of the matrix
    x_right <- x %>%
        mutate(phred = round(-10*log10(1-prob)), phred = ifelse(is.infinite(phred), 100, phred),
               dat=paste0(genotype, ':', phred)) %>%
        select(sample_id, variant_id, dat) %>%
        unique() %>%
        pivot_wider(values_from=dat, names_from=sample_id) %>%
        rename(ID='variant_id')
    
    x_right[is.na(x_right)] <- './.:0'
    
    # Merge
    x <- full_join(x_left, x_right)
    
    return(x)
    
}

main()
