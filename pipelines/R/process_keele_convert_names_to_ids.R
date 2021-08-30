
main <- function(cmdArgs=commandArgs(T)) {
    
    input_file <- cmdArgs[1]
    converstion_table <- cmdArgs[2]
    output_file <- cmdArgs[3]
    
    #input_file <- '/scratch/users/paedugar/cnv_gtex_project/expression_matrices/keele_kidney.txt'
    #converstion_table <- '/scratch/users/paedugar/cnv_gtex_project/auxiliary_files/Mus_musculus.NCBIM37.67.id_converstion_table.txt'
    #output_file <- 'deleteme.txt'
    
    exp_mat <- read.table(input_file, sep='\t', stringsAsFactors=F, header=T, check.names=F)
    conversion_vector <- read.table(converstion_table, sep='\t', stringsAsFactors=F, header=T)
    conversion_vector <- conversion_vector[!duplicated(conversion_vector[,2]),]
    conversion_vector <- setNames(conversion_vector[,1], conversion_vector[,2])
    
    exp_mat[,1] <- conversion_vector[exp_mat[,1]]
    exp_mat <- exp_mat[!is.na(exp_mat[,1]),]
    
    write.table(exp_mat, output_file, sep='\t', quote=F, row.names=F, col.names=T)
    
}

main()
