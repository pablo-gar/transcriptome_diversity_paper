library('stringr')

main <- function(cmdArgs=commandArgs(T)) {
    
    out_file <- cmdArgs[1]
    mat_files <- cmdArgs[-1]
    
    results <- data.frame(id='', name='', path=mat_files, group='nagarajan', subgroup='nagarajan', method='nagarajan', count_type='', tissue='cell_line_MCF7', stringsAsFactors=F)
    for (i in 1:length(mat_files)) {
        
        id = str_remove(basename(mat_files[i]), '.txt')
        splitted <- str_split(id, '_')[[1]]
        
        results[i, c('id', 'name', 'count_type')] <- c(id, id, toupper(splitted[2]))
        
    }
    
    write.table(results, out_file, sep='\t', quote=F, row.names=F)
    
    
}

main()
