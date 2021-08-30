library('stringr')

main <- function(cmdArgs=commandArgs(T)) {
    
    out_file <- cmdArgs[1]
    mat_files <- cmdArgs[-1]
    
    results <- data.frame(id='', name='', path=mat_files, group='sarantopoulou', subgroup='sarantopoulou', method='', count_type='', tissue='liver', stringsAsFactors=F)
    for (i in 1:length(mat_files)) {
        
        id = str_remove(basename(mat_files[i]), '.txt')
        splitted <- str_split(id, '_')[[1]]
        
        results[i, c('id', 'name', 'method', 'count_type')] <- c(id, id, splitted[2], toupper(splitted[3]))
        
    }
    
    write.table(results, out_file, sep='\t', quote=F, row.names=F)
    
    
}

main()
