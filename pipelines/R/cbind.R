library('dplyr')
library('readr')

main <- function(cmdArgs = commandArgs(T)) {
    
    out_file <- cmdArgs[1]
    files <- cmdArgs[-1]
    
    tables <- list()
    
    for(f in files) 
        tables[[f]] <- read_tsv(f)
    
    tables <- do.call(bind_cols, tables)
        
    write_tsv(tables, out_file)
}

main()
