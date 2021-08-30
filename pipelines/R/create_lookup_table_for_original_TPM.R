library('dplyr')
library('readr')

main <- function(cmdArgs = commandArgs(T)) {
    
    out_file <- cmdArgs[1]
    files <- cmdArgs[-1]
    
    tables <- list()
    
    for(f in files) 
        tables[[f]] <- read_tsv(f, col_types = cols(.default = "c"))
    
    tables <- do.call(bind_rows, tables)
    tables <- tables[tables$count_type=='TPM',]
    tables <- tables[!is.na(tables$count_type),]
        
    write_tsv(tables, out_file)
}

main()
