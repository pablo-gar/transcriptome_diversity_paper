source('../../R/misc.R')

cmdArgs <- commandArgs(T)
output_file <- cmdArgs[1]
table_files <- cmdArgs[-1]

output_table <- concatenate_table_files2(table_files, id_names=gsub('.txt', '', basename(table_files)))

write_tsv(output_table, output_file)

