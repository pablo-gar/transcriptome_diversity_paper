f Creates a sample table with all the files that exist in the bam folder

library('rjson')
source('../../R/gtex.R', chdir=T)
CONFIG <- fromJSON(file='../../config.json')

# Read first line of gene expression table to get 
samples <- unlist(strsplit(readLines(GTEX_CON$expresionAllTissues, n=3)[3], '\t'))[-(1:2)]

# Get the sraIds
conversion_table <- read.table(file.path(GTEX_CON$root, GTEX_CON$sraTable_v8), heade=T, sep='\t', stringsAsFactors=F)
conversion_table <- conversion_table[,c('Run', 'body_site', 'submitted_subject_id', 'Sample.Name')]
conversion_table <- conversion_table[!duplicated(conversion_table$Sample.Name),]
rownames(conversion_table) <- conversion_table$Sample.Name


files <- list.files(file.path(CONFIG$projectDir, CONFIG$bam_files$dir), full.names=T, recursive=T, pattern='bam$')
group <- basename(dirname(files))
sample_id <- gsub('(SRR.+)_.+', '\\1', basename(files))
names(files) <- sample_id


out <- data.frame(sample_id=conversion_table[samples, 'Run'], 
                  group=conversion_table[samples, 'body_site'], 
                  path=files[conversion_table[samples, 'Run']], 
                  gtex_id_long= samples,
                  gtex_id=conversion_table[samples,'submitted_subject_id'], stringsAsFactors=F)

out$group <- code_friendly_tissues(out$group)

write.table(out, file.path(CONFIG$projectDir, CONFIG$auxiliary_Files$dir, 'sample_table_simple_calling_for_expression.txt'), sep='\t', quote=F, row.names=F, col.names=T)
