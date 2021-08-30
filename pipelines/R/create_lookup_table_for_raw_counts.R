library('stringr')

cmdArgs <- commandArgs(T)
root_dir <- cmdArgs[1]
exp_dir <- cmdArgs[2]
aux_dir <- cmdArgs[3]
dat_dir <- cmdArgs[4]
out_file <- cmdArgs[5]

current_files <- list.files(file.path(root_dir, exp_dir), full.names=T)
current_files_relative <- file.path(exp_dir, basename(current_files))


dats <- data.frame(id=gsub('\\..+', '', basename(current_files)), stringsAsFactors=F)
dats$name <-  dats$id
dats$path <- current_files_relative
dats$group <- str_split_fixed(dats$id, '_', 2)[,1]
dats$subgroup <- NA
dats$method <- NA
dats$count_type <- 'read_counts'
dats$tissue <- NA
dats$gtf_file <- NA


# Eliminate datasets that we won't use, i.e. those that do not have raw counts
dats <- dats[dats$group != 'arora',]
dats <- dats[!grepl('TCGA', dats$group),]
dats <- dats[!grepl('tmm', dats$id),]
dats <- dats[!grepl('tpm', dats$id),]
dats <- dats[!grepl('cov', dats$id),]
dats <- dats[!grepl('temp', dats$id),]

# remove tbi for keele
dats <- dats[!grepl('\\.tbi$', dats$path),]
dats <- dats[!grepl('\\.bed.gz$', dats$path),]

# For lin et al select lin_counts
dats <- dats[dats$id != 'lin',]

# For gtex keep raw counts
dats <- dats[!(dats$group=='gtex' & !grepl('counts', dats$id)),]

# Remove the word "counts" from the ids and names
dats$id <- gsub('_counts', '', dats$id)
dats$name <- gsub('_counts', '', dats$name)

# Add gtf info
dats$gtf_file [dats$group == 'gtex'] <- file.path(dat_dir, 'gtex', 'gencode.v26.GRCh38.genes.gtf')
dats$gtf_file [dats$group == 'keele'] <- file.path(aux_dir, 'Mus_musculus.NCBIM37.67.gtf.gz')
dats$gtf_file [dats$group == 'lin'] <- file.path(dat_dir, 'lin', 'GSE60314_dme5_57_ERCC.gtf.gz')
dats$gtf_file [dats$group == 'nagarajan'] <- file.path(aux_dir, 'Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz')
dats$gtf_file [grepl('sarantopoulou', dats$group)] <- file.path(aux_dir, 'Mus_musculus.NCBIM37.67.gtf.gz')

write.table(dats, out_file, sep='\t', row.names=F, quote=F, col.names=T)
