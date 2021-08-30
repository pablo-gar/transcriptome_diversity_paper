rm(list=ls())
library(SummarizedExperiment)
#library(Hmisc)
#library(ggplot2)
#library(pheatmap)
#library(RColorBrewer)

# folder where S3BUCKET data and github directory are stored. eg: ~/Downloads
bigdir = '/scratch/users/paedugar/cnv_gtex_project/expression_datasets/arora/'
# github directory eg: ~/Downloads/UncertaintyRNA
git_dir = file.path(bigdir,  "UncertaintyRNA")
# S3 bucket directory eg: ~/Downloads/OriginalTCGAGTExData
s3_dir = file.path(bigdir,  "OriginalTCGAGTExData")
# when you run our RMD files, all results will be stored here. 
# This will essentially remake the "data" subfolder from github repo.
# eg:~/Downloads/data
results_dir = file.path(bigdir, "data")
if(!file.exists( file.path(s3_dir, "SE_objects"))){
  stop("Please go through vignette 3 & 4 to make SE objects or download from S3 bucket")
}
if(!file.exists( file.path( results_dir))){
   system(paste0("mkdir ", results_dir))
}
if(!file.exists( file.path( results_dir, "pca_data"))){
   system(paste0("mkdir ", file.path(results_dir, "pca_data")))
}

maindir = file.path(results_dir, "pca_data")

tcga_gdc <- get(load( file.path( s3_dir, "SE_objects","tcga_gdc_log2_TPM.RData")))
tcga_mskcc_norm <- get(load( file.path( s3_dir, "SE_objects", 
                                        "tcga_mskcc_norm_log2_TPM.RData")))
tcga_mskcc_batch <- get(load( file.path( s3_dir, "SE_objects", 
                                         "tcga_mskcc_batch_log2_TPM.RData")))
tcga_recount2 <- get(load( file.path( s3_dir, "SE_objects", 
                                      "tcga_recount2_log2_TPM.RData")))
tcga_xena <- get(load( file.path( s3_dir, "SE_objects", 
                                  "tcga_xena_log2_TPM.RData")))
tcga_piccolo <- get(load( file.path( s3_dir, "SE_objects",
                                     "tcga_piccolo_log2_TPM.RData")))
gdc_mat = assay(tcga_gdc)
mskcc_norm_mat=assay(tcga_mskcc_norm)
mskcc_batch_mat=assay(tcga_mskcc_batch)
piccolo_mat=assay(tcga_piccolo)
recount2_mat=assay(tcga_recount2)
xena_mat= assay(tcga_xena)

final_all = cbind(gtex_v6_mat,  mskcc_norm_mat, mskcc_batch_mat,
                  recount2_mat, xena_mat)

gtex_tpm = prcomp(t(final_all))

saveRDS(gtex_tpm, file=file.path(results_dir, "pca_data", 'gtex_TPM.Rmd'))

# GTEx
gtex_v6 <- get(load( file.path( s3_dir, "SE_objects","gtex_v6_log2_TPM.RData")))
gtex_mskcc_norm <- get(load( file.path( s3_dir, "SE_objects",
                                        "gtex_mskcc_norm_log2_TPM.RData")))
gtex_mskcc_batch <- get(load( file.path( s3_dir, "SE_objects",
                                         "gtex_mskcc_batch_log2_TPM.RData")))
gtex_recount2 <- get(load( file.path( s3_dir, "SE_objects", 
                                      "gtex_recount2_log2_TPM.RData")))
gtex_xena <- get(load( file.path( s3_dir, "SE_objects",
                                  "gtex_xena_log2_TPM.RData")))
gtex_v6_mat = assay(gtex_v6)
mskcc_norm_mat=assay(gtex_mskcc_norm)
mskcc_batch_mat=assay(gtex_mskcc_batch)
recount2_mat=assay(gtex_recount2)
xena_mat= assay(gtex_xena)
