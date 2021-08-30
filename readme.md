This repository contains all the software and scripts used for the following study.

**Garcia-Nieto PE, Wang B, Fraser HB.** ***Transcriptome diversity is a systematic source of bias in RNA-sequencing data.*** 

The contents of this repository encompass:

- A snakeame pipeline to download all the data used in the study from the original sources.
- A snamake pipeline to process the data and reporoduce the analyses from the study.
    
***This document is not intended to be a thorough description of the methods or results. It is a guide that serves as a reproducibility reference for the code used in the aforementioned study.***

For a complete description of the methods and results please refer to the publication.



# Configuring pipelines

## Requirements
- Unix -- all software was run on bash in Ubuntu 14.04
- snakemake -- for all pipelines
- R 3.4:
	- broom
	- dplyr
	- edgeR
	- GenomicFeatures
	- ggplot2
	- peer
	- purr
	- readr
	- rjson
	- stringr
	- SummarizedExperiment
	- tidyr


## Setting up working path

Modify the file `config.json`, all buckets should be self-explanatory. The file is setup as it was used for the original study, which was run on the Sherlock cluster at Stanford University.

The following field indicates the root of where the data and analyses will be stored when running the pipelines. Please modfiy to your convenience:

```json
"projectDir": "/scratch/users/paedugar/transcriptome_diversity"
```

## Parallelizing pipelines

The pipelines are designed to run jobs in parallel. There are some extra files and scripts to run snakemake on a SLURM cluster:


- The file `cluster.json` contains rule-specific specifications and can be used on snakemake.
- `submit.py` is a custom submission script to be used for job submission on snakemake.
- `jobState` is a bash script to be used by snakemake to check for sate of jobs, **this script should be added to the bash `$PATH`**


# Downloading data from original sources


The following snakemake pipeline will download and processes the data from original sources and it will make it ready to be processed.

Linear execution:

```bash
cd pipelines
snakemake --snakefile Snakefile_download_datasets.smk --printshellcmds --keep-going --restart-times 2 
```

Parallel execution in a SLURM cluster:

```bash
cd pipelines
snakemake --snakefile Snakefile_download_datasets.smk --restart-times 2 --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```

-- 

**DISCLAIMER**: Some of the original sources are not available as of June 2021, The processed data has been deposited in google drive at:

[https://drive.google.com/file/d/1fSWt2ToSac5-usb8L0CpVzyw6eaIDgxe/view?usp=sharing](https://drive.google.com/file/d/1fSWt2ToSac5-usb8L0CpVzyw6eaIDgxe/view?usp=sharing)

# Running analysis from paper

The following snakemake pipeline will execute the analyses presented in the paper. The data has to be first downloaded (see above).

Linear execution:

```bash
cd pipelines
snakemake --snakefile Snakefile_main.smk --printshellcmds --keep-going --restart-times 2 
```

Parallel execution in a SLURM cluster:

```bash
cd pipelines
snakemake --snakefile Snakefile_main.smk --restart-times 2 --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```
 
