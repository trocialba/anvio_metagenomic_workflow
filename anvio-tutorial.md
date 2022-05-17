# Running anvio metagenomic workflow
In this tutorial it is provided an overview of the anvi’o workflow for the analysis of assembly-based shotgun metagenomic data as it is shown here: https://merenlab.org/2016/06/22/anvio-tutorial-v2/
In this tutorial we will use 8 test samples, located in 00_RAW folder

## Installation

The following softwares need to be installed: anvi'o, bowtie2, iu-filter-quality-minoche, samtools, megahit, anvi-init-bam, see here: https://anvio.org/install/


## 1. Quality filtering
1.1. First, create a TAB-delimited table samples.txt where are our raw R1 and R2 files (see samples.txt file)
1.2. Create a directory for Quality Filtered R1, R2 and then use iu-gen-configs to create config files for ilumina-utils in it:

mkdir 01_QC
iu-gen-configs samples.txt -o 01_QC
ls 01_QC/

`mkdir 01_QC`
