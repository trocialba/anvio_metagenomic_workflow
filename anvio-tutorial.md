# Running anvio metagenomic workflow
In this tutorial it is provided an overview of the anvi’o workflow for the analysis of assembly-based shotgun metagenomic data as it is shown here: https://merenlab.org/2016/06/22/anvio-tutorial-v2/
In this tutorial we will use 8 test samples, located in 00_RAW folder

## Installation

The following softwares need to be installed: anvi'o, bowtie2, iu-filter-quality-minoche, samtools, megahit, anvi-init-bam, see here: https://anvio.org/install/


## 1. Quality filtering
1.1. First, create a TAB-delimited table samples.txt where are our raw R1 and R2 files (see samples.txt file)

1.2. Create a directory for Quality Filtered R1, R2 and then use iu-gen-configs to create config files for ilumina-utils in it:

```
mkdir 01_QC
iu-gen-configs samples.txt -o 01_QC
ls 01_QC/
```

1.3. Quality filtering for each sample:

```
iu-filter-quality-minoche 01_QC/Sample_01.ini
```

or for all samples at once:

```
for ini in 01_QC/*.ini; do iu-filter-quality-minoche $ini; done
```


1.4. QC check output

```
grep 'total pairs passed' 01_QC/*STATS.txt
```

## 2. Co-assembly

2.1. Create two environment variables:

```
R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

echo $R1s

echo $R2s
```


2.2. Run MEGAHIT:

```
megahit -1 $R1s -2 $R2s --min-contig-len $MIN_CONTIG_SIZE -m 0.85 -o 02_ASSEMBLY/ -t $NUM_THREADS

mkdir 03_CONTIGS
```

