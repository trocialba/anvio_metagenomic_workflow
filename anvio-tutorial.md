# Running anvio metagenomic workflow
In this tutorial it is provided an overview of the anvi’o workflow for the analysis of assembly-based shotgun metagenomic data as it is shown here: https://merenlab.org/2016/06/22/anvio-tutorial-v2/
In this tutorial we will use 8 test samples, located in 00_RAW folder

## Installation

The following softwares need to be installed: anvi'o, bowtie2, iu-filter-quality-minoche, samtools, megahit, anvi-init-bam, see here: https://anvio.org/install/

## Preparation
Assembly and mapping are key steps for most assembly-based, genome-resolved metagenomic studies, and there are many ways to accomplish each of these steps. That’s why, the anvi’o metagenomic workflow only starts once you have your contigs and BAM files available. One way how to get contigs and BAM files is shown here: https://merenlab.org/tutorials/assembly-based-metagenomics/

### 1. Quality filtering
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

### 2. Co-assembly

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
Replace $MIN_CONTIG_SIZE and $NUM_THREADS with actual values you want (Ex. 1000 for $MIN_CONTIG_SIZE, and 40 for $NUM_THREADS)

We have our contigs.fa under the directory 03_CONTIGS/.


### 3. Mapping

3.1. Create a new directory for mapping results and build an index for our contigs:

```
mkdir 04_MAPPING
bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
```

3.2. Get indexed BAM files for all samples:

```
for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    
    R1s=`ls 01_QC/$sample*QUALITY_PASSED_R1* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls 01_QC/$sample*QUALITY_PASSED_R2* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    
    bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 $R1s -2 $R2s --no-unal -S 04_MAPPING/$sample.sam
    samtools view -F 4 -bS 04_MAPPING/$sample.sam > 04_MAPPING/$sample-RAW.bam
    anvi-init-bam 04_MAPPING/$sample-RAW.bam -o 04_MAPPING/$sample.bam
    rm 04_MAPPING/$sample.sam 04_MAPPING/$sample-RAW.bam
done
```


