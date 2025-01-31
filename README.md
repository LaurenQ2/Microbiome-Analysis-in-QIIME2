## This tutorial presents a QIIME2 data analysis pipeline with 16S (bacterial) and ITS (fungal) sequences from a soil sampling study.

* **Featured commands:** 
  + Importing, 
  + Filtering with DADA2, 
  + Removing rare ASVs, 
  + Phylogeny with MAFFT and FastTree, 
  + 16S taxonomy with Greengenes2,
  + ITS taxonomy with UNITE,
  + Taxonomy barplots
  + Core diversity metrics

* Additional diversity visualizations were produced using R and the `QIIME2R` suite of packages and tools. Additional .R code will be linked in a future update.

#### First, let's set the stage. Load the QIIME2 package.
``` bash
#!/bin/bash
module load qiime2/2024.5
```
&nbsp;
> Sequence data was received from the sequencing lab trimmed, demultiplexed, and zipped. In these first couple steps we'll unzip, import the files, and run a quality check

``` bash
unzip -d $PWD/raw_data/fastq_trimmed.zip
gunzip $PWD/raw_data/unZipped/filepath_R1.fastq.gz > $PWD/fastQ/filepath_R1.fq
```
&nbsp;
> In order to import these sequencing files we need to know whether they have a Phred quality score of 33 or 64. This code will check the Phred quality score for you:
``` bash
vsearch --fastq_chars $PWD/fastQ/filepath_R1.fastq
```
