## This tutorial presents a QIIME2 data analysis pipeline that can be used with 16S (bacterial) and ITS (fungal) sequences from a soil sampling study.
### 16S is used as an example in all commands where the same code is used for both 16S and ITS
&nbsp;
#### Commands with separate 16S and ITS pipelines are DADA2 and taxonomic classification 
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
#### Sequence data was received from the sequencing lab trimmed, demultiplexed, and zipped. In these first couple steps we'll unzip, import the files, and run a quality check

``` bash
unzip -d $PWD/raw_data/fastq_trimmed.zip
gunzip $PWD/raw_data/unZipped/filepath_R1.fastq.gz > $PWD/fastQ/filepath_R1.fq
```
&nbsp;
> In order to import these sequencing files we need to know whether they have a Phred quality score of 33 or 64. This code will check the Phred quality score:
``` bash
vsearch --fastq_chars $PWD/fastQ/filepath_R1.fastq
```
> Import the files.
``` bash
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $PWD/fastQ/metadata_16S.tsv \
--output-path $PWD/fastQ/16S_paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2
```
> Summarize the imported file stats. \
> These .qzv visualization files can be viewed in the [QIIME2 viewer.](https:/view.qiime2.org)
``` bash
 qiime demux summarize \
--i-data $PWD/fastQ/16S_paired-end-demux.qza \
--o-visualization $PWD/fastQ/16S_paired-end-demux.qzv
```
#### DADA2 16S Bacterial Filtering
``` bash
qiime dada2 denoise-paired \
--i-demultiplexed-seqs $PWD/fastQ/16S_paired-end-demux.qza \
--p-trim-left-f 5 \
--p-trim-left-r 5 \
--p-trunc-len-f 200 \
--p-trunc-len-r 100 \
--p-n-threads 2 \
--o-table $PWD/analysis/16S/16S_table.qza \
--o-representative-sequences $PWD/analysis/16S/16S_rep-seqs.qza \
--o-denoising-stats $PWD/analysis/16S/16S_denoising-stats.qza
```
#### DADA2 ITS Fungal Filtering
``` bash
qiime dada2 denoise-paired \
--i-demultiplexed-seqs $PWD/fastQ/ITS_paired-end-demux.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-max-ee-f 3 \
--p-max-ee-r 5 \
--p-trunc-q 2 \
--p-pooling-method pseudo \
--p-n-threads 2 \
--o-table $PWD/analysis/ITS/ITS_table.qza \
--o-representative-sequences $PWD/analysis/ITS/ITS_rep-seqs.qza \
--o-denoising-stats $PWD/analysis/ITS/ITS_denoising-stats.qza
```
&nbsp;
> Create data visualizations from the output files generated in the last step.
``` bash
qiime metadata tabulate \
--m-input-file $PWD/analysis/16S/16S_denoising-stats.qza \
--o-visualization $PWD/analysis/16S/16S_denoising-stats.qzv

qiime feature-table summarize \
--i-table $PWD/analysis/16S/16S_table.qza \
--o-visualization $PWD/analysis/16S/16S_table-summary.qzv
```
&nbsp;
#### Filter out rare ASVs
``` bash
qiime feature-table filter-features \
--i-table $PWD/analysis/16S/16S_table.qza \
--p-min-frequency 5 \
--o-filtered-table $PWD/analysis/16S/filtered/16S_filtered-table.qza

qiime feature-table filter-seqs \
--i-data $PWD/analysis/16S/16S_rep-seqs.qza \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--o-filtered-data $PWD/analysis/16S/filtered/16S_filtered-rep-seqs.qza

qiime feature-table summarize \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--o-visualization $PWD/analysis/16S/filtered/16S_filtered-table-summary.qzv
```
&nbsp;
#### Phylogeny
``` bash
qiime alignment mafft \
--i-sequences $PWD/analysis/16S/filtered/16S_filtered-rep-seqs.qza \
--o-alignment $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned.qza

qiime alignment mask \
--i-alignment $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned.qza \
--o-masked-alignment $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked.qza

qiime phylogeny fasttree \
--i-alignment $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked.qza \
--o-tree $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree.qza

qiime phylogeny midpoint-root \
--i-tree $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree.qza \
--o-rooted-tree $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree_rooted.qza
```
> Repeat commands above with ITS files.
&nbsp;
#### Taxonomy

&nbsp;
