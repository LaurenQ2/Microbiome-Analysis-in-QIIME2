## This tutorial presents a QIIME2 data analysis pipeline that can be used with 16S (bacterial) and ITS (fungal) sequences from a soil sampling study.

### Commands with separate 16S and ITS pipelines are: DADA2 and taxonomic classification. 
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
  
#### Reference code provided in this ReadMe.md is a truncated version of the full code; [that file can be found](https://github.com/LaurenQ2/Microbiome-Analysis-in-QIIME2/blob/main/microbiome-analysis-w-qiime2.sh)  in the main directory.  
#### Note on the reference code: 16S is used as an example in all commands where the same code is used for both 16S and ITS.
&nbsp;
#### First, let's set the stage. Load the QIIME2 package.
``` bash
#!/bin/bash
module load qiime2/2024.5
```
#### Sequence data was received from the sequencing lab trimmed, demultiplexed, and zipped. In these first couple steps we'll unzip, import the files, and run a quality check

``` bash
unzip -d $PWD/raw_data/fastq_trimmed.zip
gunzip $PWD/raw_data/unZipped/filepath_R1.fastq.gz > $PWD/fastQ/filepath_R1.fq
```
> In order to import these sequencing files we need to know whether they have a Phred quality score of 33 or 64. This code will check the Phred quality score:
``` bash
vsearch --fastq_chars $PWD/fastQ/filepath_R1.fastq
```
> Link forward and reverse read paths with other sample metadata in a metadata.tsv file, then import the files.
``` bash
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $PWD/fastQ/metadata_16S.tsv \
--output-path $PWD/fastQ/16S_paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2
```
> Summarize the imported file stats. \
> These .qzv visualization files can be viewed in the [QIIME2 visualizer.](https://view.qiime2.org)
``` bash
 qiime demux summarize \
--i-data $PWD/fastQ/16S_paired-end-demux.qza \
--o-visualization $PWD/fastQ/16S_paired-end-demux.qzv
```
&nbsp;
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
&nbsp;
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
* **Steps:** 
  + Sequence alignment with MAFFT,   
  + filtering with mask 
  + Phylogeny built with FastTree,
  + Root the tree, download the .qzv file, and generate a phylogeny visualization by going to the website <https://itol.embl.de/upload.cgi>.

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
&nbsp;
#### Rarefaction curves
* Set max-depth to the maximum number of sequences in your least populated sample.
* Remove metadata file in second command below to separate individual sample curves in [the QIIME2 visualizer.](https://view.qiime2.org/)
``` bash
qiime diversity alpha-rarefaction \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--p-max-depth 47385 \
--p-steps 20 \
--i-phylogeny $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--o-visualization $PWD/analysis/16S/rarefaction/16S_rarefaction_curves.qzv

qiime diversity alpha-rarefaction \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--p-max-depth 47385 \
--p-steps 20 \
--i-phylogeny $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--o-visualization $PWD/analysis/16S/rarefaction/16S_rarefaction_curves_each_curve.qzv
```
&nbsp;
#### Taxonomy 16S 
* Classifier used was Greengenes2
* Reference databases downloaded at <http://ftp.microbio.me/greengenes_release/2022.10/>
* Specific reference databases used: `2022.10.backbone.full-length.fna.qza` and `2022.10.backbone.tax.qza`
> First, install the Greengenes2 program.
``` bash
module load qiime2/2024.5
pip install q2-greengenes2
```
> Next, train the classifier.
``` bash
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads $PWD/analysis/2022.10.backbone.full-length.fna.qza \
--i-reference-taxonomy $PWD/analysis/2022.10.backbone.tax.qza \
--o-classifier $PWD/analysis/2022.10.backbone.full-length.nb.qza
```
> Then, run the classifier on the input reads from your study.
``` bash
 qiime feature-classifier classify-sklearn \
--i-classifier $PWD/analysis/2022.10.backbone.full-length.nb.qza \
--i-reads $PWD/analysis/16S/filtered/16S_filtered-rep-seqs.qza \
--o-classification $PWD/analysis/16S/taxonomy/16S_taxonomy.qza
```
> Export the output file to look at the classifications and confidence scores.
> Output file will be named 'taxonomy.tsv'.
``` bash
qiime tools export \
--input-path $PWD/analysis/16S/taxonomy/16S_taxonomy.qza \
--output-path $PWD/analysis/16S/taxonomy
```
* Create a barplot visualization of your taxonomy data.
* This barplot, and the exportable .csv file, is helpful for identifying differences in  taxonomic distribution between samples.
``` bash
 qiime taxa barplot \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--i-taxonomy $PWD/analysis/16S/taxonomy/16S_taxonomy.qza \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--o-visualization $PWD/analysis/16S/taxonomy/16S_taxa-barplot.qzv
```
&nbsp;
#### Taxonomy ITS 
* Classifier used was UNITE
* Reference databases downloaded at <https://doi.org/10.15156/BIO/2959339>
* Specific pre-trained reference database used: `unite_ver10_dynamic_all_04.04.2024-Q2-2024.5.qza`
 

