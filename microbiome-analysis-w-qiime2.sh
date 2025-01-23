#!/bin/bash

#unzip data
unzip -d $PWD/raw_data/fastq_trimmed.zip

#move unzipped .gz files to separate directory called "unZipped"
#create fastQ directory in parent directory
#run command below for all .gz files
gunzip $PWD/raw_data/unZipped/filepath_R1.fastq.gz > $PWD/fastQ/filepath_R1.fq

#determine Phred quality
module load qiime2/2024.5
vsearch --fastq_chars $PWD/fastQ/filepath_R1.fastq

#create separate .tsv metadata files for 16S and ITS
#create the following columns in row 1: sample-id, forward-absolute-filepath, reverse-absolute-filepath
#additional columns in row 1 can be added as required by the study
#identify data types in row 2, for first three columns these are: #q2:types, categorical, categorical (other data types may be 'numerical')
#this metadata file connects forward (R1) and reverse (R2) reads to other critical metadata

#import 16S samples into qiime2
#Metadata file functions as the manifest 'input-path' in this command
#Phred quality score used in input-format will be either 'Phred33V2' or 'Phred64V2'
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $PWD/fastQ/metadata_16S.tsv \
--output-path $PWD/fastQ/16S_paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

#repeat command for ITS files
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $PWD/fastQ/metadata_ITS.tsv \
--output-path $PWD/fastQ/ITS_paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

#summarize sequence data
#create .qzv visualization file
qiime demux summarize \
--i-data $PWD/fastQ/16S_paired-end-demux.qza \
--o-visualization $PWD/fastQ/16S_paired-end-demux.qzv

qiime demux summarize \
--i-data $PWD/fastQ/ITS_paired-end-demux.qza \
--o-visualization $PWD/fastQ/ITS_paired-end-demux.qzv

#download .qzv files
#open .qzv files in the QIIME2 visualizer (https://view.qiime2.org/)
#review the 'Interactive Quality Plot' tab to identify low quality regions within the sequences for removal

#DADA2 16S
#create a new directory called 'Analysis' in the parent directory
module load qiime2/2024.5

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

#visualize denoising stats by creating a .qzv
#open .qzv files in the QIIME2 visualizer (https://view.qiime2.org/)
qiime metadata tabulate \
--m-input-file $PWD/analysis/16S/16S_denoising-stats.qza \
--o-visualization $PWD/analysis/16S/16S_denoising-stats.qzv

#DADA2 summary
qiime feature-table summarize \
--i-table $PWD/analysis/16S/16S_table.qza \
--o-visualization $PWD/analysis/16S/16S_table-summary.qzv

#DADA2 ITS
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

#visualize denoising stats by creating a .qzv
qiime metadata tabulate \
--m-input-file $PWD/ITS/ITS_denoising-stats.qza \
--o-visualization $PWD/analysis/ITS/ITS_denoising-stats.qzv

#DADA2 summary
qiime feature-table summarize \
--i-table $PWD/analysis/ITS/ITS_table.qza \
--o-visualization $PWD/analysis/ITS/ITS_table-summary.qzv

#filter out rare ASVs
#create new subdirectory in 16S and ITS folders called 'filtered'
module load qiime2/2024.5

qiime feature-table filter-features \
--i-table $PWD/analysis/16S/16S_table.qza \
--p-min-frequency 5 \
--o-filtered-table $PWD/analysis/16S/filtered/16S_filtered-table.qza

qiime feature-table filter-features \
--i-table $PWD/analysis/ITS/ITS_table.qza \
--p-min-frequency 7 \
--o-filtered-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza

#filtered rep-seqs
qiime feature-table filter-seqs \
--i-data $PWD/analysis/16S/16S_rep-seqs.qza \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--o-filtered-data $PWD/analysis/16S/filtered/16S_filtered-rep-seqs.qza

qiime feature-table filter-seqs \
--i-data $PWD/analysis/ITS/ITS_rep-seqs.qza \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza \
--o-filtered-data $PWD/analysis/ITS/filtered/ITS_filtered-rep-seqs.qza

#Create a visualization of the Feature Table
qiime feature-table summarize \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--o-visualization $PWD/analysis/16S/filtered/16S_filtered-table-summary.qzv

qiime feature-table summarize \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza \
--o-visualization $PWD/analysis/ITS/filtered/ITS_filtered-table-summary.qzv

#phylogeny
#steps: MAFFT, filtering, FastTree, rooted tree
#create new subdirectory in 16S and ITS folders called 'phylogeny'
#Use rooted tree .qzv to generate a phylogeny visualization, 
#by going to the website, https://itol.embl.de/upload.cgi
module load qiime2/2024.5

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

#repeat for ITS files
qiime alignment mafft \
--i-sequences $PWD/analysis/ITS/filtered/ITS_filtered-rep-seqs.qza \
--o-alignment $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned.qza

qiime alignment mask \
--i-alignment $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned.qza \
--o-masked-alignment $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked.qza

qiime phylogeny fasttree \
--i-alignment $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked.qza \
--o-tree $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked_tree.qza

qiime phylogeny midpoint-root \
--i-tree $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked_tree.qza \
--o-rooted-tree $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked_tree_rooted.qza

#rarefaction curves
#set max-depth to the max number of sequences in the least populated sample
#create new subdirectory in 16S and ITS folders called 'rarefaction'
module load qiime2/2024.5

qiime diversity alpha-rarefaction \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--p-max-depth 47385 \
--p-steps 20 \
--i-phylogeny $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--o-visualization $PWD/analysis/16S/rarefaction/16S_rarefaction_curves.qzv

qiime diversity alpha-rarefaction \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza \
--p-max-depth 89899 \
--p-steps 20 \
--i-phylogeny $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--m-metadata-file $PWD/fastQ/metadata_ITS.tsv \
--o-visualization $PWD/analysis/ITS/rarefaction/ITS_rarefaction_curves.qzv

#remove metadata file to separate individual sample curves in the .qzv visualizer
qiime diversity alpha-rarefaction \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--p-max-depth 47385 \
--p-steps 20 \
--i-phylogeny $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--o-visualization $PWD/analysis/16S/rarefaction/16S_rarefaction_curves_each_curve.qzv

qiime diversity alpha-rarefaction \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza \
--p-max-depth 89899 \
--p-steps 20 \
--i-phylogeny $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--o-visualization $PWD/analysis/ITS/rarefaction/ITS_rarefaction_curves.qzv

#taxonomy
#16S classifier used was Greengenes2
#reference databases downloaded at http://ftp.microbio.me/greengenes_release/2022.10/
#reference databases uploaded to $PWD/analysis directory:
#2022.10.backbone.full-length.fna.qza and 2022.10.backbone.tax.qza
#install Greengenes2
module load qiime2/2024.5
pip install q2-greengenes2

#train classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads $PWD/analysis/2022.10.backbone.full-length.fna.qza \
--i-reference-taxonomy $PWD/analysis/2022.10.backbone.tax.qza \
--o-classifier $PWD/analysis/2022.10.backbone.full-length.nb.qza

#run classifier
#create new subdirectory in $PWD/analysis/16S called 'taxonomy' 
qiime feature-classifier classify-sklearn \
--i-classifier $PWD/analysis/2022.10.backbone.full-length.nb.qza \
--i-reads $PWD/analysis/16S/filtered/16S_filtered-rep-seqs.qza \
--o-classification $PWD/analysis/16S/taxonomy/16S_taxonomy.qza

#export the output file to look at the classifications and confidence scores
#output file will be named 'taxonomy.tsv'
qiime tools export \
--input-path $PWD/analysis/16S/taxonomy/16S_taxonomy.qza \
--output-path $PWD/analysis/16S/taxonomy

#barplot visualization
qiime taxa barplot \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--i-taxonomy $PWD/analysis/16S/taxonomy/16S_taxonomy.qza \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--o-visualization $PWD/analysis/16S/taxonomy/16S_taxa-barplot.qzv

#ITS pre-trained classifier used was UNITE
#reference reads from UNITE downloaded at https://doi.org/10.15156/BIO/2959339
#reference database uploaded to $PWD/analysis directory:
#unite_ver10_dynamic_all_04.04.2024-Q2-2024.5.qza 
#create new subdirectory in $PWD/analysis/ITS called 'taxonomy'
qiime feature-classifier classify-sklearn \
--i-classifier $PWD/analysis/unite_ver10_dynamic_all_04.04.2024-Q2-2024.5.qza \
--i-reads $PWD/analysis/ITS/filtered/ITS_filtered-rep-seqs.qza \
--o-classification $PWD/analysis/ITS/taxonomy/ITS_taxonomy.qza 

#export a 'taxonomy.tsv' file to look at classifications and confidence scores
qiime tools export \
--input-path $PWD/analysis/ITS/taxonomy/ITS_taxonomy.qza \
--output-path $PWD/analysis/ITS/taxonomy

#barplot visualization
qiime taxa barplot \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza \
--i-taxonomy $PWD/analysis/ITS/taxonomy/ITS_taxonomy.qza \
--m-metadata-file $PWD/fastQ/metadata_ITS.tsv \
--o-visualization $PWD/analysis/ITS/taxonomy/ITS_taxa-barplot.qzv 

#taxonomy based filtering
#remove mitochondria and chloroplasts
module load qiime2/2024.5
qiime taxa filter-table \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table.qza \
--i-taxonomy $PWD/analysis/16S/taxonomy/16S_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table $PWD/analysis/16S/filtered/16S_filtered-table-no-contam.qza

qiime taxa filter-seqs \
--i-sequences $PWD/analysis/16S/filtered/16S_filtered-rep-seqs.qza \
--i-taxonomy $PWD/analysis/16S/taxonomy/16S_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences $PWD/analysis/16S/filtered/16S_filtered-rep-seqs-no-contam.qza

#check sequence counts 
qiime feature-table summarize \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table-no-contam.qza \
--o-visualization $PWD/analysis/16S/filtered/16S_filtered-table-no-contam-summary.qzv

#repeat commands for ITS files
qiime taxa filter-table \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table.qza \
--i-taxonomy $PWD/analysis/ITS/taxonomy/ITS_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table $PWD/analysis/ITS/filtered/ITS_filtered-table-no-contam.qza

qiime taxa filter-seqs \
--i-sequences $PWD/analysis/ITS/filtered/ITS_filtered-rep-seqs.qza \
--i-taxonomy $PWD/analysis/ITS/taxonomy/ITS_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences $PWD/analysis/ITS/filtered/ITS_filtered-rep-seqs-no-contam.qza

qiime feature-table summarize \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table-no-contam.qza \
--o-visualization $PWD/analysis/ITS/filtered/ITS_filtered-table-no-contam-summary.qzv

#diversity
#use max-depth from rarefaction step for p-sampling depth
#create new subdirectories in $PWD/analysis/16S and $PWD/analysis/ITS called 'diversity'
module load qiime2/2024.5

qiime diversity core-metrics-phylogenetic \
--i-table $PWD/analysis/16S/filtered/16S_filtered-table-no-contam.qza \
--i-phylogeny $PWD/analysis/16S/phylogeny/16S_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--p-sampling-depth 47385 \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--output-dir $PWD/analysis/16S/diversity

qiime diversity core-metrics-phylogenetic \
--i-table $PWD/analysis/ITS/filtered/ITS_filtered-table-no-contam.qza \
--i-phylogeny $PWD/analysis/ITS/phylogeny/ITS_filtered-rep-seqs-aligned_masked_tree_rooted.qza \
--p-sampling-depth 89899 \
--m-metadata-file $PWD/fastQ/metadata_ITS.tsv \
--output-dir $PWD/analysis/ITS/diversity

#export shannon vector values
#this command will create a subdirectory within $PWD/analysis/16S/diversity 
#named "shannon-diversity" containing a file of per sample shannon values  
#named 'alpha-diversity.tsv' 
qiime tools export \
--input-path $PWD/analysis/16S/diversity/shannon_vector.qza \
--output-path $PWD/analysis/16S/diversity/shannon-diversity

qiime tools export \
--input-path $PWD/analysis/ITS/diversity/shannon_vector.qza \
--output-path $PWD/analysis/ITS/diversity/shannon-diversity

#export faith pd vector values
#this command will create a subdirectory within $PWD/analysis/16S/diversity 
#named "faith-pd" containing a file of per sample Faith PD values  
#named 'alpha-diversity.tsv'
qiime tools export \
--input-path $PWD/analysis/16S/diversity/faith_pd_vector.qza \
--output-path $PWD/analysis/16S/diversity/faith-pd

qiime tools export \
--input-path $PWD/analysis/ITS/diversity/faith_pd_vector.qza \
--output-path $PWD/analysis/ITS/diversity/faith-pd

#export evenness vector values
#this command will create a subdirectory within $PWD/analysis/16S/diversity 
#named "evenness" containing a file of per sample evenness values  
#named 'alpha-diversity.tsv'
qiime tools export \
--input-path $PWD/analysis/16S/diversity/evenness_vector.qza \
--output-path $PWD/analysis/16S/diversity/evenness

qiime tools export \
--input-path $PWD/analysis/ITS/diversity/evenness_vector.qza \
--output-path $PWD/analysis/ITS/diversity/evenness

#create boxplots comparing metadata categories
#for this tutorial we will compare “rootstock” and “scion” categories 
module load qiime2/2024.5

#Shannon diversity boxplots
qiime diversity alpha-group-significance \
--i-alpha-diversity $PWD/analysis/16S/diversity/shannon_vector.qza \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--o-visualization $PWD/analysis/16S/diversity/shannon-diversity/16S_shannon-compare-groups.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity $PWD/analysis/ITS/diversity/shannon_vector.qza \
--m-metadata-file $PWD/fastQ/metadata_ITS.tsv \
--o-visualization $PWD/analysis/ITS/diversity/shannon-diversity/ITS_shannon-compare-groups.qzv

#Faith PD boxplots
qiime diversity alpha-group-significance \
--i-alpha-diversity $PWD/analysis/16S/diversity/faith_pd_vector.qza \
--m-metadata-file $PWD/fastQ/metadata_16S.tsv \
--o-visualization $PWD/analysis/16S/diversity/faith-pd/16S_faith-pd-compare-groups.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity $PWD/analysis/ITS/diversity/faith_pd_vector.qza \
--m-metadata-file $PWD/fastQ/metadata_ITS.tsv \
--o-visualization $PWD/analysis/ITS/diversity/faith-pd/ITS_faith-pd-compare-groups.qzv

#Additional visualizations, including Beta diversity PCoA plots, produced in RStudio