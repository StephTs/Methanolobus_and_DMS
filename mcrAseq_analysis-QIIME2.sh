##### QIIME2 - Analysis of Baltic _mcrA_ amplicon sequences #####
# Author: Stephania L. Tsola
# Date: May 2022
# QIIME2 verion 2021.11


conda activate qiime2-2021.11

# 1) Import files to QIIME2
# Makes one single qiime file with all the sequences
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path home/mcra \
--source-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path mcrA_demux-paired-end.qza

# 2) Visualise merged file (the file created helps to check quality control and understand where to trim and truncate your samples in the DADA2 step (step 3)
qiime demux summarize \
--i-data mcrA_demux-paired-end.qza \
--o-visualization mcrA_demux.qzv

# 3) Sequence quality control
qiime dada2 denoise-paired \
--i-demultiplexed-seqs mcrA_demux-paired-end.qza \
--o-table mcrA_table.qza \
--o-representative-sequences mcrA_rep-seqs.qza \
--o-denoising-stats mcrA_den-stats.qza \
--p-trim-left-f 18 \
--p-trim-left-r 9 \
--p-trunc-len-f 301 \
--p-trunc-len-r 271 \
--verbose

# 4) Visualize the denoising stats
qiime metadata tabulate \
--m-input-file mcrA_den-stats.qza \
--o-visualization mcrA_den-stats.qzv

# 5) Feature table and feature data visualization
qiime feature-table summarize \
--i-table mcrA_table.qza \
--o-visualization mcrA_table.qzv \
--m-sample-metadata-file mcrA-metadata.tsv

qiime feature-table mcrA_tabulate-seqs \
--i-data mcrA_rep-seqs.qza \
--o-visualization mcrA_rep-seqs.qzv

# 8) Train a classifier using the mcrA database (fasta and taxonomy files)
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  NoBacUncArc5Env-seqs-cleaned.qza \
  --i-reference-taxonomy NoBacUncArc5Env-tax.qza \
  --o-classifier mcrA-classifier.qza

# 9) Use the classifier and visualize the resulting taxonomic assignments
qiime feature-classifier classify-sklearn \
--i-classifier mcrA-classifier.qza \
--i-reads mcrA_rep-seqs.qza \
--o-classification mcrA_taxonomy.qza \
  
qiime metadata tabulate \
--m-input-file mcrA_taxonomy.qza \
--o-visualization mcrA_taxonomy.qzv

# 10) View the taxonomic composition of the samples with interactive bar plots
qiime taxa barplot \
--i-table mcrA_table.qza \
--i-taxonomy mcrA_taxonomy.qza \
--m-metadata-file mcrA-metadata.tsv \
--o-visualization mcrA_taxa-bar-plots.qzv
 
