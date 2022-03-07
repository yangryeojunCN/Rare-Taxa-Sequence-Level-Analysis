    cd /mnt/c/...

    qiime tools import \
       --input-path clean_data/manifest \
       --output-path DaDa2_results/demux.qza \
       --type SampleData[PairedEndSequencesWithQuality] \
       --input-format PairedEndFastqManifestPhred33

    cd /mnt/c/...
    qiime demux summarize \
      --i-data DaDa2_results/demux.qza \
      --o-visualization DaDa2_results/demux.qzv

    time qiime dada2 denoise-paired \
      --i-demultiplexed-seqs DaDa2_results/demux.qza \
      --p-n-threads 0 \
      --p-trunc-len-f 270 \
      --p-trunc-len-r 210 \
      --p-max-ee-f 2 \
      --p-max-ee-r 2 \
      --o-table DaDa2_results/table.qza \
      --o-representative-sequences DaDa2_results/rep-seqs.qza \
      --o-denoising-stats DaDa2_results/denoising-stats.qza

    qiime metadata tabulate \
      --m-input-file DaDa2_results/denoising-stats.qza \
      --o-visualization DaDa2_results/denoising-stats.qzv

    qiime tools export \
       --input-path DaDa2_results/denoising-stats.qzv \
       --output-path DaDa2_results/denoising-stats

     qiime feature-table summarize \
       --i-table DaDa2_results/table.qza \
       --o-visualization DaDa2_results/table.qzv \
       --m-sample-metadata-file clean_data/sample-metadata.tsv
     qiime tools export \
       --input-path DaDa2_results/table.qzv \
       --output-path DaDa2_results/Feature_table_summary
      
     qiime feature-table tabulate-seqs \
       --i-data DaDa2_results/rep-seqs.qza \
       --o-visualization DaDa2_results/rep-seqs.qzv

    qiime feature-table rarefy \
      --i-table DaDa2_results/table.qza \
      --p-sampling-depth 22627 \
      --o-rarefied-table DaDa2_results/table-rarefied.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-338-806-nb-classifier-q2_202102.qza \
      --i-reads DaDa2_results/rep-seqs.qza \
      --o-classification DaDa2_results/taxonomy-silvaV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-515-806-nb-classifier-q2_202102.qza \
      --i-reads DaDa2_results/rep-seqs.qza \
      --o-classification DaDa2_results/taxonomy-silvaV4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-338-806-nb-classifier-q2_202102.qza \
      --i-reads DaDa2_results/rep-seqs.qza \
      --o-classification DaDa2_results/taxonomy-MiDASV3V4.qza

       time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-515-806-nb-classifier-q2_202102.qza \
      --i-reads DaDa2_results/rep-seqs.qza \
      --o-classification DaDa2_results/taxonomy-MiDASV4.qza

    qiime tools export \
      --input-path DaDa2_results/table.qza \
      --output-path DaDa2_results/R
    
    biom convert \
    -i DaDa2_results/R/feature-table.biom \
    -o DaDa2_results/R/feature_table.tsv \
    --to-tsv

    cd DaDa2_results/R
    sed -i '1d' feature_table.tsv
    sed -i 's/#OTU ID/Feature ID/' feature_table.tsv

    qiime tools export \
      --input-path DaDa2_results/table-rarefied.qza \
      --output-path DaDa2_results/R

    biom convert \
      -i DaDa2_results/R/feature-table.biom \
      -o DaDa2_results/R/feature_table_rarefied.tsv \
      --to-tsv

    cd DaDa2_results/R
    sed -i '1d' feature_table_rarefied.tsv
    sed -i 's/#OTU ID/Feature ID/' feature_table_rarefied.tsv
    
    qiime tools export \
      --input-path DaDa2_results/rep-seqs.qza \
      --output-path DaDa2_results/R
  
    qiime tools export \
      --input-path DaDa2_results/taxonomy-silvaV3V4.qza \
      --output-path DaDa2_results/R

    cd DaDa2_results/R
    mv taxonomy.tsv taxonomy-silva-V3V4.tsv

    qiime tools export \
      --input-path DaDa2_results/taxonomy-silvaV4.qza \
      --output-path DaDa2_results/R

    cd DaDa2_results/R
    mv taxonomy.tsv taxonomy-silva-V4.tsv

    cd /mnt/c/Omics/qiime2/Rare_taxa

    qiime tools export \
      --input-path DaDa2_results/taxonomy-MiDASV3V4.qza \
      --output-path DaDa2_results/R

    cd DaDa2_results/R
    mv taxonomy.tsv taxonomy-MiDAS-V3V4.tsv

    qiime tools export \
      --input-path DaDa2_results/taxonomy-MiDASV4.qza \
      --output-path DaDa2_results/R

    cd DaDa2_results/R
    mv taxonomy.tsv taxonomy-MiDAS-V4.tsv
