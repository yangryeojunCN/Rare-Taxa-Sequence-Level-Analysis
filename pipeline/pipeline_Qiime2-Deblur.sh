cd /mnt/c/...
    
    qiime tools import \
      --input-path clean_data/manifest \
      --output-path Deblur_results/demux.qza \
      --type SampleData[JoinedSequencesWithQuality] \
      --input-format SingleEndFastqManifestPhred33

    qiime demux summarize \
      --i-data Deblur_results/demux.qza \
      --o-visualization Deblur_results/demux.qzv

    qiime quality-filter q-score \
      --i-demux Deblur_results/demux.qza \
      --o-filtered-sequences Deblur_results/demux-filtered.qza \
      --o-filter-stats Deblur_results/demux-filter-stats.qza
      
    qiime tools export \
      --input-path Deblur_results/demux-filter-stats.qza \
      --output-path Deblur_results/denoising-stats

    time qiime deblur denoise-16S \
      --i-demultiplexed-seqs Deblur_results/demux-filtered.qza \
      --p-trim-length 250 \
      --p-sample-stats \
      --o-representative-sequences Deblur_results/rep-seqs.qza \
      --o-table Deblur_results/table.qza \
      --o-stats Deblur_results/stats.qza

    qiime deblur visualize-stats \
      --i-deblur-stats Deblur_results/stats.qza \
      --o-visualization Deblur_results/stats.qzv

    qiime tools export \
      --input-path Deblur_results/stats.qzv \
      --output-path Deblur_results/denoising-stats

    qiime feature-table summarize \
      --i-table Deblur_results/table.qza \
      --o-visualization Deblur_results/table.qzv \
      --m-sample-metadata-file clean_data/sample-metadata.tsv
    
    qiime tools export \
      --input-path Deblur_results/table.qzv \
      --output-path Deblur_results/Feature_table_summary

    qiime feature-table tabulate-seqs \
      --i-data Deblur_results/rep-seqs.qza \
      --o-visualization Deblur_results/rep-seqs.qzv

    qiime feature-table rarefy \
      --i-table Deblur_results/table.qza \
      --p-sampling-depth 14446 \
      --o-rarefied-table Deblur_results/table-rarefied.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-338-806-nb-classifier-q2_202102.qza \
      --i-reads Deblur_results/rep-seqs.qza \
      --o-classification Deblur_results/taxonomy-silvaV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-515-806-nb-classifier-q2_202102.qza \
      --i-reads Deblur_results/rep-seqs.qza \
      --o-classification Deblur_results/taxonomy-silvaV4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier RDP-18-99-338-806-nb-classifier-q2_202102.qza \
      --i-reads Deblur_results/rep-seqs.qza \
      --o-classification Deblur_results/taxonomy-rdpV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier RDP-18-99-515-806-nb-classifier-q2_202102.qza \
      --i-reads Deblur_results/rep-seqs.qza \
      --o-classification Deblur_results/taxonomy-rdpV4.qza

      time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-338-806-nb-classifier-q2_202102.qza \
      --i-reads Deblur_results/rep-seqs.qza \
      --o-classification Deblur_results/taxonomy-MiDASV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-515-806-nb-classifier-q2_202102.qza \
      --i-reads Deblur_results/rep-seqs.qza \
      --o-classification Deblur_results/taxonomy-MiDASV4.qza

    qiime tools export \
      --input-path Deblur_results/table.qza \
      --output-path Deblur_results/R
    
    biom convert \
    -i Deblur_results/R/feature-table.biom \
    -o Deblur_results/R/feature_table.tsv \
    --to-tsv

    cd Deblur_results/R
    sed -i '1d' feature_table.tsv
    sed -i 's/#OTU ID/Feature ID/' feature_table.tsv

    qiime tools export \
      --input-path Deblur_results/table-rarefied.qza \
      --output-path Deblur_results/R

    biom convert \
      -i Deblur_results/R/feature-table.biom \
      -o Deblur_results/R/feature_table_rarefied.tsv \
      --to-tsv

    cd Deblur_results/R
    sed -i '1d' feature_table_rarefied.tsv
    sed -i 's/#OTU ID/Feature ID/' feature_table_rarefied.tsv
 
    qiime tools export \
      --input-path Deblur_results/rep-seqs.qza \
      --output-path Deblur_results/R

    qiime tools export \
      --input-path Deblur_results/taxonomy-silvaV3V4.qza \
      --output-path Deblur_results/R

    cd Deblur_results/R
    mv taxonomy.tsv taxonomy-silva-V3V4.tsv

    qiime tools export \
      --input-path Deblur_results/taxonomy-silvaV4.qza \
      --output-path Deblur_results/R

    cd Deblur_results/R
    mv taxonomy.tsv taxonomy-silva-V4.tsv

    qiime tools export \
      --input-path Deblur_results/taxonomy-MiDASV3V4.qza \
      --output-path Deblur_results/R

    cd Deblur_results/R
    mv taxonomy.tsv taxonomy-MiDAS-V3V4.tsv

    qiime tools export \
      --input-path Deblur_results/taxonomy-MiDASV4.qza \
      --output-path Deblur_results/R

    cd Deblur_results/R
    mv taxonomy.tsv taxonomy-MiDAS-V4.tsv
