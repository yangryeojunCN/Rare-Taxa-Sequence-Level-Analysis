    cd /c/...
    cat USEARCH_results/*.merged.fq > USEARCH_results/all.merged.fq
    
    usearch -fastq_filter all.merged.fq -fastq_maxee 1.0 \
      -fastaout filtered.fa -relabel Filt

    usearch -fastx_uniques filtered.fa \
      -sizeout -relabel Uniq \
      -fastaout uniques.fa

    usearch -cluster_otus uniques.fa  \
      -otus otus.fa -uparseout uparse.txt -relabel Otu
      
    usearch -unoise3 uniques.fa  \
    -zotus zotus.fa -tabbedout unoise3.txt

    usearch -otutab all.merged.fq -otus otus.fa -otutabout otutab_raw.txt
    usearch -otutab all.merged.fq -zotus zotus.fa -otutabout zotutab_raw.txt
  
    awk '!/^>/ { printf "%s", $0; n = "\n" } 
    /^>/ { print n $0; n = "" }
    END { printf "%s", n }
    ' otus.fa > otus_Qiime2.fa

    awk '!/^>/ { printf "%s", $0; n = "\n" } 
    /^>/ { print n $0; n = "" }
    END { printf "%s", n }
    ' zotus.fa > zotus_Qiime2.fa
   
    cd /mnt/c/...

    time qiime tools import \
      --input-path USEARCH_results/zotus_Qiime2.fa \
      --output-path USEARCH_results/zotus_rep-seqs.qza \
      --type 'FeatureData[Sequence]'
  
    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-338-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/zotus_rep-seqs.qza \
      --o-classification USEARCH_results/zotus_taxonomy-silvaV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-515-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/zotus_rep-seqs.qza \
      --o-classification USEARCH_results/zotus_taxonomy-silvaV4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-338-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/zotus_rep-seqs.qza \
      --o-classification USEARCH_results/zotus_taxonomy-MiDASV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-515-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/zotus_rep-seqs.qza \
      --o-classification USEARCH_results/zotus_taxonomy-MiDASV4.qza

    biom convert \
    -i zotutab_raw.txt \
    -o zotutab_raw.biom \
    --to-hdf5 \
    --table-type="OTU table"

    qiime tools import \
      --input-path zotutab_raw.biom \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV210Format \
      --output-path zotutab_raw.qza

    qiime feature-table rarefy \
      --i-table zotutab_raw.qza \
      --p-sampling-depth 31586 \
      --o-rarefied-table zotutab_raw_rarefied.qza
      
    qiime tools export \
      --input-path zotutab_raw_rarefied.qza \
      --output-path R

    biom convert \
      -i R/feature-table.biom \
      -o R/feature_table_rarefied.tsv \
      --to-tsv

    cd R
    sed -i '1d' feature_table_rarefied.tsv
    sed -i 's/#OTU ID/Feature ID/' feature_table_rarefied.tsv
   
    qiime tools export \
      --input-path zotus_taxonomy-silvaV3V4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-silva-V3V4.tsv

    qiime tools export \
      --input-path zotus_taxonomy-silvaV4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-silva-V4.tsv
    
    qiime tools export \
      --input-path zotus_taxonomy-MiDASV3V4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-MiDAS-V3V4.tsv

    cd /mnt/c/Omics/qiime2/Rare_taxa/USEARCH_results

    qiime tools export \
      --input-path zotus_taxonomy-MiDASV4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-MiDAS-V4.tsv
    
    time qiime tools import \
      --input-path USEARCH_results/otus_Qiime2.fa \
      --output-path USEARCH_results/otus_rep-seqs.qza \
      --type 'FeatureData[Sequence]'

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-338-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/otus_rep-seqs.qza \
      --o-classification USEARCH_results/otus_taxonomy-silvaV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier silva-138-99-515-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/otus_rep-seqs.qza \
      --o-classification USEARCH_results/otus_taxonomy-silvaV4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-338-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/otus_rep-seqs.qza \
      --o-classification USEARCH_results/otus_taxonomy-MiDASV3V4.qza

    time qiime feature-classifier classify-sklearn \
      --i-classifier MiDAS-4.8.1-515-806-nb-classifier-q2_202102.qza \
      --i-reads USEARCH_results/otus_rep-seqs.qza \
      --o-classification USEARCH_results/otus_taxonomy-MiDASV4.qza

    biom convert \
    -i otutab_raw.txt \
    -o otutab_raw.biom \
    --to-hdf5 \
    --table-type="OTU table"

    qiime tools import \
      --input-path otutab_raw.biom \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV210Format \
      --output-path otutab_raw.qza

    qiime feature-table rarefy \
      --i-table otutab_raw.qza \
      --p-sampling-depth 30948 \
      --o-rarefied-table otutab_raw_rarefied.qza
    
    qiime tools export \
      --input-path otutab_raw_rarefied.qza \
      --output-path R

    biom convert \
      -i R/feature-table.biom \
      -o R/feature_table_rarefied.tsv \
      --to-tsv

    cd R
    sed -i '1d' feature_table_rarefied.tsv
    sed -i 's/#OTU ID/Feature ID/' feature_table_rarefied.tsv
  
    qiime tools export \
      --input-path otus_taxonomy-silvaV3V4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-silva-V3V4.tsv

    qiime tools export \
      --input-path otus_taxonomy-silvaV4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-silva-V4.tsv

    qiime tools export \
      --input-path otus_taxonomy-MiDASV3V4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-MiDAS-V3V4.tsv
 
    qiime tools export \
      --input-path otus_taxonomy-MiDASV4.qza \
      --output-path R

    cd R
    mv taxonomy.tsv taxonomy-MiDAS-V4.tsv                                                                                                    v taxonomy-MiDAS-V4.tsv
    
