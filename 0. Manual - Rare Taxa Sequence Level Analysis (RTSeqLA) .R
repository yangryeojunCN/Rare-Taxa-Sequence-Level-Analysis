rm(list=ls())

# Run these commands one by one to install the required packages.
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("stringr")
install.packages("stringdist")
install.packages("tidyverse")

# Be sure these packages were installed.
library(Biostrings)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(stringdist)
library(tidyverse)

# Set the directory where you want to work.
setwd("C:/Omics/qiime2/data") 

# Store the feature tables, representative sequences and taxonomy into the workspace.
table.DaDa2 <- read.table("DaDa2_feature_table_rarefied.tsv", header=T, row.names= 1, sep="\t")
table.Deblur <- read.table("Deblur_feature_table_rarefied.tsv", header=T, row.names= 1, sep="\t")
table.UNOISE3 <- read.table("UNOISE_feature_table_rarefied.tsv", header=T, row.names= 1, sep="\t")
table.UPARSE <- read.table("UPARSE_feature_table_rarefied.tsv", header=T, row.names= 1, sep="\t")

seq.DaDa2 <- as.data.frame(readDNAStringSet("DaDa2_dna-sequences.fasta"))
seq.Deblur <- as.data.frame(readDNAStringSet("Deblur_dna-sequences.fasta"))
seq.UNOISE3 <- as.data.frame(readDNAStringSet("UNOISE_zotus_Qiime2.fa"))
seq.UPARSE <- as.data.frame(readDNAStringSet("UPARSE_otus_Qiime2.fa"))

tax.DaDa2 <- read.table("DaDa2_taxonomy-MiDAS-V4.tsv", header=T,  sep="\t") 
tax.Deblur <- read.table("Deblur_taxonomy-MiDAS-V4.tsv", header=T,  sep="\t") 
tax.UNOISE3 <- read.table("UNOISE_taxonomy-MiDAS-V4.tsv", header=T,  sep="\t") 
tax.UPARSE <- read.table("UPARSE_taxonomy-MiDAS-V4.tsv", header=T,  sep="\t") 

# Note: the specific amplicon-region reference database should be the same as database in taxonomy classification.
ref_region <- as.data.frame(readDNAStringSet("database/MiDAS-4.8.1-515-806-sequences.fasta"))
tax_region  <- read.table("database/MiDAS-4.8.1-515-806-classifier.tsv", header=T,  sep="\t") 

# 1.Defining the rare taxa
source("1.Defining_the_rare_taxa.R")

    ## Two question will then appear: 
    ## - Please enter the number of the samples (e.g. 3):
    ## - Please enter the occurrence of the features among all samples (e.g. 2):

# 2.Trimming representative sequences
source("2.Trimming_representative_sequences_of_rare_taxa.R")

    ## A question will then appear: 
    ## - Please enter the trimming sequence length:
    ## The trimming sequence length is the same as truncation parameter in Deblur pipeline.

# 3.Matching sequences from reference database.
source("3.Matching_sequences_from_reference_database.R")

# 4.Calculating hamming distance.
source("4.Calculating_hamming_distance.R")
P  
   ## Result can be found in folder "hamming"
