#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH --mem=5000                         # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o final_gene_count_mat_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e final_gene_count_mat_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script create the final gene count matrix that will be used for analyses

# parameters
# none

# returns
# gene count matrix for all samples

# create necessary directories if they do not already exist
if [ ! -d "../data/sequencing/gene_count_matrices_for_concat_custom/" ]
  then
    mkdir ../data/sequencing/gene_count_matrices_for_concat_custom/
fi

if [ ! -d "../data/sequencing/gene_count_matrices_for_concat_custom/final_counts_dir_custom/" ]
  then
    mkdir ../data/sequencing/gene_count_matrices_for_concat_custom/final_counts_dir_custom/
fi

# copy gene counts
cp ../data/sequencing/210616Chi/D21-*/*_counts_custom.tsv \
../data/sequencing/gene_count_matrices_for_concat_custom/

# copy gene names (same order for all, so only need to copy one)
cp ../data/sequencing/210616Chi/D21-216001-5171G/210616Chi_D21-216001_gene_names_custom.tsv \
../data/sequencing/gene_count_matrices_for_concat_custom/210616Chi_D21-216001_gene_names_custom.tsv

# get gene lengths
cat ../data/sequencing/210616Chi/D21-216001-5171G/210616Chi_D21-216001_idxstats_custom.tsv | grep -v "*" | cut -f2 > \
../data/sequencing/gene_count_matrices_for_concat_custom/final_counts_dir_custom/GENE_LENGTHS_custom.tsv

# final count matrix
paste ../data/sequencing/gene_count_matrices_for_concat_custom/210616Chi_D21-216001_gene_names_custom.tsv \
../data/sequencing/gene_count_matrices_for_concat_custom/*_counts_custom.tsv > \
../data/sequencing/gene_count_matrices_for_concat_custom/final_counts_dir_custom/gene_mat_all_no_len_custom.tsv

