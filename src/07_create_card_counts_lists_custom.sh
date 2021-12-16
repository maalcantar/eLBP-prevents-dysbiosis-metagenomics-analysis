#!/bin/bash

#SBATCH -c 6                               # Request six cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH --mem=10000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o batch_count_card_mapped_reads_custom_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_count_card_mapped_reads_custom_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script generates a gene count matrix for reads that mapped against sequences
# in the CARD database

# parameters
# batch of sequences for which you want to create gene count matrix. for example, argument '21600' will create matrices for all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# matrix of gene counts for each sample -- also, gene list and count list

# load SAMtools
module load c3ddb/samtools/1.6

# look inside every directory containing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "extracting mapped reads fastq files in : $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do
      # extract reads mapped to CARD database
      # if both reads in a pair map, count as 1. if only 1 read in a pair maps, also count as 1
      # adapted from: https://github.com/karkman/MetagenomeCourse2019/tree/master/Day1#humann2
      samtools view -@ 6 -h ${F}_host_removed_mapped_unmapped_card_custom.bam | awk '$7!="=" || ($7=="=" && and($2,0x40)) {print $0}' \
        | samtools view -Su  - \
        | samtools sort -@ 6 -o ${F}_host_removed_mapped_to_card_custom_sorted.bam

      # creating index file
      echo "creating index file reads"
      samtools index ${F}_host_removed_mapped_to_card_custom_sorted.bam

      # find reads mapped to each gene
      echo "calculating gene counts"
      samtools idxstats ${F}_host_removed_mapped_to_card_custom_sorted.bam > ${F}_idxstats_custom.tsv

      prefix=${F}
      seq_num=${prefix:(-4)}

      # extract gene names (column 1) and corresponding counts (column 3)
      cat ${F}_idxstats_custom.tsv | grep -v "*" | cut -f1 > ${F}_gene_names_custom.tsv
      cat ${F}_idxstats_custom.tsv | grep -v "*" | cut -f3 > ${F}_counts_custom.tsv

      # create gene count matrix
      paste ${F}_gene_names_custom.tsv ${F}_counts_custom.tsv > ../data/sequencing/gene_count_matrices/210616Chi_D21-${seq_num}_genemat_custom.tsv
    done

done

