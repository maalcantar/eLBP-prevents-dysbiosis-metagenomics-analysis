#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 3-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=55000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o count_non_host_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e count_non_host_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script count the number of non-host reads in all metagenomics fastq files

# parameters
# none

# returns
# text files with number of forward and reverse reads, and a text file with corresponding sequence ids

# creating empty text files that will store the number of reads and sequence
# number corresponding to those number of reads
echo "number of forward reads" > ../data/sequencing/fastq_read_num_F_host_removed.txt
echo "number of reverse reads" > ../data/sequencing/fastq_read_num_R_host_removed.txt
echo "sequence number" > ../data/sequencing/fastq_read_num_order_host_removed.txt

# look inside every directory containing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "count non-host reads for : $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do

      prefix=${F}
      seq_num=${prefix:(-4)}
      # write number of reads and sequence number to text file
      echo $(cat ${prefix}_1_host_removed.fastq | wc -l)/4|bc >> ../data/sequencing/fastq_read_num_F_host_removed.txt
      echo $(cat ${prefix}_2_host_removed.fastq| wc -l)/4|bc >> ../data/sequencing/fastq_read_num_R_host_removed.txt
      echo 200302Chi_D20-${seq_num} >> ../data/sequencing/fastq_read_num_order_host_removed.txt
    done

done


