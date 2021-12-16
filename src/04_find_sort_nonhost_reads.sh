#!/bin/bash

#SBATCH -c 6                               # Request six cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 3-00:00                         # Runtime in D-HH:MM format
#SBATCH --mem=75000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o batch_find_sort_nonhost_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_find_sort_nonhost_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script extracts all reads that did NOT map to the mouse reference genome,
# and sorts those reads in a manner that is compatible with Bedtools (used in next step)

# parameters
# batch of sequences to filter and sort. for example, argument '21600' will filter and sort all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# sorted non-host sequences

# load SAMtools
module load c3ddb/samtools/1.6

# look inside every directory containing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "finding unmapped and sorting fastq files in : $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do
      # find unmapped reads
      samtools view -b -f 12 -F 256 -@ 6 ${F}_mapped_and_unmapped.bam > ${F}_bothEndsUnmapped.bam

      # sorting unmapped reads so paired end files match
      samtools sort -n -@ 6 ${F}_bothEndsUnmapped.bam -o ${F}_bothEndsUnmapped_sorted.bam

      # remove intermediate file
      rm ${F}_bothEndsUnmapped.bam
    done

done

