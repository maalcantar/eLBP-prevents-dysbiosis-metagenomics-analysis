#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 3-00:00                         # Runtime in D-HH:MM format
#SBATCH --mem=75000                          # Memory total in MB (for all cores)
#SBATCH -p sched_mem1TB_centos7                           # Partition to run in
#SBATCH -o batch_bam_to_fq_no_host_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_bam_to_fq_no_host_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script converts sorted bam files to fastq files.
# the ultimate result is fastq files with non-host reads

# parameters
# batch of sequences to be converted to fastq. for example, argument '21600' will convert all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# fastq file with non-host reads

# load bedtools
module load c3ddb/bedtools/2.28.0

# look inside every directory containing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "converting bam to fastq files for : $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do

      # Convert bam to fastq using bedtools
      bedtools bamtofastq -i ${F}_bothEndsUnmapped_sorted.bam -fq ${F}_1_host_removed.fastq -fq2 ${F}_2_host_removed.fastq

    done

done

