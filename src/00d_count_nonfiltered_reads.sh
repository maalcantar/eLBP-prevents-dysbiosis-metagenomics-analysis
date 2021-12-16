#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 3-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=50000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o count_nonfiltered_reads_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e count_nonfiltered_reads_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script counts number of reads in each sequencing file before any filtering has occured

# parameters
# none

# returns
# text files with i) number of forward reads ii) number of reverse reads and iii) corresponding sequence id

# creating empty text files that will store the number of reads and corresponding sequence id
echo "number of forward reads" > ../data/sequencing/fastq_read_num_nonfiltered_F.txt
echo "number of reverse reads" > ../data/sequencing/fastq_read_num_nonfiltered_R.txt
echo "sequence number" > ../data/sequencing/fastq_read_num_order_nonfiltered.txt

# look inside every directory containing sequence reads
for dir in ../data/sequencing/210616Chi/D21* ; do
    echo "counting fastq files in : $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do
      prefix=${F}
      seq_num=${prefix:(-4)}

      # write number of reads and sequence number to text file 
      echo $(cat ${prefix}_1_sequence.fastq | wc -l)/4|bc >> ../data/sequencing/fastq_read_num_nonfiltered_F.txt
      echo $(cat ${prefix}_2_sequence.fastq | wc -l)/4|bc >> ../data/sequencing/fastq_read_num_nonfiltered_R.txt
      echo 210616Chi_D21-${seq_num} >> ../data/sequencing/fastq_read_num_order_nonfiltered.txt
    done

done


