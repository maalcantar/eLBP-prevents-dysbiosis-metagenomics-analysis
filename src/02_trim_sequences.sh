#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-22:00                         # Runtime in D-HH:MM format
#SBATCH --mem=50000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o batch_trim_sequences_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_trim_sequences_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# quality trim sequences using sickle (https://github.com/najoshi/sickle)

# parameters
# batch of sequences to be quality filtered. for example, argument '21600' will quality filter all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# trimmed sequences

# look inside every directory containing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "trimming sequences in: $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do
      # quality trim at Q>20 and 30nt minimum length
      ~/sickle/sickle-master/sickle pe -f ${F}_1_sequence_cutadapt.fastq \
      -r ${F}_2_sequence_cutadapt.fastq -q 20 -l 30 -t sanger \
      -o ${F}_1_sequence_cutadapt_trimmed.fastq -p ${F}_2_sequence_cutadapt_trimmed.fastq \
       -s ${F}_trimmed_singles_file.fastq

    done

done

