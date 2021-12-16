#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=80000                          # Memory total in MB (for all cores)
#SBATCH -p defq                          # Partition to run in
#SBATCH -o batch_remove_adapters_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_remove_adapters_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# removes nextera transposase adapter sequences from metagenomic sequencing data using cutadapt (https://cutadapt.readthedocs.io/en/stable/installation.html)

# parameters
# batch of sequences to be adapter trimmed. for example, argument '21600' will trim all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# adapter filtered sequences

# need to load miniconda since cutadapt is a python package
module load c3ddb/miniconda/3.7
source activate meta_seq # environment containing cutadapt package

# look inside every directory containing sequencing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "removing adapters from: $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do
      # trim adapters and remove sequences with length <36nt (minimum length chosen based on the fact that there was non-neglible nextera adapter sequences as early as 36nt)
      cutadapt --minimum-length=36 \
      -a CTGTCTCTTATA -A CTGTCTCTTATA \
      -o ${F}_1_sequence_cutadapt.fastq -p ${F}_2_sequence_cutadapt.fastq  \
      ${F}_1_sequence.fastq ${F}_2_sequence.fastq

    done

done

