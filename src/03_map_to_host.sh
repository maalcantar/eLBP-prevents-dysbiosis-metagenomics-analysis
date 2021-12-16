#!/bin/bash

#SBATCH -c 6                               # Request six cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=80000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o batch_map_to_host_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_map_to_host_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script uses bowtie to map sequence reads onto the mouse genome in order
# to filter host reads

# parameters
# batch of sequences to be mapped to mouse genome. for example, argument '21600' will map all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# sam and bam file of both mapped and unmapped reads

# load bowtie2 and SAMtools
module load c3ddb/bowtie2/2.2.6
module load c3ddb/samtools/1.6

# look inside every directory containing reads
for dir in ../data/sequencing/210616Chi/D21-$1* ; do
    echo "mapping fastq files in : $dir/"

    # get file name prefix
    find $dir -type f -name  "*_1_sequence.fastq" | \
    sed 's/_1_sequence.fastq$//' | \

    while read F; do
      # map fastq files to mouse genome
      # convert sam to binary bam file using samtools
      bowtie2 -x ../data/reference_genomes/mouse/mouse_db -p 6 \
      -1 ${F}_1_sequence_cutadapt_trimmed.fastq \
      -2 ${F}_2_sequence_cutadapt_trimmed.fastq | samtools view -bS - > ${F}_mapped_and_unmapped.bam
    done

done

