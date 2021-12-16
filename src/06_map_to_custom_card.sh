#!/bin/bash

#SBATCH -c 6                               # Request six cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 2-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=10000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o batch_map_to_custom_card_21618-5232G_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e batch_map_to_custom_card_21618-5232G_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script uses bowtie to map sequence reads onto a custom card
# database. this custom card database contains ARGs from CARD protein homolog model v3.0.9 and
# B-lactamase + chloramphenicol resistance genes used in engineered probiotics from this study

# parameters
# batch of sequences to be mapped. for example, argument '21600' will map all sequences whose id begin with the number 21600
# (i.e., D21-1216001 through D21-1216009)

# returns
# sam and bam file of mapped reads

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
      # map fastq files to custom card database
      bowtie2 -x ../data/card/card_custom_db/nucleotide_fasta_protein_homolog_model_custom_db \
      -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p 6 \
      -1 ${F}_1_host_removed.fastq -2 ${F}_2_host_removed.fastq -S ${F}_host_removed_mapped_unmapped_card_custom.sam

      echo "converting to bam"
      # convert sam to binary bam file using samtools
      samtools view -bS -@ 6 ${F}_host_removed_mapped_unmapped_card_custom.sam > ${F}_host_removed_mapped_unmapped_card_custom.bam

      rm ${F}_host_removed_mapped_unmapped_card_custom.sam

    done

done

