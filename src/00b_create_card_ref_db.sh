#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=50000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o create_card_db_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e create_card_db_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script uses bowtie to create an index database of card antibiotic resistance genes
# this script assumes you have the 'nucleotide_fasta_protein_homolog_model.fasta' file in a directory called '../data/card'
# Note: we worked with Comprehensive Antibiotic Resistance Database (CARD) protein homolog model version v3.0.9.
# which can be downloaded at: https://card.mcmaster.ca/download/0/broadstreet-v3.0.9.tar.bz2
# this script creates an index using only the native sequences found in that database, as opposed to script '00c_create_card_custom_ref_db.sh' 
# which creates an index with the custom card database, which includes synthetic B-lactamase sequences and chloramphenicol resistance sequences from our engineered probiotics

# parameters
# none

# returns
# sequence index of card arg database

# load bowtie
module load c3ddb/bowtie2/2.2.6

# build genome index from card database
bowtie2-build ../data/card/nucleotide_fasta_protein_homolog_model.fasta \
../data/card/nucleotide_fasta_protein_homolog_model_db

