#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH --mem=5000                          # Memory total in MB (for all cores)
#SBATCH -p defq                           # Partition to run in
#SBATCH -o create_custom_card_db_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e create_custom_card_db_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script uses bowtie to create an index database of the custom card database
# this script assumes you have the 'nucleotide_fasta_protein_homolog_model.fasta'
# file in a directory called '../data/card/card_custom_db/' (see description for '00b_create_card_ref_db.sh' for more details) 
# Note: this script works with the custom card database, which includes synthetic B-lactamase sequences used with
# engineered probiotic, as opposed to script '00b_create_card_ref_db.sh' which works with the native CARD database 

# load bowtie
module load c3ddb/bowtie2/2.2.6

# parameters
# none

# returns
# sequence index of card arg database

# build genome index from card database
bowtie2-build ../data/card/card_custom_db/nucleotide_fasta_protein_homolog_model_custom_Llactis.fasta \
../data/card/card_custom_db/nucleotide_fasta_protein_homolog_model_custom_db

