#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p defq                           # Partition to run in
#SBATCH --mem=50000                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=alcantar@mit.edu   # Email to which notifications will be sent

# this script uses bowtie to create an index database of the mouse genome
# this is required when mapping reads onto the mouse (i.e., host) genome in order to filter host reads
# the only requirement is that the mouse genome be present in the '../data/reference_genomes' directory, in a sub-directory termed 'mouse'

# parameters
# none

# returns
# sequence index of mouse genome

# load bowtie2
module load c3ddb/bowtie2/2.2.6

# mouse genome is called 'genome.fasta' in this scenario
# reference genome can be downloaded from: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.26_GRCm38.p6/
# from that link, the actual genome file that we used is called: GCF_000001635.26_GRCm38.p6_genomic.fna.gz
bowtie2-build ../data/reference_genomes/mouse/genome.fasta ../data/reference_genomes/mouse/mouse_db
