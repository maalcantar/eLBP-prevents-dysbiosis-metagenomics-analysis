# Engineered live biotherapeutic product (eLBP) prevents antibiotic-induced dysbiosis: metagenomics analysis

This repository contains code needed to reproduce metagenomics (i.e., shotgun sequencing) analyses described in “An engineered live biotherapeutic for the prevention of antibiotic-induced dysbiosis” (Cubillos-Ruiz et al. 2022). 

# Installation & requirements  

This repository, including all code needed to reproduce analyses, can be installed using:

~~~
git clone https://github.com/maalcantar/engineered-probiotic-prevents-dysbiosis-metagenomics-analysis.git
cd engineered-probiotic-prevents-dysbiosis-metagenomics-analysis
pip install -r requirements_python.txt #requirements for Jupyter Notebooks 
~~~

R requirements:
* broom v0.7.9 
* cowplot v1.1.1 
* ggplot2 v3.3.5 
* multcomp v1.4-17 
* phyloseq v1.30.0 
* plyr v1.8.6
* rmarkdown v2.10
* tidyverse v1.3.1

Additional requirements: 
* Bowtie2 v2.2.6
* Samtools v1.6

The majority of the metagenomic sequence processing and analyses were conducted using the Commonwealth Computational Cloud for Data Driven Biology ([C3DDB](https://www.mghpcc.org/c3ddb/)) cluster.

# Directory structure

### source code

All code is in  <code>src/</code>, which contains a combination of bash, python (in Jupyter Notebooks), and R scripts. The numbering at the beginning of each file name indicates the order in which that script should be run. 

### data

All data files are found in and/or will be written to <code>data/</code>

* <code>data/card/</code>
  * contains both the native CARD database (protein homolog model version v3.0.9.) and the custom CARD database (which includes SpTEM1 B-lactamase and chloramphenicol resistance used in this study). 
* <code>data/custom_card_db_data/</code>
  * folder where ARG analyses (based in custom card database) are written to.
* <code>data/reference_genomes/</code>
  * directory where mouse reference genome should be stored.
  mouse reference genome can be downloaded from: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.26_GRCm38.p6/
   * the actual genome file is called: GCF_000001635.26_GRCm38.p6_genomic.fna.gz
* * <code>data/sequencing/</code> 
  * contains sequencing data and results for some analysis applied to sequencing data.

Metagenomics sequencing data is not included in this repository due to their large file sizes; however, all data are publicly available under NCBI Bioproject PRJNA803721.

### figures

Figures created with R can be found in figs/. The majority of figures shown in the publication were recreated using PRISM.

Please feel free to reach out if you have any questions about implementation or reproducing analyses! (alcantar [at] mit [dot] edu).
