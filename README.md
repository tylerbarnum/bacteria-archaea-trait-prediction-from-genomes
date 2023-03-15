# Prediction of traits (optimum growth temperature, salinity, etc.) for bacteria and archaea from genome properties

Microorganisms have evolved to live in diverse environmental conditions. One of the most defining traits of a microorganism is the environment where it grows best, usually determined through experimentation. When experiments cannot easily be performed - such as when the organism is not easily cultured in the laboratory - it is of great interest to predict these optimum growth conditions from genomes. (Now, obtaining a genome from environmental samples is most often easier than obtaining a strain in the laboratory). 

Predicting of microbial traits often starts with the functions encoded by genes in the genome. However, research has established that the composition of genomes, such as the frequencies of nucleotides, codons, and amino acids, is also under selection by environmental conditions such as temperature and salinity. Here, I evaluate the ability of select features based on those properties to predict growth conditions as reported in BacDive.

**WARNING: WORK IN PROGRESS**

Tools include:

1. Programmatic access of BacDive data on growth conditions
2. Derivation of genome properties based on protein amino acid and [LATER] nucleotide sequences. 
3. [LATER] Feature engineering and test/train split (at high taxonomic levels) to support genome classification
4. [LATER] Modeling training and evaluation

## Prerequisites

Developed with python 3.10. 

Install required packages:

```shell
pip install -r requirements.txt
```

Download genomic data and metadata from the Genome Taxonomy Database

```shell
GENOME_DIR=/path/to/data
cd $GENOME_DIR
# Metadata for bacteria and archaea
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tar.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tar.gz
# Genomes and proteins (nt and amino acid)
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_proteins_nt_reps.tar.gz
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
```

Register with BacDive for access to the BacDive API at `https://api.bacdive.dsmz.de/`