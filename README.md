# Prediction of traits (optimum growth temperature, salinity, etc.) for bacteria and archaea from genome properties

**WARNING: WORK IN PROGRESS**

Goal: Bacteria and archaea have evolved to live in diverse environmental conditions. One of the most defining traits of these microorganisms is the environment where it grows best. This information is usually determined through experimentation. Often experiments cannot easily be performed, like when the organism is not easily cultured in the laboratory. With today's sequencing and computational infrastructure, it is usually easier to obtain a genome from environmental samples than to obtain a strain in the laboratory for experimentation. Therefore, it is of great value to predict optimum growth conditions from genomes.

Approach: Predicting of microbial traits often starts with the functions encoded by genes in the genome. (For example, https://diaspora-project.de/). However, research has established that the composition of genomes, such as the frequencies of nucleotides, codons, and amino acids, is also under selection by environmental conditions such as temperature and salinity. Selection can act on intracellular proteins but should be greatest for proteins secreted into the extracellular environment. **Here, I evaluate the ability of select features based on those properties to predict growth conditions.** Proteins likely to be extracellular and soluble are identified with a fast heuristic (orders of magnitude faster than Signalp v.6, but half of identified proteins are false positives for secretion) to increase predictive power. Later work will include prediction from functions encoded in the genome with an emphasis on feature engineering. Machine learning methods are more powerful with more informative derived features and less prone to overfitting with fewer irrelevant features. Feature selection grounded by expertise in microbial genomics should increase the quality and scope of predictions.

Data: Predictions are trained on data from BacDive, which enables programmatic access to the DSMZ strain collection. Genomes used for training and, later, predictions are from the Genome Taxonomy Database (GTDB). Currently, BacDive has the most comprehensive set of information on microbial traits, and GTDB has the highest quality phylogenetic representation of microbial taxonomy. Note that huge portions of the microbial taxonomy have not been sampled. Predictive models in microbiology should work for new or poorly sampled groups of microorganisms. To be confident in the model's ability to do so, training and testing of the model is carried out with data held out at predefined taxonomic levels rather than a random hold out set. 

Contents of this repo allow:

1. Programmatic access of BacDive data on growth conditions
2. Derivation of genome properties based on protein amino acid and [LATER] nucleotide sequences. 
3. Feature engineering to support genome classification
4. [LATER] Modeling training and evaluation with test/train split at high taxonomic levels
5. [LATER] Prediction of optimum growth conditions using a pretrained model

# Use
## Installation

Developed with python 3.10. 

Install required packages:

```shell
# activate your virtual environment
pip install -r requirements.txt
```

## Prediction

To use a pre-trained model to predict traits for genomes ... [WORK IN PROGRESS]

## Training

To train a model from data from BacDive and the Genome Taxonomy Database, follow these steps. For speed, this repo was node coded to use other datasets, so modifications to code would be needed to do so.

### 1. Download genome data

Download genomic data and metadata from the Genome Taxonomy Database. Anticipate long downloads as the files are very large (~100G total uncompressed).

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

### 2. Download BacDive data

Register with BacDive for access to the BacDive API at `https://api.bacdive.dsmz.de/`

Run the `query_bacdive` module (runtime: about 1 sec per 100 entries):

```shell
python3 src/query_bacdive.py -c .bacdive_credentials -max 171000 -o data/bacdive_data.json
```

### 3. Creating the feature table

[LATER] _Note to self: convert module into script_

```shell
# add instructions here
```

Features can be inspected using a Dash app. Currently variables are hard-coded, but this will run
when called from the main directory if a subdirectory data has features file `./data/gtdb_v207_features_20230407.tsv`:

```shell
python3 src/app_feature_analysis.py
```

### 4. Training the model

To be continued
