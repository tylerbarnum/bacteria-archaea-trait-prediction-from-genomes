Prediction of traits (optimum growth temperature, salinity, etc.) for bacteria and archaea from genomes

# Prerequisites

Install required packages:

```shell
pip install -r requirements.txt
```

Download genomic data and metadata from the Genome Taxonomy Database

```shell
# Metadata for bacteria and archaea
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tar.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tar.gz
# Genomes and proteins (nt and amino acid)
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_proteins_nt_reps.tar.gz
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz
```

Register with BacDive for access to the BacDive API at `https://api.bacdive.dsmz.de/`