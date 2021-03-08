# HGT_trees
Workflow implemented using Snakemake for running phylogenetic trees to assess the evolutionary history of proteins of interest. Assesses the taxonomy of closely related sequences to determine the "source" of potential HGTs.

## Workflow outline

1. Sequences of interest are searched against the NR NCBI database
2. Hits are collated and 

## Installation

### Clone workflow
git clone https://github.com/jennahd/HGT_trees.git path/to/workdir
cd path/to/workdir

### Install Snakemake using conda

### Additional dependencies

### Download databases

## Useing pipeline

### Edit config file as needed
vim config/example_config.yaml
mv config/example_config.yaml config/config.yaml

Can change the name of the config file to any name of your choice. In general it's good practice to have a different config file for each run you do of the pipeline and to save them for reproducibility.

### Edit dataset_characteristics.tsv
This file adds annotations to the final pdf tree for each of your proteins of interest

edit and rename the example file accordingly:
config/example_dataset_characteristics.tsv

The file is tab-seperated and the first field needs to correspond to the protein_ID. You can have as many additional tabs with information as you would like after the first one.

### Edit names_focal_map.tsv
This file adds species names and colours to the final pdfs, and is also used to indicate which species are "focal" and should be used for selecting subtrees.

edit and rename the example file accordingly:
example_names_focal_map.tsv

The file is tab-seperated and should include the following fields:

1. The species_ID found in the header of the fasta files
2. The species full name to be used in pdfs
3. Where the species is "focal" or "other" (i.e., whether you want it used for selecting subtrees or not)
4. The species taxonomic group, a label with this information will be added to the pdfs
5. The colour you would like the taxon labelled in the tree pdfs (could for example correspond to the taxonomic group) 

### Add faa files
Add fasta files with protein sequences of interest to a folder in the working directory where you want results written to

Each fasta file should be named:
protein_ID_1.faa
protein_ID_2.faa
protein_ID_3.faa etc.

Inside the fasta file, the headers should have the species ID seperated from the accession and description with a pipe "|":

>Species_ID_1|accession_description
MHTEFF...

This is so that the species names and colours can be added to the final pdfs, and information about whether a species is "focal" (i.e., should be used for selecting subtrees)

### Execute workflow
activate environment:
conda activate snakemake_HGT_trees

Use a dry run first to check that everything is working:
snakemake -n --configfile config.yaml 

Run:
snakemake --cores #_cores --configfile config.yaml
