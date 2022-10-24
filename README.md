# HGT_trees
Workflow implemented using Snakemake for running phylogenetic trees to assess the evolutionary history of proteins of interest. Assesses the taxonomy of closely related sequences to determine the "source" of potential HGTs.

### Workflow outline

1. Sequences of interest are searched against the NR NCBI database
2. Hits are collated and unique sequences downloaded.
3. Hits are reduced based on sequence identity (based on the given percent cut-off, 0.7 to 0.8 recommended)
4. Taxonomy is added to sequence headers
5. Initial alignment and trimming step
6. Removal of short sequences (< 40% of the alignment length)
7. Initial fasttree and identification of sequences on long branches (i.e., outliers)
8. Second alignment, trimming, and fasttree step with the long branching sequences removed
9. Idenficiation of subtrees within the larger tree that include the "focal" sequences of interest
10. Selection of up to 3 subtrees with the most number of "focal" sequences (this number is currently hardcoded at the top of the snakefile, but can be changed manually) within the min and max bounds of numbers of taxa, alongside appropriate outgroup sequences (based on midpoint rooting of the intial large tree)
11. Alignment, trimming, and fasttree of each of the subtrees
12. Interproscan domains for subtree sequences determined and pfam domains will be included in the final tree pdfs
13. In the final tree pdf, all subtrees are coloured according to phylum (same colours across proteins included in one run), and rooted by the outgroup, with species names added and coloured according to a mapping file. The "focal" group of sequences of interest is indicated, and pfam domains are mapped to each sequence included.
14. The taxonomy of the group sister to the "focal" group and that subtending it "nested" (support >= 0.7), is then investigated at domain, superphylum, and phylum levels, and a "source" identified at majority cutoffs of 50, 75, 90, and 100 % for each.
15. The results are plotted across all proteins included in the same run.

If you would like to run more robust ML trees (for example with IQ-TREE), do so from the final trimmed alginments and then run the following bash script to create the pdfs and determine HGT source taxonomy:

scripts/run_labelling_HGT_sources.sh

*Here support needs to be >= 80 (for ultrafast bootstraps)*

## Installation

### Clone workflow
git clone https://github.com/jennahd/HGT_trees.git path/to/workdir
cd path/to/workdir

### Install Snakemake using conda
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

### Additional dependencies

**Software**
*The following should be in your path*
- diamond: https://github.com/bbuchfink/diamond
- cd-hit: https://github.com/weizhongli/cdhit
- mafft: https://mafft.cbrc.jp/alignment/software/linux.html
- trimal: http://trimal.cgenomics.org/downloads
- interproscan: https://github.com/ebi-pf-team/interproscan
- FastTree: https://anaconda.org/bioconda/fasttree
- seqtk: https://github.com/lh3/seqtk

**Python 3 modules**
- Biopython
- ete3
- palettable
- numpy
- statsmodels
- argparse
- random

**R packages**
- ggforce
- gridExtra
- ggpubr
- dplyr

### Download databases
mkdir workflow/databases

*Place the following files where you'd like on your machine and add a soft link to the database folder*

**taxonomy - only "fullnamelineage.dmp" needed, other files can be removed**
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -zxvf new_taxdump.tar.gz
```

**all nr proteins AND diamond nr blast database**
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz --db nr.dmnd
```

**NCBI nr database**
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/nr.*.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/nr.*.tar.gz.md5
for i in *.md5 ; do md5sum -c $i ; done
for i in *.tar.gz ; do tar -zxvf $i ; done
```

## Useing pipeline

### Edit config file as needed
```
vim config/example_config.yaml
mv config/example_config.yaml config/config.yaml
```

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

">Species_ID_1|accession_description MHTEFF..."

This is so that the species names and colours can be added to the final pdfs, and information about whether a species is "focal" (i.e., should be used for selecting subtrees)

### Execute workflow
activate environment:
```
conda activate snakemake_HGT_trees
```

Use a dry run first to check that everything is working:
```
snakemake -n --configfile config.yaml 
```

Run:
```
snakemake --cores #_cores --configfile config.yaml
```
