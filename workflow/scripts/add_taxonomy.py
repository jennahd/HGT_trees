#!/usr/bin/env python

##Modules##
from Bio import SeqIO
from Bio.Seq import Seq

##Functions##
def get_taxids(prot):
    with open(prot, 'r') as f:
        taxid_dict = {}
        for line in f:
            accession = line.split("\t")[0:1][0]
            taxid = line.split("\t")[1:2][0].strip("\n")
            taxid_dict[accession] = taxid
        return taxid_dict

def get_taxonomy(taxid_dict, tax):
    with open(tax, 'r') as f:
        taxonomy_dict = {}
        for line in f:
            taxid = line.split("\t")[0:1][0]
            taxonomy = line.split("\t")[4:5][0].replace("cellular organisms; ", "").replace("; ", "@").replace(" ", "_").replace("(","").replace(")","").replace(":","_").replace("\\","_").replace("/","_")
            if taxid in taxid_dict.values():
                taxonomy_dict[taxid] = taxonomy
        return taxonomy_dict

def replace_fasta_headers(faa, taxid_dict, taxonomy_dict, out):
    output_file = open(out, "w")
    taxonomy = "na"
    for seq_record in SeqIO.parse(faa, "fasta"):
        accession = seq_record.description.split(" ",1)[0]
        description = seq_record.description.split(" ",1)[1].split("]",1)[0].replace(" ","_").replace("(","").replace(")","").replace(":","_").replace("\\","_").replace(";","_").replace("/","_") + "]"
        if accession in taxid_dict.keys():
            taxid = taxid_dict[accession]
            if taxid in taxonomy_dict.keys():
                taxonomy = taxonomy_dict[taxid]
        else:
            taxonomy = "na"
        new_seq_header = accession + "|" + description + "|" + taxonomy
        new_seq_header.replace(":","_").replace(";","_").replace("(","_").replace(")","_")
        new_seq_header = new_seq_header
        output_file.write(">%s\n%s\n" % (
        new_seq_header,
        seq_record.seq))

##Implementation##
tax=snakemake.input[0]
prot=snakemake.input[1]
faa=snakemake.input[2]
out=snakemake.output[0]

#get taxids
taxid_dict = get_taxids(prot)

#get taxonomy
taxonomy_dict = get_taxonomy(taxid_dict, tax)

#Replace fasta headers
replace_fasta_headers(faa, taxid_dict, taxonomy_dict, out)
