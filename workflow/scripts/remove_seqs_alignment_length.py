#!/usr/bin/env python
#Modules
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Functions
def parse_mapping_file(map):
    focal_list = []
    with open(map, 'r') as f:
        for line in f:
            ID = line.split("\t")[0:1][0]
            type = line.split("\t")[2:3][0]
            if type == "focal":
                focal_list.append(ID)
    return focal_list

def read_alignment(alignment):
    alignment = AlignIO.read(open(alignment), "fasta")
    alignment_length = alignment.get_alignment_length()
    print("Alignment length: %i" % alignment_length)
    return alignment, alignment_length

def find_seqs_that_pass_length_cutoff(alignment, perc, alignment_length, focal_list):
    passed_seq_identifiers = []
    number_filled_needed = int(round(float(alignment_length)*float(perc)/100))
    print("Number of filled positions (not gaps) needed to pass length cutoff: %s" \
    % number_filled_needed)
    print(len(alignment))
    for seq in alignment:
        ID = seq.id.split("|")[0:1][0]
        position = 0
        number_gaps = 0
        number_filled = 0
        while position < alignment_length:
            if seq[position] == "-":
                number_gaps += 1
            else:
                number_filled +=1
            position += 1
        if number_filled >= number_filled_needed:
            passed_seq_identifiers.append(seq.id)
        elif ID in focal_list:
            passed_seq_identifiers.append(seq.id)
    print(" %s sequences passed the gap length filter" % len(passed_seq_identifiers))
    return passed_seq_identifiers

def output_list_passed_seqs(alignment, passed_seq_identifiers, out):
    out = open(out, "w")
    seq_list = []
    for seq in alignment:
        if seq.id in passed_seq_identifiers:
            seq_list.append(seq)
    SeqIO.write(seq_list, out, "fasta")
    out.close()

#Implementation
alignment = snakemake.input[0]
map = snakemake.input[1]
out = snakemake.output[0]
perc = snakemake.params[0]

#Get list of focal taxa
focal_list = parse_mapping_file(map)
#parse alginment and get alignment length
alignment, alignment_length = read_alignment(alignment)
#find sequences that pass the length cutoff
passed_seq_identifiers = find_seqs_that_pass_length_cutoff(alignment, perc, alignment_length, focal_list)
#output passed sequeces in a fasta file
output_list_passed_seqs(alignment, passed_seq_identifiers,out)
