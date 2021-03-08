#!/usr/bin/env python

##Modules##
from ete3 import *
#import statsmodels.stats.api as sms
import numpy as np

##Functions##
def parse_mapping_file(map):
    species_name_dict = {}
    ID_colour_dict = {}
    focal_list = []
    with open(map, 'r') as f:
        for line in f:
            ID = line.split("\t")[0:1][0]
            name = line.split("\t")[1:2][0]
            type = line.split("\t")[2:3][0]
            font_colour = line.split("\t")[4:5][0].strip("\n")
            species_name_dict[ID] = name
            ID_colour_dict[ID] = font_colour
            if type == "focal":
                focal_list.append(ID)
    return species_name_dict, ID_colour_dict, focal_list

def get_cluster_annotation(info, cluster_name):
    annotation = None
    with open(info, 'r') as f:
        for line in f:
            cluster = line.split("\t")[0:1][0]
            description = line.replace("\t", "\n")
            if cluster == cluster_name:
                annotation = description
    return annotation

def set_tree_characteristics(tree):
    #initialize tree and modify general tree characteristics
    t = Tree(tree)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    t.ladderize(direction=1)
    style = NodeStyle()
    style["size"] = 0
    style["vt_line_width"] = 0.5
    style["hz_line_width"] = 0.5
    for node in t.traverse():
        node.img_style = style
    return t

def get_branch_lengths(t):
    #Get branch lengths and upper bound of interquartile range
    branch_lengths = []
    for node in t.traverse():
        if node.is_leaf():
            branch_lengths.append(node.dist)
    branch_lengths = sorted(branch_lengths)
    q1, q3= np.percentile(branch_lengths,[25,75])
    iqr = q3 - q1
    upper_bound = q3 +(1.5 * iqr)
    #quantile98 = sms.DescrStatsW(branch_lengths).quantile(0.98, return_pandas = False)
    return upper_bound

def remove_long_branching(t, upper_bound, species_name_dict, focal_list):
    #Remove long branching taxa with outlier terminal branch lengths
    long_branching = []
    keep_taxa = []
    for node in t.traverse():
        if node.is_leaf() and node.dist <= upper_bound:
            keep_taxa.append(node.name)
        elif node.is_leaf() and node.dist > upper_bound:
            ID = node.name.split("|")[0:1][0]
            if ID in focal_list:
                keep_taxa.append(node.name)
            else:
                long_branching.append(node.name)
    return keep_taxa, long_branching

def get_tree_labels(t, ID_colour_dict, species_name_dict, long_branching):
    #Add species names, and removed long branching taxa
    for node in t.traverse():
        if node.is_leaf():
            ID = node.name.split("|")[0:1][0]
            if node.name in long_branching:
                name_face = TextFace("LONG BRANCH REMOVED--" + node.name, ftype='helveticae', fsize=8, fgcolor="red", bold=True)
                node.add_face(name_face, column=0, position='branch-right')
            elif ID in species_name_dict.keys():
                name_face = TextFace(node.name + "@" + species_name_dict[ID], ftype='helveticae', fsize=8, fgcolor=ID_colour_dict[ID], bold=True)
                node.add_face(name_face, column=0, position='branch-right')
            else:
                name_face = TextFace(node.name, ftype='helveticae', fsize=7, fgcolor="black")
                node.add_face(name_face, column=0, position='branch-right')
    return t

def output_tree(t, out, annotation):
	#Output tree
	ts = TreeStyle()
	ts.show_leaf_name = False
	ts.show_branch_support = True
	ts.scale = 1500
	ts.title.add_face(TextFace(annotation,fsize=20,bold=True),0)
	t.render(out,tree_style=ts, units="mm", h=400)

################################################################################
#Implementation

#Input files
map = snakemake.input.map
info = snakemake.input.info
tree = snakemake.input.tree

#output files
out = snakemake.output.out
long = snakemake.output.long
keep = snakemake.output.keep

#Parameters
cluster_name = snakemake.params.cluster_name
#focal_taxonomy =snakemake.params.focal_taxonomy

#Parse mapping file
species_name_dict, group_colour_dict, focal_list = parse_mapping_file(map)

#Extract annotation
annotation = get_cluster_annotation(info, cluster_name)

#Tree characteristics
t = set_tree_characteristics(tree)

#Find and remove long-branching taxa
upper_bound = get_branch_lengths(t)
keep_taxa, long_branching = remove_long_branching(t, upper_bound, species_name_dict, focal_list)

#Make pdf of long-branch seqs removed
print("The number of taxa with outlier terminal branch lengths removed is %s out of %s taxa" %(len(long_branching), len(t)))

#Output sequences keeping
with open(keep, "w") as f:
    for taxon in keep_taxa:
        f.write("%s\n" %taxon)
#Output long-branching seqs
with open(long, "w") as f:
    for taxon in long_branching:
        f.write("%s\n" %taxon)

#Get tree with labels and output
t = get_tree_labels(t, group_colour_dict, species_name_dict, long_branching)
output_tree(t, out, annotation)
