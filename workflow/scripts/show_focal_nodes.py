#!/usr/bin/env python

##Modules##
import random
from ete3 import *
import numpy as np

##Functions##
def parse_mapping_file(map):
    #make dictionary for replacing IDs with species names and labelling groups with colours
    species_name_dict = {}
    group_dict = {}
    focal_list = []
    with open(map, 'r') as f:
        for line in f:
            ID = line.split("\t")[0:1][0]
            name = line.split("\t")[1:2][0]
            type = line.split("\t")[2:3][0]
            group = line.split("\t")[3:4][0]
            colour = line.split("\t")[4:5][0].strip("\n")
            species_name_dict[ID] = name
            group_dict[ID] = (group, colour)
            if type == "focal":
                focal_list.append(ID)
    return species_name_dict, group_dict, focal_list

def get_cluster_annotation(info, cluster_name):
    annotation = None
    with open(info, 'r') as f:
        for line in f:
            cluster = line.split("\t")[0:1][0]
            description = line.replace("\t", "\n")
            if cluster == cluster_name:
                annotation = description
    return annotation

def find_all_focal_taxa_in_tree(tree, focal_list, focal_taxonomy):
    t = Tree(tree)
    t.unroot()
    focal_taxa_in_tree = []
    focal_taxonomy_taxa = []
    for node in t.traverse():
        if node.is_leaf():
            ID = node.name.split("|")[0:1][0]
            if ID in focal_list:
                focal_taxa_in_tree.append(node.name)
            elif focal_taxonomy in node.name:
                focal_taxonomy_taxa.append(node.name)
    return t, focal_taxa_in_tree, focal_taxonomy_taxa

def find_focal_only_nodes(t, focal_taxa_in_tree):
    focal_only_nodes = []
    candidates = []
    for node in t.traverse():
        names = node.get_leaf_names()
        focal = len([i for i in focal_taxa_in_tree if i in names])
        if focal > 0 and focal == len(names):
                focal_only_nodes.append(node)
    return focal_only_nodes

def find_focal_only_unique_nodes(t, focal_only_nodes):
    focal_only_nodes_unique = []
    for node in focal_only_nodes:
        if any(x in node.get_ancestors() for x in focal_only_nodes) == False:
            focal_only_nodes_unique.append(node)
    return focal_only_nodes_unique

def set_root(t, focal_only_nodes_unique):
    largest_node = focal_only_nodes_unique[0]
    for node in focal_only_nodes_unique:
        if len(node) > len(largest_node):
            largest_node = node
    far_node = None
    far_node_dist = 0
    try:
        t.set_outgroup(largest_node)
        for node in t.traverse():
            if node not in focal_only_nodes_unique and node != t:
                if t.get_distance(largest_node, node, topology_only=True) > far_node_dist and len(node) >= 2:
                    far_node = node
                    far_node_dist = t.get_distance(largest_node, node, topology_only=False)
        t.set_outgroup(far_node)
    except:
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
    return t

def find_nodes_of_interest(t, focal_only_nodes_unique, focal_taxa_in_tree, focal_taxonomy_taxa):
    focal_nodes = []
    for node in focal_only_nodes_unique:
        no_tax = False
        while no_tax == False:
            try:
                check = node.up
                up_names = [i for i in check.get_leaf_names() if i not in node.get_leaf_names()]
                if len([i for i in up_names if i in (focal_taxa_in_tree + focal_taxonomy_taxa)]) >= len(up_names)/4:
                    node = check
                else:
                    check2 = node.up.up
                    up_names = [i for i in check2.get_leaf_names() if i not in node.get_leaf_names()]
                    if len([i for i in up_names if i in (focal_taxa_in_tree + focal_taxonomy_taxa)]) >= len(up_names)/4:
                        node = check2
                    else:
                        no_tax = True
            except:
                no_tax = True
        if node not in focal_nodes:
            focal_nodes.append(node)
    return focal_nodes

def find_unique_focal_nodes(t, focal_nodes, focal_taxa_in_tree, focal_taxonomy_taxa):
    focal_nodes_unique = []
    for node in focal_nodes:
        if any(x in node.get_ancestors() for x in focal_nodes) == False:
            length = len([i for i in node.get_leaf_names() if i in (focal_taxa_in_tree + focal_taxonomy_taxa)])
            focal_nodes_unique.append((node, length))
    focal_nodes_unique = [i for i in focal_nodes_unique if i[1] >= 3]
    if len(focal_nodes_unique) > 0:
        focal_nodes_unique.sort(key=lambda k: k[1], reverse=True)
        focal_nodes_unique = list(list(zip(*focal_nodes_unique))[0])
    return focal_nodes_unique

def get_focal_taxa(focal_nodes, focal_taxa_in_tree, focal_taxonomy_taxa):
    focal_taxa_list = []
    taxonomy_list = []
    for node in focal_nodes:
        for taxon in node.get_leaf_names():
            if taxon in focal_taxa_in_tree:
                focal_taxa_list.append(taxon)
            elif taxon in focal_taxonomy_taxa:
                taxonomy_list.append(taxon)
    return focal_taxa_list, taxonomy_list

def set_tree_style(t):
    t.ladderize(direction=1)
    style = NodeStyle()
    style["size"] = 0
    style["vt_line_width"] = 0.5
    style["hz_line_width"] = 0.5
    for node in t.traverse():
        node.img_style = style
    return t

def label_node(t, focal_nodes):
    #label subtree
    nodestl = NodeStyle()
    nodestl["bgcolor"] = "#e6e6e6"
    nodestl["size"] = 0
    nodestl["vt_line_width"] = 0.5
    nodestl["hz_line_width"] = 0.5
    for node in focal_nodes:
        node.set_style(nodestl)
    return t

def label_tree_taxa(t, group_dict, species_name_dict, focal_taxonomy):
    #Add species names, label focal taxa, outgroup taxa
    for node in t.traverse():
        if node.is_leaf():
            ID = node.name.split("|")[0:1][0]
            if ID in species_name_dict.keys():
                new_name = node.name + "@" + species_name_dict[ID]
                name_face = TextFace(new_name, ftype='helveticae', fsize=8, fgcolor="black", bold=True)
                group_face = TextFace("@" + group_dict[ID][0], ftype='helveticae', fsize=8, fgcolor=group_dict[ID][1], bold=True)
                node.add_face(name_face, column=0, position='branch-right')
                node.add_face(group_face, column=1, position='branch-right')
            elif focal_taxonomy in node.name:
                name_face = TextFace(node.name, ftype='helveticae', fsize=8, fgcolor="black", bold=True)
                group_face = TextFace("@" + focal_taxonomy, ftype='helveticae', fsize=8, fgcolor="#636363", bold=True)
                node.add_face(name_face, column=0, position='branch-right')
                node.add_face(group_face, column=1, position='branch-right')
            else:
                name_face = TextFace(node.name, ftype='helveticae', fsize=7, fgcolor="black")
                node.add_face(name_face, column=0, position='branch-right')
    return t

def output_tree(t, overview, annotation):
    #Output tree
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.scale = 1500
    ts.title.add_face(TextFace(annotation,fsize=20,bold=True),0)
    t.render(overview,tree_style=ts, units="mm", h=400)

def output_taxa(taxa, focal_nodes, focal_taxa_list, taxonomy_list):
    with open(taxa, "w") as f:
        count = 0
        for node in focal_nodes:
            for taxon in node:
                if taxon.name in focal_taxa_list:
                    f.write("%s\tinput\t%s\n" %(count, taxon.name))
                elif taxon.name in taxonomy_list:
                    f.write("%s\ttaxonomy\t%s\n" %(count, taxon.name))
            count += 1

################################################################################
#Implementation

#input files
map = snakemake.input.map
info = snakemake.input.info
tree = snakemake.input.tree

#output files
overview = snakemake.output.overview
taxa = snakemake.output.taxa

#parameters
cluster_name = snakemake.params.cluster_name
focal_taxonomy = snakemake.params.focal_taxonomy

#Parse mapping file
species_name_dict, group_dict, focal_list = parse_mapping_file(map)

#Extract annotation
annotation = get_cluster_annotation(info, cluster_name)

#make list of focal taxa in tree
t, focal_taxa_in_tree, focal_taxonomy_taxa = find_all_focal_taxa_in_tree(tree, focal_list, focal_taxonomy)

#Find nodes including focal taxa of interest
focal_only_nodes = find_focal_only_nodes(t, focal_taxa_in_tree)
focal_only_nodes_unique = find_focal_only_unique_nodes(t, focal_only_nodes)
t = set_root(t, focal_only_nodes_unique)
focal_nodes = find_nodes_of_interest(t, focal_only_nodes_unique, focal_taxa_in_tree, focal_taxonomy_taxa)
focal_nodes_unique = find_unique_focal_nodes(t, focal_nodes, focal_taxa_in_tree, focal_taxonomy_taxa)

if len(focal_nodes_unique) == 0:
    #If there is no focal node identified output tree without, for example if there are no clades with at least 3 taxa from the focal list or with focal taxonomy
    t = set_tree_style(t)
    t = label_tree_taxa(t, group_dict, species_name_dict, focal_taxonomy)
    output_tree(t, overview, annotation)
    open(taxa,'w')
else:
    #Find focal taxa
    focal_taxa_list, taxonomy_list = get_focal_taxa(focal_nodes_unique, focal_taxa_in_tree, focal_taxonomy_taxa)
    #Output tree
    t = set_tree_style(t)
    t = label_node(t, focal_nodes_unique)
    t = label_tree_taxa(t, group_dict, species_name_dict, focal_taxonomy)
    output_tree(t, overview, annotation)
    output_taxa(taxa, focal_nodes_unique, focal_taxa_list, focal_taxonomy_taxa)
