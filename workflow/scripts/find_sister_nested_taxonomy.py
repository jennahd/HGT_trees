#!/usr/bin/env python

#General modules
from ete3 import *
import numpy as np
from collections import Counter

##Functions##

def open_tree(tree):
    tree = Tree(tree)
    tree.ladderize(direction=1)
    return tree

def get_focal_list(map):
    focal_list = []
    with open(map, 'r') as f:
        for line in f:
            ID = line.split("\t")[0:1][0]
            type = line.split("\t")[2:3][0]
            if type == "focal":
                focal_list.append(ID)
    return focal_list

def get_focal_taxa_dict(tree):
    focal_taxa_dict = {}
    for node in tree.traverse():
        if node.is_leaf():
            if "FOCAL" in node.name:
                num = node.name.split("__")[0:1][0].split("_")[1:2][0]
                if num not in focal_taxa_dict.keys():
                    focal_taxa_dict[num] = []
                    focal_taxa_dict[num].append(node.name)
                else:
                    focal_taxa_dict[num].append(node.name)
    return focal_taxa_dict

def get_outgroup_list(tree, focal_list, focal_taxonomy):
	outgroup_list = []
	focal_taxonomy_taxa = []
	focal_taxa_in_tree = []
	for node in tree.traverse():
		if node.is_leaf():
			if "OUTGROUP" in node.name:
				outgroup_list.append(node.name)
			elif focal_taxonomy in node.name:
				focal_taxonomy_taxa.append(node.name)
			ID = node.name.split("|")[0:1][0]
			if "FOCAL" in ID or "OUTGROUP" in ID:
				ID = ID.split("__")[1:]
				ID = ", ".join(ID)
			if ID in focal_list:
				focal_taxa_in_tree.append(node.name)
	return outgroup_list, focal_taxonomy_taxa, focal_taxa_in_tree

def initial_root_tree(tree, focal_taxa_dict):
    tree.unroot()
    focal_taxa = [item for sublist in list(focal_taxa_dict.values()) for item in sublist]
    non_focal_taxa = np.setdiff1d(tree.get_leaf_names(), focal_taxa)
    start_out = None
    start_out_size = 0
    for node in tree.traverse():
        if len([i for i in node.get_leaf_names() if i in non_focal_taxa]) == len(node.get_leaf_names()) and len(node) >= 2:
            if len(node) > start_out_size:
                start_out = node
                start_out_size = len(node)
    if start_out != None:
        tree.set_outgroup(start_out)
    else:
         R = tree.get_midpoint_outgroup()
         tree.set_outgroup(R)
    return tree

def get_focal_node_dict(tree, focal_taxa_dict):
    focal_node_dict = {}
    for num in focal_taxa_dict.keys():
        if len(focal_taxa_dict[num]) > 1:
            node = tree.get_common_ancestor(focal_taxa_dict[num])
            focal_node_dict[num] = node
        else:
            node = (tree.get_leaves_by_name(focal_taxa_dict[num][0]))[0].up
            focal_node_dict[num] = node
    return tree, focal_node_dict

def check_up_nodes(tree, focal_node_dict, focal_taxonomy_taxa, focal_taxa_in_tree):
    focal_node_dict_up = {}
    for num, node in focal_node_dict.items():
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
        focal_node_dict_up[num] = node
    return focal_node_dict_up

def set_outgroup(tree, outgroup_list):
	out_ancestor = tree.get_common_ancestor(outgroup_list)
	out_nodes_list = out_ancestor.get_descendants()
	out_nodes_list.append(out_ancestor)
	out_nodes_dict = {}
	out_node = None
	for node in out_nodes_list:
		out_taxa = 0
		other_taxa = 0
		for leaf in node.get_leaf_names():
			if "OUTGROUP" in leaf:
				out_taxa += 1
			else:
				other_taxa += 1
		if other_taxa == 0:
			out_nodes_dict[node] = out_taxa
	out_node = sorted(out_nodes_dict, key= lambda k: len(k), reverse=True)[0]
	tree.set_outgroup(out_node)
	return tree

def set_far_node(tree, focal_node_dict, cluster_num):
	far_node = None
	far_node_dist = 0
	if focal_node_dict[cluster_num] != tree:
		tree.set_outgroup(focal_node_dict[cluster_num])
		for node in tree.traverse():
			dist = tree.get_distance(focal_node_dict[cluster_num], node, topology_only=True)
			if node not in focal_node_dict[cluster_num] and node != tree and dist > far_node_dist and len(node) >= 2:
				far_node = node
				far_node_dist = dist
	if far_node != None:
		tree.set_outgroup(far_node)
	else:
		R = tree.get_midpoint_outgroup()
		tree.set_outgroup(R)
	return tree

def find_sister_group(tree, focal_node):
	group_names = []
	sister_node = None
	No_sister = False
	try:
		sister_node = focal_node.up
		while sister_node.support < 0.7 and sister_node != tree:
			sister_node = sister_node.up
		group_names = [i for i in sister_node.get_leaf_names() if i not in focal_node.get_leaf_names()]
	except:
		No_sister = True
	return sister_node, group_names, No_sister

def get_taxonomy_list(group_node, selected_names):
	tax_list = []
	for name in selected_names:
		try:
			taxonomy = name.split("|")[2:3][0].split("@")[0:3]
			taxonomy = "@".join(taxonomy).strip("\n")
			if taxonomy != "na" or taxonomy != "":
				tax_list.append(taxonomy)
		except:
			pass
	return tax_list

def select_taxonomy(tax_list, cutoff):
	domain, domain_num = Counter([i.split("@")[0] for i in tax_list]).most_common(1)[0]
	superphylum, superphylum_num = Counter(["@".join(i.split("@", 2)[:2]) for i in tax_list]).most_common(1)[0]
	phylum, phylum_num = Counter(tax_list).most_common(1)[0]
	total = len(tax_list)
	majority_tax = ""
	if phylum_num/total >= cutoff:
		majority_tax = phylum.split("@")[0] + "\t" + phylum.split("@")[1] + "\t" + phylum.split("@")[2]
		if phylum.split("@")[2] == "":
			majority_tax = phylum.split("@")[0] + "\t" + phylum.split("@")[1] + "\t" + phylum.split("@")[1]
	elif superphylum_num/total >= cutoff:
		majority_tax = superphylum.split("@")[0] + "\t" + superphylum.split("@")[1] + "\t" + superphylum.split("@")[1]
	elif domain_num/total >= cutoff:
		majority_tax = domain.split("@")[0] + "\t" + domain.split("@")[0] + "\t" + domain.split("@")[0]
	else:
		majority_tax = "cellular_organisms\tcellular_organisms\tcellular_organisms"
	return majority_tax

################################################################################
#Implementation
cutoffs = [0.5, 0.75, 0.90, 1]

#Files
#input
map = snakemake.input.map
tree = snakemake.input.tree

#output
out = snakemake.output.out

#params
focal_taxonomy = snakemake.params.focal_taxonomy
cluster_name = snakemake.params.cluster_name
cluster_num = snakemake.params.num

try:
    tree = open_tree(tree)
except:
    #If there is no tree, then open an empty file
    open(out, "w")
else:
    #Parse mapping file to obtain list of all focal taxa
    focal_list = get_focal_list(map)
    #Extract outgroup taxa from tree, taxa that are focal based on the input mapping file, and with the focal taxonomy
    outgroup_list, focal_taxonomy_taxa, focal_taxa_in_tree = get_outgroup_list(tree, focal_list, focal_taxonomy)
    #Make dictionary of taxa labelled as "FOCAL" in the tree
    focal_taxa_dict = get_focal_taxa_dict(tree)
    #Root initially by the largest node with only "FOCAL"-labelled taxa as outgroup
    tree = initial_root_tree(tree, focal_taxa_dict)
    #Make dictionary of these nodes
    tree, focal_node_dict = get_focal_node_dict(tree, focal_taxa_dict)
    #Re-root the tree, using the outgroup if present, and if not using the furthest node from the focal node
    if len(outgroup_list) > 0:
        tree = set_outgroup(tree, outgroup_list)
    else:
        tree = set_far_node(tree, focal_node_dict, cluster_num)
    #Check nodes up in the tree structure for taxa with the focal taxonomy or from the focal list, and change focal node as needed
    focal_node_dict_up = check_up_nodes(tree, focal_node_dict, focal_taxonomy_taxa, focal_taxa_in_tree)
    with open(out, "w") as f:
        for num, focal_node in focal_node_dict_up.items():
            #get sister and nested nodes and members
            sister_node, sister_names, No_sister = find_sister_group(tree, focal_node)
            nested_node, nested_names, No_nested = find_sister_group(tree, sister_node)
            #get sister and members taxonomy
            sister_tax_list = get_taxonomy_list(sister_node, sister_names)
            nested_tax_list = get_taxonomy_list(nested_node, nested_names)
            #select taxonomy based on percentage cutoffs
            for cutoff in cutoffs:
                if len(sister_tax_list) > 0:
                    sister_majority_tax = select_taxonomy(sister_tax_list, cutoff)
                else:
                    sister_majority_tax = "de-novo\tde-novo\tde-novo"
                if len(nested_tax_list) > 0:
                    nested_majority_tax = select_taxonomy(nested_tax_list, cutoff)
                else:
                    if len(sister_tax_list) == 0 and len(nested_tax_list) == 0:
                        nested_majority_tax = "de-novo\tde-novo\tde-novo"
                    else:
                        nested_majority_tax = "none\tnone\tnone"
                f.write("%s\t%s\t%s\t%s\t%s\n" %(cluster_name, cluster_num, cutoff, sister_majority_tax, nested_majority_tax))
