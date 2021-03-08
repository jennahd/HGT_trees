#!/usr/bin/env python

##Modules##
from ete3 import *
from Bio import SeqIO
from palettable.tableau import Tableau_20
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

def pick_phyla_colours(phyla):
	#Colour most common phyla found across the entire input dataset
	phyla_colour_dict = {}
	colour_list = Tableau_20.hex_colors
	colour_list.remove("#D62728")
	colour_list.remove("#AEC7E8")
	with open(phyla, 'r') as f:
		bact_num = 0
		for line in f:
			domain = line.split(" ")[2:3][0].split("@")[0:1][0].strip("\n")
			group = line.split(" ")[2:3][0].split("@")[1:2][0].strip("\n")
			if domain == "Archaea" or group == "Archaea":
				phyla_colour_dict[group] = "#AEC7E8"
			elif domain == "Eukaryota" or group == "Eukaryota":
				phyla_colour_dict[group] = "#D62728"
			elif domain == "Bacteria":
				if bact_num < 18:
					colour = colour_list[bact_num]
					phyla_colour_dict[group] = colour
					bact_num += 1
				else:
					break
	return phyla_colour_dict

def get_cluster_annotation(info, cluster_name):
    annotation = None
    with open(info, 'r') as f:
        for line in f:
            cluster = line.split("\t")[0:1][0]
            description = line.replace("\t", "\n")
            if cluster == cluster_name:
                annotation = description
    return annotation

def open_tree(tree):
    t = Tree(tree)
    t.ladderize(direction=1)
    return t

def get_outgroup_list(t, focal_list, focal_taxonomy):
	outgroup_list = []
	focal_taxonomy_taxa = []
	focal_taxa_in_tree = []
	for node in t.traverse():
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

def make_focal_taxa_dict(t):
	focal_taxa_dict = {}
	for node in t.traverse():
		if node.is_leaf():
			if "FOCAL" in node.name:
				num = node.name.split("__")[0:1][0].split("_")[1:2][0]
				if num not in focal_taxa_dict.keys():
					focal_taxa_dict[num] = []
					focal_taxa_dict[num].append(node.name)
				else:
					focal_taxa_dict[num].append(node.name)
	return focal_taxa_dict

def initial_root_tree(t, focal_taxa_dict):
    t.unroot()
    focal_taxa = [item for sublist in list(focal_taxa_dict.values()) for item in sublist]
    non_focal_taxa = np.setdiff1d(t.get_leaf_names(), focal_taxa)
    start_out = None
    start_out_size = 0
    for node in t.traverse():
        if len([i for i in node.get_leaf_names() if i in non_focal_taxa]) == len(node.get_leaf_names()) and len(node) >= 2:
            if len(node) > start_out_size:
                start_out = node
                start_out_size = len(node)
    if start_out != None:
        t.set_outgroup(start_out)
    else:
         R = t.get_midpoint_outgroup()
         t.set_outgroup(R)
    return t

def make_focal_node_dict(t, focal_taxa_dict):
	focal_node_dict = {}
	for num in focal_taxa_dict.keys():
		if len(focal_taxa_dict[num]) > 1:
			node = t.get_common_ancestor(focal_taxa_dict[num])
			focal_node_dict[num] = node
		else:
			node = (t.get_leaves_by_name(focal_taxa_dict[num][0]))[0].up
			focal_node_dict[num] = node
	return focal_node_dict

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

def root_tree_no_out(t, focal_node):
	far_node = None
	far_node_dist = 0
	if focal_node != t:
		t.set_outgroup(focal_node)
		for node in t.traverse():
			dist = t.get_distance(focal_node, node, topology_only=True)
			if node not in focal_node and node != t and dist > far_node_dist and len(node) >= 2:
				far_node = node
				far_node_dist = dist
	if far_node != None:
		t.set_outgroup(far_node)
	else:
		R = t.get_midpoint_outgroup()
		t.set_outgroup(R)
	return t

def root_tree_out(t, outgroup_taxa):
	#if possible, root with outgroup taxa, otherwise midpoint
	t.unroot()
	#t.set_outgroup(focal_node)
	out_ancestor = t.get_common_ancestor(outgroup_taxa)
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
	t.set_outgroup(out_node)
	return t

def general_tree_settings(t):
	t.ladderize(direction=1)
    #Modify general tree characteristics
	style = NodeStyle()
	style["size"] = 0
	style["vt_line_width"] = 1
	style["hz_line_width"] = 1
	for node in t.traverse():
		node.img_style = style
	return t

def background_colours(t, node):
	nodestl = NodeStyle()
	nodestl["size"] = 0
	nodestl["vt_line_width"] = 1
	nodestl["hz_line_width"] = 1
	nodestl["bgcolor"] = "#ededed"
	node.set_style(nodestl)
	return t

def colour_node_names(t, species_name_dict, focal_taxonomy, phyla_colour_dict, group_dict):
    #Replace IDs with species names and colour taxa according to phyla dictionary
	for node in t.traverse():
		if node.is_leaf():
			ID = node.name.split("|")[0:1][0]
			if "FOCAL" in ID:
				ID = ID.split("__")[1:]
				ID = ", ".join(ID)
			#Add species names
			if ID in species_name_dict.keys():
				new_name = node.name + "@" + species_name_dict[ID]
				name_face = TextFace(new_name, ftype='helveticae', fsize=9, fgcolor="black", bold=True)
				group_face = TextFace("@" + group_dict[ID][0], ftype='helveticae', fsize=9, fgcolor=group_dict[ID][1], bold=True)
				node.add_face(name_face, column=0, position='branch-right')
				node.add_face(group_face, column=1, position='branch-right')
			elif focal_taxonomy in node.name:
				name_face = TextFace(node.name, ftype='helveticae', fsize=9, fgcolor="black", bold=True)
				group_face = TextFace("@" + focal_taxonomy, ftype='helveticae', fsize=9, fgcolor="#636363", bold=True)
				node.add_face(name_face, column=0, position='branch-right')
				node.add_face(group_face, column=1, position='branch-right')
			#Add phyla colours
			else:
				if "@" in node.name:
					try:
						taxa_id = node.name.split("|")[2:3][0].split("@")[1:2][0]
					except:
						name_face = TextFace(node.name, ftype='helveticae', fsize=9, fgcolor="black", bold=True)
						node.add_face(name_face, column=0, position='branch-right')
					else:
						if taxa_id in phyla_colour_dict.keys():
							name_face = TextFace(node.name, ftype='helveticae', fsize=9, fgcolor=phyla_colour_dict[taxa_id], bold=True)
							node.add_face(name_face, column=0, position='branch-right')
						else:
							name_face = TextFace(node.name, ftype='helveticae', fsize=9, fgcolor="black", bold=True)
							node.add_face(name_face, column=0, position='branch-right')
				else:
					name_face = TextFace(node.name, ftype='helveticae', fsize=9, fgcolor="black", bold=True)
					node.add_face(name_face, column=0, position='branch-right')
	return t

def parse_interproscan_output(t, domains):
    domain_dict = {}
    domain_desc_dict = {}
    colours = ["#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
    "#6a3d9a","#ffff99","#b15928","#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
    "#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928", "#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
    "#ffff99","#b15928","#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
    "#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928", "#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
    "#ffff99","#b15928","#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
    "#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928", "#ffffff","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99"]
    col_num = 0
    with open(domains, 'r') as f:
        for line in f:
            taxon = line.split("\t")[0:1][0]
            domain_type = line.split("\t")[3:4][0]
            domain = line.split("\t")[4:5][0]
            desc = line.split("\t")[5:6][0]
            start = line.split("\t")[6:7][0]
            stop = line.split("\t")[7:8][0]
            if domain_type == "Pfam":
                if taxon not in domain_dict.keys():
                    domain_dict[taxon] = []
                    domain_dict[taxon].append((domain, start, stop))
                else:
                    domain_dict[taxon].append((domain, start, stop))
                if domain not in domain_desc_dict.keys():
                    if col_num <= 45:
                        domain_desc_dict[domain] = (colours[col_num], desc)
                        col_num += 1
                    else:
                        domain_desc_dict[domain] = ("ffffff", desc)
    return domain_dict, domain_desc_dict

def get_tree_names(t, domain_dict):
	tree_names_dict = {}
	names = t.get_leaf_names()
	for taxon in domain_dict.keys():
		if taxon not in names:
			for name in names:
				if name in taxon:
					tree_names_dict[taxon] = name
		else:
			tree_names_dict[taxon] = taxon
	return tree_names_dict

def faa_dict(faa):
	faa_dict = {}
	for seq_record in SeqIO.parse(faa, "fasta"):
		taxon = seq_record.description
		seq = str(seq_record.seq)
		faa_dict[taxon] = seq
	return faa_dict

def map_domains_to_tree(t, faa_dict, domain_dict, domain_desc_dict, tree_names_dict):
    for taxon in domain_dict.keys():
        domains = sorted(domain_dict[taxon], key=lambda x: x[1])
        motifs = []
        for domain in domains:
            label = "arial|10|black|" + domain[0]
            motif = [int(domain[1]), int(domain[2]), "[]", None, 10, "black", domain_desc_dict[domain[0]][0], label]
            motifs.append(motif)
        try:
            seqFace = SeqMotifFace(faa_dict[taxon], motifs=motifs, seq_format="-")
            (t & tree_names_dict[taxon]).add_face(seqFace, 0, "aligned")
        except:
            pass
    return t

def output_tree(t, out, annotation, domain_desc_dict):
	domain_description = ""
	for domain, info in domain_desc_dict.items():
		domain_description += domain + ": " + info[1] + "\n"
	#print(domain_description)
	#Output tree
	ts = TreeStyle()
	ts.show_leaf_name = False
	ts.show_branch_support = True
	ts.scale = 1500
	ts.title.add_face(TextFace(annotation,fsize=20,bold=True),0)
	ts.aligned_header.add_face(TextFace(domain_description,fsize=18,bold=True), 0)
	t.render(out,tree_style=ts, units="mm", h=400)

################################################################################
#Implementation

#Files
#input
map = snakemake.input.map
info = snakemake.input.info
phyla = snakemake.input.phyla
domains = snakemake.input.domains
tree = snakemake.input.tree
faa = snakemake.input.faa
#output
out = snakemake.output.out
#params
cluster_name = snakemake.params.cluster_name
focal_taxonomy = snakemake.params.focal_taxonomy
num = snakemake.params.num

#Parse chlamydiae species names and colours
species_name_dict, group_dict, focal_list = parse_mapping_file(map)

#Get list of phyla and assign colours
phyla_colour_dict = pick_phyla_colours(phyla)

#Extract annotation
annotation = get_cluster_annotation(info, cluster_name)

#Parse/label tree
#If there is no input tree, make a dummy file for the snakemake workflow
try:
    t = open_tree(tree)
except:
    open(out,'w')
else:
    #Extract outgroup taxa from tree, taxa that are focal based on the input mapping file, and with the focal taxonomy
    outgroup_taxa, focal_taxonomy_taxa, focal_taxa_in_tree = get_outgroup_list(t, focal_list, focal_taxonomy)
    #Make dictionary of taxa labelled as "FOCAL" in the tree
    focal_taxa_dict = make_focal_taxa_dict(t)
    #Root initially by the largest node with only "FOCAL"-labelled taxa as outgroup
    t = initial_root_tree(t, focal_taxa_dict)
    #Make dictionary of these nodes
    focal_node_dict = make_focal_node_dict(t, focal_taxa_dict)
    #Re-root the tree, using the outgroup if present, and if not using the furthest node from the focal node
    if len(outgroup_taxa) == 0:
        t = root_tree_no_out(t, focal_node_dict[num])
    else:
        t = root_tree_out(t, outgroup_taxa)
    #Check nodes up in the tree structure for taxa with the focal taxonomy or from the focal list, and change focal node as needed
    focal_node_dict_up = check_up_nodes(t, focal_node_dict, focal_taxonomy_taxa, focal_taxa_in_tree)
    #Set tree settings
    t = general_tree_settings(t)
    #Colour the background of the focal node
    for num in focal_node_dict_up:
        t = background_colours(t, focal_node_dict_up[num])
    #Colour taxa names, both accoriding to the mapping file, and phyla colour assignments
    t = colour_node_names(t, species_name_dict, focal_taxonomy, phyla_colour_dict, group_dict)
    #Parse Pfam domain annotations and make dictionary connected to tree names
    domain_dict, domain_desc_dict = parse_interproscan_output(t, domains)
    tree_names_dict = get_tree_names(t, domain_dict)
    #Parse protein faa file use to make tree for placement of domains
    faa_dict = faa_dict(faa)
    #Map protein domains onto tree
    t = map_domains_to_tree(t, faa_dict, domain_dict, domain_desc_dict, tree_names_dict)
    #Write pdf of tree
    output_tree(t, out, annotation, domain_desc_dict)
