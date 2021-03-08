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

def make_focal_taxa_dict(taxa):
    focal_taxa_dict = {}
    with open(taxa, 'r') as f:
        for line in f:
            num = line.split("\t")[0:1][0]
            type = line.split("\t")[1:2][0]
            ID = line.split("\t")[2:3][0].strip("\n")
            if num not in focal_taxa_dict.keys():
                focal_taxa_dict[num] = []
                focal_taxa_dict[num].append(ID)
            else:
                focal_taxa_dict[num].append(ID)
    return focal_taxa_dict

def initial_root_tree(tree, focal_taxa_dict):
    t = Tree(tree)
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
    return t, focal_node_dict

def root_tree(t, focal_node):
    far_node = None
    far_node_dist = 0
    try:
        t.set_outgroup(focal_node)
        for node in t.traverse():
            dist = t.get_distance(focal_node, node, topology_only=False)
            if dist > far_node_dist and len(node) >= 2 and node not in focal_node:
                far_node = node
                far_node_dist = dist
        t.set_outgroup(far_node)
    except:
         R = t.get_midpoint_outgroup()
         t.set_outgroup(R)
    return t

def subtree_candidate_selection(t, focal_node):
    subtree_candidates = {}
    ancestral_nodes = focal_node.get_ancestors()
    ancestral_nodes.remove(t)
    for node in ancestral_nodes:
        num_taxa = len(node.get_leaf_names())
        count = 0
        step_node = focal_node
        while node != step_node:
            original_length = len(step_node)
            step_node = step_node.up
            new_length = len(step_node)
            if new_length - original_length > 1:
                count += 1
        if count >= 3 and node.support >= 0.7:
            subtree_candidates[node] = (node.support, node.dist, num_taxa, count)
        if len(subtree_candidates) == 0:
            if count >= 3:
                subtree_candidates[node] = (node.support, node.dist, num_taxa, count)
    if len(subtree_candidates) == 0:
        subtree_candidates[t] = (t.support, t.dist, len(t), 0)
    return subtree_candidates

def refine_subtree_candidates(subtree_candidates, focal_taxa, max_taxa, min_taxa):
    refined_subtree_candidates = {}
    min_limit = (min_taxa + len(focal_taxa))
    max_limit = (max_taxa + len(focal_taxa))
    refined_subtree_candidates = {k: v for k, v in subtree_candidates.items() if v[2] >= min_limit and v[2] <= max_limit}
    if len(refined_subtree_candidates) == 0:
        node_dict = {k: v for k, v in subtree_candidates.items() if v[2] >= min_limit}
        sorted_num = sorted(node_dict.items(), key= lambda k: k[1][2])
        refined_subtree_candidates[sorted_num[0][0]] = sorted_num[0][1]
    return refined_subtree_candidates

def select_subtree(refined_subtree_candidates):
    subtree = None
    sorted_support = sorted(refined_subtree_candidates.items(), key= lambda k: k[1][0], reverse=True)
    subtree = sorted_support[0][0]
    subtree_taxa = subtree.get_leaf_names()
    return subtree, subtree_taxa

def find_distant_nodes(subtree, subtree_taxa, focal_node):
    distant_nodes = {}
    removed_distant_taxa = []
    for node in subtree.traverse():
        dist = subtree.get_distance(focal_node, node, topology_only=False)
        count = 0
        common_node = node.get_common_ancestor(focal_node)
        step_node = node
        while common_node != step_node:
            step_node = step_node.up
            count += 1
        if count >= 6 and node.support >= 0.7 and len(node) >= 5 and node not in focal_node:
            distant_nodes[node] = (node.support, dist, len(node), count)
    return distant_nodes

def unique_distant_nodes(distant_nodes):
    unique_distant_nodes = {}
    for node in distant_nodes.keys():
        length = 0
        for node2 in distant_nodes.keys():
            if node in node2.get_descendants():
                length +=1
        if length == 0:
            unique_distant_nodes[node] = distant_nodes[node]
    unique_distant_nodes = sorted(unique_distant_nodes.items(), key= lambda k: k[1][3], reverse=True)
    return unique_distant_nodes

def remove_distant_nodes(subtree_taxa, unique_distant_nodes, max_taxa, removed_taxa, focal_taxa):
    max_limit = (max_taxa + len(focal_taxa))
    count = 0
    while (len(subtree_taxa) - len(removed_taxa)) > max_limit:
        if count < len(unique_distant_nodes):
            node = unique_distant_nodes[count][0]
            for taxon in node.get_leaf_names():
                removed_taxa.append(taxon)
            count += 1
        else:
            break
    return removed_taxa

def get_outgroup_seqs(subtree, subtree_taxa):
    outgroup_taxa = []
    upnode = subtree
    while len(outgroup_taxa) < 20:
        upnode = upnode.up
        new_taxa = []
        if upnode != None:
            for leaf in upnode:
                if leaf.name not in subtree_taxa:
                    if leaf.name not in outgroup_taxa:
                        new_taxa.append(leaf.name)
            if len(new_taxa) > 50:
                new_taxa = random.sample(new_taxa, k=50)
        else:
            break
        outgroup_taxa.extend(new_taxa)
    return outgroup_taxa

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

def label_subtree(t, subtree):
    #label subtree
    nodestl = NodeStyle()
    nodestl["bgcolor"] = "#cdced1"
    nodestl["size"] = 0
    nodestl["vt_line_width"] = 0.5
    nodestl["hz_line_width"] = 0.5
    subtree.set_style(nodestl)
    return t

def label_refined_candidates(t, subtree, candidates):
    description = "Support: %s\nDistance: %s\nTaxa: %s\nUp: %s" %(candidates[subtree][0], candidates[subtree][1], candidates[subtree][2], candidates[subtree][3])
    node_face = TextFace(description, ftype='helveticae', fsize=10, bold=True)
    node_face.background.color = "white"
    node_face.margin_bottom = 2
    node_face.border.width = 1
    subtree.add_face(node_face, column=0, position = "branch-top")
    return t

def label_tree_taxa(t, group_dict, species_name_dict, focal_taxa, focal_taxonomy, outgroup_taxa, removed_taxa):
    #Add species names, label focal taxa, outgroup taxa
    for node in t.traverse():
        if node.is_leaf():
            ID = node.name.split("|")[0:1][0]
            if node.name in outgroup_taxa:
                name_face = TextFace(node.name, ftype='helveticae', fsize=8, fgcolor="green", bold=True)
                node.add_face(name_face, column=0, position='branch-right')
            elif node.name in removed_taxa:
                name_face = TextFace(node.name, ftype='helveticae', fsize=8, fgcolor="blue", bold=True)
                node.add_face(name_face, column=0, position='branch-right')
            elif ID in species_name_dict.keys():
                new_name = node.name + "@" + species_name_dict[ID]
                if node.name in focal_taxa:
                    new_name = "FOCAL__" + new_name
                name_face = TextFace(new_name, ftype='helveticae', fsize=8, fgcolor="black", bold=True)
                group_face = TextFace("@" + group_dict[ID][0], ftype='helveticae', fsize=8, fgcolor=group_dict[ID][1], bold=True)
                node.add_face(name_face, column=0, position='branch-right')
                node.add_face(group_face, column=1, position='branch-right')
            elif focal_taxonomy in node.name:
                new_name = node.name
                if node.name in focal_taxa:
                    new_name = "FOCAL__" + new_name
                name_face = TextFace(new_name, ftype='helveticae', fsize=8, fgcolor="black", bold=True)
                group_face = TextFace("@" + focal_taxonomy, ftype='helveticae', fsize=8, fgcolor="#636363", bold=True)
                node.add_face(name_face, column=0, position='branch-right')
                node.add_face(group_face, column=1, position='branch-right')
            else:
                name_face = TextFace(node.name, ftype='helveticae', fsize=7, fgcolor="black")
                node.add_face(name_face, column=0, position='branch-right')
    return t


def label_tree_taxa_simple(t, group_dict, species_name_dict, focal_taxa, focal_taxonomy):
    #Add species names, label focal taxa, outgroup taxa
    for node in t.traverse():
        if node.is_leaf():
            ID = node.name.split("|")[0:1][0]
            if ID in species_name_dict.keys():
                new_name = node.name + "@" + species_name_dict[ID]
                if node.name in focal_taxa:
                    new_name = "FOCAL__" + new_name
                name_face = TextFace(new_name, ftype='helveticae', fsize=8, fgcolor="black", bold=True)
                group_face = TextFace("@" + group_dict[ID][0], ftype='helveticae', fsize=8, fgcolor=group_dict[ID][1], bold=True)
                node.add_face(name_face, column=0, position='branch-right')
                node.add_face(group_face, column=1, position='branch-right')
            elif focal_taxonomy in node.name:
                new_name = node.name
                if node.name in focal_taxa:
                    new_name = "FOCAL__" + new_name
                name_face = TextFace(new_name, ftype='helveticae', fsize=8, fgcolor="black", bold=True)
                group_face = TextFace("@" + focal_taxonomy, ftype='helveticae', fsize=8, fgcolor="#636363", bold=True)
                node.add_face(name_face, column=0, position='branch-right')
                node.add_face(group_face, column=1, position='branch-right')
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

#input files
map = snakemake.input.map
info = snakemake.input.info
tree = snakemake.input.tree
taxa = snakemake.input.taxa

#output files
overview = snakemake.output.out
out = snakemake.output.out
seqs = snakemake.output.seqs
focal = snakemake.output.focal
outgroup = snakemake.output.outgroup
removed = snakemake.output.removed
summary = snakemake.output.summary

#parameters
cluster_name = snakemake.params.cluster_name
num = snakemake.params.num
min_taxa = int(snakemake.params.min_taxa)
max_taxa = int(snakemake.params.max_taxa)
focal_taxonomy = snakemake.params.focal_taxonomy

#Parse mapping file
species_name_dict, group_dict, focal_list = parse_mapping_file(map)

#Extract annotation
annotation = get_cluster_annotation(info, cluster_name)

#Make dictionary of focal_node taxa and nodes
focal_taxa_dict = make_focal_taxa_dict(taxa)

#Make output files
#If the number of taxa in the tree is less than or equal to max_taxa output all
try:
    t = initial_root_tree(tree, focal_taxa_dict)
except:
    open(out,'w')
    open(focal,'w')
    open(seqs,'w')
    open(outgroup,'w')
    open(removed,'w')
    open(summary, 'w')
else:
    t, focal_node_dict = make_focal_node_dict(t, focal_taxa_dict)
    if len(t) <= max_taxa:
        if int(num) == 0:
            focal_taxa = [item for sublist in list(focal_taxa_dict.values()) for item in sublist]
            non_focal_taxa = np.setdiff1d(t.get_leaf_names(), focal_taxa)
            t = root_tree(t, focal_node_dict["0"])
            t = set_tree_style(t)
            t = label_node(t, focal_node_dict.values())
            t = label_tree_taxa_simple(t, group_dict, species_name_dict, focal_taxa, focal_taxonomy)
            output_tree(t, out, annotation)
            with open(focal, "w") as f:
                for taxon in focal_taxa:
                    f.write("%s\n" %taxon)
            with open(seqs, "w") as f:
                for taxon in non_focal_taxa:
                    f.write("%s\n" %taxon)
            open(outgroup,'w')
            open(removed,'w')
            #Print subtree characteristics
            print("SMALL TREE - ALL TAXA INCLUDED")
            print("Cluster name: %s" %(cluster_name))
            print("Number of taxa in tree: %s" %(len(t)))
            print("Number of focal taxa in tree: %s" %(len(focal_taxa)))
            with open(summary, "w") as f:
                f.write("%s\t%s\t%s\t%s\t%s\t0\t0" %(cluster_name, num, len(t), len(t), len(focal_taxa)))
        else:
            open(out,'w')
            open(focal,'w')
            open(seqs,'w')
            open(outgroup,'w')
            open(removed,'w')
            open(summary, 'w')
    else:
        try:
            len(focal_node_dict[num])
        except:
            open(out,'w')
            open(focal,'w')
            open(seqs,'w')
            open(outgroup,'w')
            open(removed,'w')
            open(summary, 'w')
        else:
            #Make list of focal taxa in node of interest
            focal_taxa = focal_taxa_dict[num]
            #Select subtree and subtree taxa
            t = root_tree(t, focal_node_dict[num])
            subtree_candidates = subtree_candidate_selection(t, focal_node_dict[num])
            refined_subtree_candidates = refine_subtree_candidates(subtree_candidates, focal_taxa, max_taxa, min_taxa)
            subtree, subtree_taxa = select_subtree(refined_subtree_candidates)
            #Get outgroup seqs
            outgroup_taxa = get_outgroup_seqs(subtree, subtree_taxa)
            #For large subtrees, remove well-supported clades away from taxa of interest
            removed_taxa = []
            if len(subtree) > (max_taxa + len(focal_taxa)):
                distant_nodes = find_distant_nodes(subtree, subtree_taxa, focal_node_dict[num])
                unique_distant_nodes = unique_distant_nodes(distant_nodes)
                removed_taxa = remove_distant_nodes(subtree_taxa, unique_distant_nodes, max_taxa, removed_taxa, focal_taxa)
            subtree_taxa = np.setdiff1d(subtree_taxa, removed_taxa)
            #Get tree with labels and output
            t = set_tree_style(t)
            t = label_node(t, focal_node_dict.values())
            t = label_subtree(t, subtree)
            t = label_refined_candidates(t, subtree, refined_subtree_candidates)
            t = label_tree_taxa(t, group_dict, species_name_dict, focal_taxa, focal_taxonomy, outgroup_taxa, removed_taxa)
            #output tree pdf
            output_tree(t, out, annotation)
            #Output taxa from selected subtree (excluding focal), focal taxa, and outgroup taxa
            with open(focal, "w") as f:
                for taxon in focal_taxa:
                    f.write("%s\n" %taxon)
            with open(seqs, "w") as f:
                for taxon in subtree_taxa:
                    if taxon not in focal_taxa:
                        f.write("%s\n" %taxon)
            with open(outgroup, "w") as f:
                for taxon in outgroup_taxa:
                    f.write("%s\n" %taxon)
            if len(removed_taxa) > 0:
                with open(removed, "w") as f:
                    for taxon in removed_taxa:
                        f.write("%s\n" %taxon)
            else:
                open(removed,'w')
            #Print subtree characteristics
            print("LARGE TREE - SELECTION OF SUBTREE")
            print("Cluster name: %s" %(cluster_name))
            print("Number of taxa in initial tree: %s" %(len(t)))
            print("Number of non-focal taxa in subtree: %s" %(len(subtree_taxa) - len(focal_taxa)))
            print("Number of taxa removed from subtree: %s" %(len(removed_taxa)))
            print("Number of outgroup taxa: %s" %(len(outgroup_taxa)))
            print("Number of focal taxa in subtree: %s" %(len(focal_taxa)))
            #Write summary
            with open(summary, "w") as f:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(cluster_name, num, len(t), len(subtree), len(focal_taxa), len(outgroup_taxa), len(removed_taxa)))
