#!/bin/bash

tree_path=$1
out_path=$2
workdir=$3
snakedir=$4
map_file=$5
info_file=$6
focal_taxonomy=$7

echo "Tree path:" $tree_path
echo "Output path:"$out_path
echo "Working dir:" $workdir
echo "Snakefile dir:"$snakedir
echo "Mapping file:" $map_file
echo "Info file:" $info_file
echo "Focal taxonomy:" $focal_taxonomy

#Get old and new names
for tree in  "$tree_path"/*.treefile
do
  name=$(echo $tree | cut -d "/" -f2 | cut -d "." -f1-2)
  echo $name
  grep ">" "$workdir"/7.subtree/"$name".subtree.align.trim | cut -d ">" -f2 > "$tree_path"/"$name".trim.list
  grep ">" "$workdir"/7.subtree/"$name".subtree.align.trim | cut -d ">" -f2 | sed "s/|/_/g" | sed "s/\[/_/g" | sed "s/\]/_/g" | sed "s/\+/_/g" | sed "s/'/_/g" | sed 's/@/_/g' > "$tree_path"/"$name".tree.list
done

#get renamed tree
for tree in  "$tree_path"/*.treefile
do
   echo $tree
   name=$(echo $tree | cut -d "/" -f2 | cut -d "." -f1-2)
   python "$snakedir"/scripts/replace_names.py \
   -f $tree \
   -o "$tree_path"/"$name".tree.list  \
   -n "$tree_path"/"$name".trim.list \
   -out "$out_path"/"$name".NH.tree
 done

#Colour subtrees
for tree in "$tree_path"/*.NH.tree
do
  echo $tree ;
  ID=$(echo $tree | rev | cut -d "/" -f1 | rev | cut -d "." -f1-2)
  echo $ID
  echo "$workdir"/6.subtree_faa/"$ID".subtree.faa
  python "$snakedir"/scripts/colour_subtree_cmdl.py \
    -tree $tree \
    -map $map_file \
    -info $info_file \
    -phyla "$workdir"/7.subtree/all_subset_phyla.list \
    -domains "$workdir"/8.subtree_interproscan/"$ID".subtree.interproscan.tsv \
    -faa "$workdir"/6.subtree_faa/"$ID".subtree.faa \
    -focal_taxonomy $focal_taxonomy \
    -out "$out_path"/"$ID".pdf
done

#Get out sister and nested taxonomy
for tree in "$tree_path"/*.NH.tree
do
  ID=$(echo $tree | rev | cut -d "/" -f1 | rev | cut -d "." -f1-2)
  python "$snakedir"/scripts/find_sister_nested_taxonomy_cmdl.py \
    -tree $tree \
    -map $map_file \
    -tax $focal_taxonomy \
    -out "$out_path"/"$ID".sister_nested_taxonomy.tsv \
    -cluster "$ID"
done

#Combine sister and nested taxonomy
cat "$out_path"/*.sister_nested_taxonomy.tsv > "$out_path"/ALL.sister_nested_taxonomy.tsv

#Make plots of taxonomy
Rscript "$snakedir"/scripts/plot_HGT_sources_cmdl.R \
  -in  "$out_path"/ALL.sister_nested_taxonomy.tsv \
  -out "$out_path"/ALL.sister_nested_taxonomy.pdf

#for pathway in * ; do echo $pathway ; ../../workflow/run_labelling_HGT_sources.sh $pathway $pathway ../$pathway ../../workflow ../../mapping_files/species_name_map.tsv ../$pathway/$pathway.tsv "Chlamydiae" ; done
