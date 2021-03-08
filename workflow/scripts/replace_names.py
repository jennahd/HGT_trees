#!/usr/bin/env python

r"""
Author: Jennah Dharamshi
Date: 170712

This script takes an input file with original names, an input file with
replacement names and a file in which you want to make the name replacement.

Names to be replaced and replacement name must be on the same line number
(one name per line) in the respective file. The script will replace whatever is
on the corresponding line number with the corresponding name.

The script only replaces strings.
"""

__author__ = 'djennah'

################################################################################
#Modules
import argparse
import os

################################################################################
#Functions
def load_name_lists(old_names, new_names):
    old_names_list = []
    new_names_list = []
    handle_old = open(old_names, 'rU')
    handle_new = open(new_names, 'rU')
    for name in handle_old:
        item = name.strip("\n")
        old_names_list.append(item)
    for name in handle_new:
        item = name.strip("\n")
        new_names_list.append(item)
    return old_names_list, new_names_list
    handle_new.close()
    handle_old.close()

def replace_old(old_names_list, new_names_list, old_file):
    name_dict = {}
    item_number = 0
    while item_number < len(old_names_list):
        name_dict[old_names_list[item_number]] = [new_names_list[item_number]]
        item_number += 1
    infile = open(old_file, 'rU').read()
    for old, new in name_dict.items():
        infile = infile.replace(old, new[0])
    return infile

################################################################################
#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='replace_names.py',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Replaces strings in a given file with new strings.')

parser.add_argument('-f', '--file_replace_names', required=True, help = 'Path \
	to the input file where strings (names) will be replaced.')

parser.add_argument('-o', '--old_names', required=True, help = 'Path to the \
	input list of old strings in the original file to be replaced.')

parser.add_argument('-n', '--new_names', required=False, help = 'Path to the \
	input list of new strings to replace old strings for the new output file.')

parser.add_argument('-out', '--output_file_replaced_names', required=True, help = \
    'Path and name for the output file to be created, which will contain the new \
    names that were replaced from the original file.')

args = parser.parse_args()

#Check files and folders
if not os.path.exists(args.file_replace_names):
	raise IOError("File does not exist at given path")
if not os.path.exists(args.old_names):
	raise IOError("File does not exist at given path")
if not os.path.exists(args.new_names):
	raise IOError("File does not exist at given path")
out_folder = (os.path.abspath(args.output_file_replaced_names)).rsplit("/",1)[0]
if not os.path.exists(out_folder):
    raise IOError("Directory to write output file does not exist")
if os.path.exists(args.output_file_replaced_names):
    print("This file already exists, file will be overwritten")

################################################################################
#Implementation

if __name__ == "__main__":
    #Load name lists
    #Parse list of old and new names
    old_names_list, new_names_list = load_name_lists(args.old_names, args.new_names)
    #Check replacements
    if len(old_names_list) == len(new_names_list):
        print("Number of sequences to be replaced %s" % len(old_names_list))
        item_number = 0
        while item_number < len(old_names_list):
            print("%s will be replaced with %s" % \
            (old_names_list[item_number], new_names_list[item_number]))
            item_number += 1
    else:
        raise ValueError("There are a different number of items to be replaced \
        in the old_names list and replacements in the new_names list.")
    #Make replacements
    new_file = replace_old(old_names_list, new_names_list, args.file_replace_names)
    #output new file
    out = open(args.output_file_replaced_names,'w')
    out.write(new_file)
    out.close
