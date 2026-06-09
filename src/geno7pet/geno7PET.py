# MIT License
#
# Copyright (c) 2025 Avril Coghlan
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#====================================================================#
# VERSION 1.0.1, 14th-Oct-2025.
#====================================================================#
import shutil
import os
import sys
import random
from collections import defaultdict
from importlib import resources
from pathlib import Path
import argparse


#====================================================================#
# define a function to make a file name for a temporary file

def make_filename(outputdir):

    found_name = 0
    filename = "none"

    while found_name == 0:
        random_number = "tmp%f" % random.random()
        poss_filename = os.path.join(outputdir, random_number);
        if not os.path.exists(poss_filename):
            filename = poss_filename
            found_name = 1

    assert(found_name == 1 and filename != 'none')

    return filename

#====================================================================#
# define subroutine to check if a SNP is found in a particular isolate's assembly, based on blast results:

def check_if_SNP_in_isolate(snp, alleles_found_set, input_assembly_file, clade, give_verbose_output):

    # now check if the alternative allele is found:
    alternative_allele_name = "%s_alternative" % (snp)
    if alternative_allele_name in alleles_found_set:
        is_alternative_allele_found = True
    else:
        is_alternative_allele_found = False

    # now check if the wildtype allele is found:
    wildtype_allele_name = "%s_wildtype" % (snp)
    if wildtype_allele_name in alleles_found_set:
        is_wildtype_allele_found = True
    else:
        is_wildtype_allele_found = False

    if give_verbose_output == 'yes':
        print("clade=",clade,"snp=",snp,"alt",is_alternative_allele_found,"wt",is_wildtype_allele_found)

    return is_alternative_allele_found

#====================================================================#
# define subroutine to read in the SNPs to use to classify in a particular clade:

def read_input_snps_to_use_for_classifying_clade(file_with_snps_for_clade, clade, classifications_for_isolate, input_assembly_file, alleles_found_set, fst_cutoff_for_clade, give_verbose_output):

    # read in the input file of SNPs to use to classify in the clade 'clade':
    fileObj = open(file_with_snps_for_clade)
    total_number_snps_for_clade = 0
    number_snps_for_clade_found = 0
    for line in fileObj:
        # "AE003852_760185" 0.948913402470011
        # "AE003852_1709267" 0.949646833368501
        line = line.rstrip()
        temp = line.split()
        if len(temp) == 2: # the first line of the file has just "FST"
           snp = temp[0] # e.g. "AE003852_760185"
           fst = float(temp[1]) # e.g. 0.948913402470011
           temp2 = snp.split('\"')
           snp = temp2[1] # e.g. AE003852_760185
           # check if the FST score is >= fst_cutoff_for_clade:
           if fst >= fst_cutoff_for_clade:
               total_number_snps_for_clade = total_number_snps_for_clade + 1
               # check if the SNP was found in the isolate:
               is_snp_found = check_if_SNP_in_isolate(snp, alleles_found_set, input_assembly_file, clade, give_verbose_output)
               if is_snp_found:
                   format_string = "clade %s - FOUND snp %s FST %f" % (clade, snp, fst)
                   number_snps_for_clade_found = number_snps_for_clade_found + 1
               else:
                   format_string = "clade %s - missing snp %s FST %f" % (clade, snp, fst)
    fileObj.close()

    # calculate the number of SNPs for this clade that were found:
    if total_number_snps_for_clade > 0:
        percent_snps_for_clade_found = (number_snps_for_clade_found * 100)/total_number_snps_for_clade
        format_string = "FINDING SNPS:__________clade %s - percent SNPs found %f" % (clade, percent_snps_for_clade_found)
        if give_verbose_output == 'yes':
            print(format_string)
        if percent_snps_for_clade_found >= 50.0:
            clade2 = "%s_%f" % (clade,percent_snps_for_clade_found)
            classifications_for_isolate.append(clade2)
    else:
        format_string = "FINDING SNPS:__________clade %s - total_number_snps_for_clade %d number_snps_for_clade_found %d" % (clade, total_number_snps_for_clade, number_snps_for_clade_found)

    return classifications_for_isolate

#====================================================================#
# define subroutine to check siblings are not in the classifications set:

def check_siblings_are_not_in_classifications(siblings, current_classifications_set):
    """check siblings are not in the classifications set:
    >>> check_siblings_are_not_in_classifications(['3.1.1.1', '3.1.1.3', '3.1.1.4', '3.1.1.5'], {'3.1.1.2.1','3.1.1.2','3.1.1','3.1','1.1'})
    'fine'
    >>> check_siblings_are_not_in_classifications(['3.1.1.1', '3.1.1.2', '3.1.1.4', '3.1.1.5'], {'3.1.1.2.1','3.1.1.2','3.1.1','3.1','1.1'})
    'warning'
    """

    for sibling in siblings: # for each of the siblings:
        if sibling in current_classifications_set:
            return('warning')

    return('fine')

#====================================================================#
# define a subroutine to remove classifications:

def remove_classifications(mylist, current_classifications_set):
    """ remove classifications from the classifications set:
    >>> remove_classifications(['3.1.1.4'], {'3.1.1.2.1', '3.1.1.2', '3.1.1.4', '3.1.1', '3.1', '1.1'})
    ['1.1', '3.1', '3.1.1', '3.1.1.2', '3.1.1.2.1']
    """

    for classification in mylist:
        if classification in current_classifications_set:
            current_classifications_set.remove(classification)

    current_classifications_list = sorted(list(current_classifications_set))

    return current_classifications_list

#====================================================================#
# define a subroutine to remove the descendants of the siblings of a node:

def remove_descendants_of_siblings(siblings, current_classifications_set, parents_dict, children_dict):
    """remove descendants of the siblings of a node:
    >>> remove_descendants_of_siblings(['3.1.0', '3.1.1', '3.1.4', '3.1.5'], {'1.1','3.1','3.1.2','3.1.2.2','3.1.1.3'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    ['1.1', '3.1', '3.1.2', '3.1.2.2']
    """

    for sibling in siblings:
        # get all the descendants of sibling:
        descendants = get_descendants(sibling, parents_dict, children_dict)
        current_classifications_list = remove_classifications(descendants, current_classifications_set)
        current_classifications_set = set(current_classifications_list)

    current_classifications_list = sorted(list(current_classifications_set))

    return current_classifications_list

#====================================================================#
# define subroutine to get the sibling nodes of node 'node' in the classification tree:

def get_siblings(node, parents_dict, children_dict):
    """get sibling nodes of a node in the classication tree
    >>> get_siblings('3.1.3', {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    ['3.1.0', '3.1.1', '3.1.2', '3.1.4', '3.1.5']
    >>> get_siblings('1.1', {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    ['1.0']
    >>> get_siblings('3.1.2.1', {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    ['3.1.2.0', '3.1.2.2']
    """

    # define a list variable for the siblings of node:
    siblings_list = list()

    # find the parent of node:
    assert(node in parents_dict)
    parent = parents_dict[node]

    # find the children of the parent:
    assert(parent in children_dict)
    childrenstring = children_dict[parent]
    children = childrenstring.split(',')
    for child in children:
        if child != node:
            siblings_list.append(child)

    return siblings_list

#====================================================================#
# define a subroutine to get the descendants of a node:

def get_descendants(mynode, parents_dict, children_dict):
    """get the descendants of a node
    >>> get_descendants('3.1.2', {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    ['3.1.2.0', '3.1.2.1', '3.1.2.2', '3.1.2.1.0', '3.1.2.1.1', '3.1.2.1.1.0', '3.1.2.1.1.1']
    """

    # define a list variable for the descendants of a node:
    descendants = list()

    # now traverse the tree to find the descendants of node 'mynode':
    # use a breadth-first search of the tree:
    # start at node 'mynode':
    myqueue = []
    myresult = {}
    myqueue.append( (mynode, 0) )
    while myqueue:
        mynode, mydist = myqueue.pop(0) # Take the first item off the list 'myqueue' (and remove it)
        myresult[mynode] = mydist
        if mynode in children_dict: # If the mynode *has* children
            childrenstring = children_dict[mynode]
            children = childrenstring.split(',')
            for child in children:
                # Get the first member of each tuple, see
                # http://stackoverflow.com/questions/12142133/how-to-get-first-element-in-a-list-of-tuples
                myqueue_members = [x[0] for x in myqueue]
                if child not in myresult and child not in myqueue_members: # Don't visit a second time
                    myqueue.append( (child, mydist+1) )
                # Add node 'child' to the descendants list:
                descendants.append(child)

    return descendants

#====================================================================#
# define subroutine to find the most precise classification:

def find_most_precise_classification(classifications_set):
    """find most precise classification:
    >>> find_most_precise_classification({'3.1.1.1', '3.1.1', '3.1', '1.1'})
    '3.1.1.1'
    >>> find_most_precise_classification({'3.1.1.1', '3.1.1', '3.1', '1.1', '3.1.1.1.2'})
    '3.1.1.1.2'
    """

    max_num_parts = 0
    most_precise_classification = None
    for classification in classifications_set:
        temp = classification.split('.')
        num_parts = len(temp)
        if num_parts > max_num_parts:
            most_precise_classification = classification
            max_num_parts = num_parts

    # check that there is just one classification with max_num_parts, i.e. we only have one most precise classification:
    num_classifications_with_max_num_parts = 0
    for classification in classifications_set:
        temp = classification.split('.')
        num_parts = len(temp)
        if num_parts == max_num_parts:
            num_classifications_with_max_num_parts = num_classifications_with_max_num_parts + 1
    assert(num_classifications_with_max_num_parts == 1)

    return most_precise_classification

#====================================================================#
# define subroutine to make sure that at each level, we have some member of the classification tree selected:

def check_have_classification_at_each_level(classifications_set, parents_dict):
    """ make sure that at each level we have some member of the tree selected.
    >>> check_have_classification_at_each_level({'3.1.1.1', '3.1.1', '3.1', '1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    'fine'
    >>> check_have_classification_at_each_level({'3.1.1.1', '3.1', '1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    'warning'
    """

    result = 'fine'

    # now traverse the tree to make sure the classifications make sense:
    # start at the most precise classification, and then work upwards towards the root node:
    most_precise_classification = find_most_precise_classification(classifications_set)
    node = most_precise_classification
    while node != 'root':
        assert(node in parents_dict)
        parent = parents_dict[node]
        if parent not in classifications_set:
            if parent != 'root':
                result = 'warning'
        node = parent

    return result

#====================================================================#
# find the child node of a particular node, that is a 'catchall' node (ending in .0):

def get_catchall_child_node(node, children_dict):
    """Find the child node of a particular node, that is a catchall node ending in .0:
    >>> get_catchall_child_node('3.1', {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    '3.1.0'
    >>> get_catchall_child_node('1.1', {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    '2.0'
    >>> get_catchall_child_node('root', {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'})
    '1.0'
    """

    catchall_child_node = None
    assert(node in children_dict)
    childrenstring = children_dict[node]
    children = childrenstring.split(',')
    for child in children:
        temp = child.split('.')
        last_number = temp[-1]
        if last_number == '0':
            assert(catchall_child_node is None)
            catchall_child_node = child
    assert(catchall_child_node is not None)

    return catchall_child_node

#====================================================================#
# now traverse the tree in a breadth-first search from the root downwards:
# if there is no classification at a certain level, then classify to .0 at that level

def assign_to_catchallnodes(current_classifications_set, children_dict, parents_dict):
    """traverse the tree in a breadth-first search from the root downwards, and if there is no classification at a certain level then classify to a .0 node at that level:
    >>> assign_to_catchallnodes({'1.1', '3.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    ['1.1', '3.1', '3.1.0']
    >>> assign_to_catchallnodes(set(), {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    ['1.0']
    >>> assign_to_catchallnodes({'1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    ['1.1', '2.0']
    """

    queue = []
    result_dict = {}
    queue.append( ('root', 0) )
    while queue:
        node, dist = queue.pop(0) # Take the first item off the list 'queue' (and remove it)
        result_dict[node] = dist
        if node in children_dict: # If the node *has* children
            childrenstring = children_dict[node]
            children = childrenstring.split(',')
            # check if at least one child of node 'node' is in current_classifications_set
            have_child_in_classifications_set = False
            for child in children:
                # Get the first member of each tuple, see
                # http://stackoverflow.com/questions/12142133/how-to-get-first-element-in-a-list-of-tuples
                queue_members = [x[0] for x in queue]
                if child not in result_dict and child not in queue_members: # Don't visit a second time
                    queue.append( (child, dist+1) )
                    # If node 'child' is in the current classification set:
                    if child in current_classifications_set:
                        have_child_in_classifications_set = True
            # Check if the node 'node' is in the current classifications set (or is the root node) and
            # we don't have any children of the node 'node' in the current classifications set:
            if node in current_classifications_set or node == 'root':
                if not have_child_in_classifications_set:
                    # need to figure out what is the child node of 'node' that is a 'catchall' node:
                    catchall_child_node = get_catchall_child_node(node, children_dict)
                    current_classifications_set.add(catchall_child_node)

    current_classifications_list = sorted(list(current_classifications_set))

    return current_classifications_list

#====================================================================#
# define subroutine to traverse the tree in a breadth-first search from the root downwards, and
# (i) check that if a node is in the classification set, none of its siblings are in the classification set (if they are, given an error)
# (ii) if a node is in the classification set, remove descendants of its siblings from the classification set

def prune_classification_set(current_classifications_set, children_dict, parents_dict):
    """ traverse the tree and check that if a node is in the classication set then its sibilings are not, and also remove descendant nodes of its siblings
    >>> prune_classification_set({'3.1.1.1', '3.1.1', '3.1', '1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    ('fine', ['1.1', '3.1', '3.1.1', '3.1.1.1'])
    >>> prune_classification_set({'3.1.1.1', '3.1.1', '3.1.2', '3.1', '1.1'}, {'root': '1.0,1.1', '1.1': '2.0,3.1', '3.1': '3.1.0,3.1.1,3.1.2,3.1.3,3.1.4,3.1.5', '3.1.1': '3.1.1.0,3.1.1.1,3.1.1.2,3.1.1.3,3.1.1.4,3.1.1.5', '3.1.2': '3.1.2.0,3.1.2.1,3.1.2.2', '3.1.1.2': '3.1.1.2.0,3.1.1.2.1', '3.1.2.1': '3.1.2.1.0,3.1.2.1.1', '3.1.2.1.1': '3.1.2.1.1.0,3.1.2.1.1.1'}, {'1.0': 'root', '1.1': 'root', '2.0': '1.1', '3.1': '1.1', '3.1.0': '3.1', '3.1.1': '3.1', '3.1.2': '3.1', '3.1.3': '3.1', '3.1.4': '3.1', '3.1.5': '3.1', '3.1.1.0': '3.1.1', '3.1.1.1': '3.1.1', '3.1.1.2': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.3': '3.1.1', '3.1.1.4': '3.1.1', '3.1.1.5': '3.1.1', '3.1.2.0': '3.1.2', '3.1.2.1': '3.1.2', '3.1.2.2': '3.1.2', '3.1.1.2.0': '3.1.1.2', '3.1.1.2.1': '3.1.1.2', '3.1.2.1.0': '3.1.2.1', '3.1.2.1.1': '3.1.2.1', '3.1.2.1.1.0': '3.1.2.1.1', '3.1.2.1.1.1': '3.1.2.1.1'})
    ('warning', ['1.1', '3.1', '3.1.1', '3.1.1.1', '3.1.2'])
    """

    result = 'fine'

    queue = []
    result_dict = {}
    queue.append( ('root', 0) )
    while queue:
        node, dist = queue.pop(0) # Take the first item off the list 'queue' (and remove it)
        result_dict[node] = dist
        if node in children_dict: # If the node *has* children
            childrenstring = children_dict[node]
            children = childrenstring.split(',')
            for child in children:
                # Get the first member of each tuple, see
                # http://stackoverflow.com/questions/12142133/how-to-get-first-element-in-a-list-of-tuples
                queue_members = [x[0] for x in queue]
                if child not in result_dict and child not in queue_members: # Don't visit a second time
                    queue.append( (child, dist+1) )
                    #----------------------------------------------------------------#
                    # If node 'child' is in the current classification set:
                    if child in current_classifications_set:
                        # Check the sibling nodes of node 'child' are not in the classifications:
                        siblings = get_siblings(child, parents_dict, children_dict)
                        result = check_siblings_are_not_in_classifications(siblings, current_classifications_set)
                        if result == 'warning':
                            current_classifications_list = sorted(list(current_classifications_set))
                            return(result, current_classifications_list)
                        # Remove the descendants of the siblings of node 'child' from the set of classifications:
                        current_classifications_list = remove_descendants_of_siblings(siblings, current_classifications_set, parents_dict, children_dict)
                        current_classifications_set = set(current_classifications_list)
                        #----------------------------------------------------------------#

    current_classifications_list = sorted(list(current_classifications_set))

    return(result, current_classifications_list)

#====================================================================#
# define subroutine to filter the classifications for the isolate to make sure they make sense in terms of the tree of clades:

def compare_classifications_to_tree(initial_classifications_for_isolate, input_assembly_file, parents_dict, children_dict, give_verbose_output):

    # initial_classifications_for_isolate are the classifications made based on which SNPs are present. They
    # don't necessarily match the classification tree. We might have classifications that disagree with each other
    # e.g. 2and3, 3, 3.1, 3.1.3, 3.1.3.1, 3.1.2.1

    # make a set of the classifications that we have already:
    current_classifications_set = set()
    for classification in initial_classifications_for_isolate:
        temp = classification.split('_')
        classification1 = temp[1] # e.g. 2and3, 3, 3.2, 3.2.1, etc.
        current_classifications_set.add(classification1)

    # if there are not any classifications left in 'current_classifications_set', then add wave 1:
    if len(current_classifications_set) == 0:
        classification1 = 'wave1_1.0_0.0'
        current_classifications_set.add(classification1)
        result = 'fine'
    else:
        # now traverse the tree in a breadth-first search from the root downwards:
        # (i) check that if a node is in the classification set, none of its siblings are in the classification set (if they are, give an error)
        # (ii) if a node is in the classification set, remove descendants of its siblings from the classifications set
        (result, current_classifications_list) = prune_classification_set(current_classifications_set, children_dict, parents_dict)
        current_classifications_set = set(current_classifications_list)
        if result == 'warning':
            return(current_classifications_set, 'warning')

        # now travel up the tree from the most precise classification upwards to the root, checking we have a classification at each level:
        # note: we check that there is just one most precise classification (one with the most dots): if there is >1, given an error
        result = check_have_classification_at_each_level(current_classifications_set, parents_dict)
        if result == 'warning':
            return(current_classifications_set, 'warning')

        # now traverse the tree in a breadth-first search from the root downwards:
        # if there is no classification at a certain level, then classify to .0 at that level
        current_classifications_list = assign_to_catchallnodes(current_classifications_set, children_dict, parents_dict)
        current_classifications_set = set(current_classifications_list)

    # print out current set of classifications at this stage:
    classifications_string = ",".join(current_classifications_set)
    if give_verbose_output == 'yes':
        print("COMPARED TO TREE: classifications_set=",classifications_string)

    return(current_classifications_set, result)

#====================================================================#

# define subroutine to read in the file with files of SNPs to use for classifying:

def read_input_snps_to_use_for_classifying(input_snps_to_use_for_classifying_file, input_assembly_file, alleles_found_set, parents_dict, children_dict, output_file, give_verbose_output):

    # make a list to store the classifications for this isolate:
    initial_classifications_for_isolate = list()
    final_classifications_for_isolate = list()

    # read in the input file with files of SNPs:
    fileObj = open(input_snps_to_use_for_classifying_file)
    data_dir = Path(os.path.dirname(input_snps_to_use_for_classifying_file))
    for line in fileObj:
        # wave1 data_for_FST_statistics_update2_wave1_V2_FST.txt 0.98
        line = line.rstrip()
        temp = line.split()
        clade = temp[0]
        file_with_snps_for_clade = data_dir / temp[1]
        fst_cutoff_for_clade = float(temp[2])
        # read in the SNPs to use to classify in this clade:
        initial_classifications_for_isolate = read_input_snps_to_use_for_classifying_clade(file_with_snps_for_clade, clade, initial_classifications_for_isolate, input_assembly_file, alleles_found_set, fst_cutoff_for_clade, give_verbose_output)
    fileObj.close()

    # filter the classifications for the isolate to make sure they make sense in terms of the tree of clades:
    (final_classifications_for_isolate, result) = compare_classifications_to_tree(initial_classifications_for_isolate, input_assembly_file, parents_dict, children_dict, give_verbose_output)

    # write out the result to the output file:
    outputfileObj = open(output_file, "w")
    if result == 'fine':
        # print out the final classifications for this isolate:
        classifications_string = ",".join(final_classifications_for_isolate)
        format_string = "FINAL RESULT: input_assembly_file %s : classifications_string %s" % (input_assembly_file, classifications_string)
        if give_verbose_output == 'yes':
            print(format_string)
        # find the most precise classification:
        most_precise_classification = find_most_precise_classification(final_classifications_for_isolate)
        format_string = "%s %s\n" % (input_assembly_file, most_precise_classification)
        outputfileObj.write(format_string)
    elif result == 'warning':
        # print out the initial classifications for this isolate, and an error message:
        classifications_string = ",".join(initial_classifications_for_isolate)
        format_string = "WARNING UNCLASSIFIABLE: input_assembly_file %s : classifications_string %s" % (input_assembly_file, classifications_string)
        if give_verbose_output:
            print(format_string)
        most_precise_classification = "unclassifiable"
        format_string = "%s %s\n" % (input_assembly_file, most_precise_classification)
        outputfileObj.write(format_string)
    else:
        sys.stderr.write('ERROR:',result)
        sys.exit(1)
    outputfileObj.close()

    return

#====================================================================#
# make a blast database for the input assembly file, and run blast using the SNP alleles as queries:

def make_blast_db_and_run_blast(input_assembly_file, snp_wildtype_and_alternative_alleles_plus_flanking_seqs_fasta_file):


    # make a temporary copy of the input assembly file, so that we can make a blast database for it:
    tmp_dir = Path(make_filename(os.getcwd()))
    tmp_fasta = make_filename(tmp_dir)
    os.mkdir(tmp_dir)
    shutil.copy(input_assembly_file, tmp_fasta)
    # cmd0 = "cp %s %s" % (input_assembly_file, tmp_assembly_file)
    # os.system(cmd0)

    # make a blast database for the input assembly file:
    cmd1 = "makeblastdb -in %s -out %s -input_type fasta -dbtype nucl" % (tmp_fasta, tmp_fasta)
    os.system(cmd1)

    # make a temporary file to put the output of blast into:
    tmp_blast_output = make_filename(os.getcwd())

    # run blast using the SNP alleles as queries:
    cmd2 = "blastn -query %s -db %s -task blastn -word_size 20 -outfmt 6 > %s" % (snp_wildtype_and_alternative_alleles_plus_flanking_seqs_fasta_file, tmp_fasta, tmp_blast_output)
    os.system(cmd2)
    shutil.rmtree(tmp_dir)

    return tmp_blast_output


#====================================================================#
# read in the input classification tree:

def read_classification_tree(classification_tree):

    parents_dict = defaultdict()
    children_dict = defaultdict()

    fileObj = open(classification_tree)
    # parent children
    # root 1.0,1.1
    # ...
    for line in fileObj:
        line = line.rstrip()
        if not line.startswith('parent'):
             temp = line.split()
             parent = temp[0]
             childrenstring = temp[1]
             children = childrenstring.split(',')
             # store in dictionaries:
             for child in children:
                assert(child not in parents_dict)
                parents_dict[child] = parent
             assert(parent not in children_dict)
             children_dict[parent] = childrenstring
    fileObj.close()

    return(parents_dict, children_dict)

#====================================================================#
# read in the blast output file, to see which alleles were found in the assembly:

def read_blast_output(tmp_blast_output):

    # make a variable to record which alleles were found:
    alleles_found_set = set()

    fileObj = open(tmp_blast_output)
    for line in fileObj:
        # Note: The columns for blast m8 format are:
        # query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        # e.g.:
        # query1  1       99.010  101     1       0       1       101     25159   25259   8.30e-46        178
        # query2  1       100.000 101     0       0       1       101     25159   25259   1.95e-47        183
        line = line.rstrip()
        temp = line.split()
        snp_allele = temp[0] # e.g. AE003852_1417110_wildtype
        pc_id = float(temp[2]) # e.g. 99.010
        if pc_id == 100.000:
            alleles_found_set.add(snp_allele)
    fileObj.close()

    return alleles_found_set

#====================================================================#

def main():

    # check the command-line arguments:
    parser = argparse.ArgumentParser(
        prog="Geno7PET",
        description="geno7PET: A tool for classifying genotypes."
    )
    parser.add_argument(
        "input_assembly_file",
        help="input assembly file, e.g. my_fastas/AUSMDU00028282.fasta"
    )
    parser.add_argument(
        "output_file",
        help="name of the output file, e.g. AUSMDU00028282.fasta.out"
    )
    parser.add_argument(
        "-r",
        "--resources",
        help="directory with the parameter files, e.g. -r geno7PET_parameter_files",
        default=resources.files("geno7pet").joinpath("resources").joinpath("geno7PET_parameter_files")
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Select verbose output",
        action="store_true",
        default=False
    )
    args = parser.parse_args()

    input_assembly_file = args.input_assembly_file
    parameter_file_dir = args.resources
    output_file = args.output_file
    give_verbose_output = args.verbose

    # find the names of the parameter files:
    input_snps_to_use_for_classifying_file = "%s/%s" % (parameter_file_dir, "snps_for_classifying") # input file with the files of SNPs to use for classifying
    snp_wildtype_and_alternative_alleles_plus_flanking_seqs_fasta_file = "%s/%s" % (parameter_file_dir, "snp_wildtype_and_alternative_alleles_plus_flanking_seqs_50bp.fasta") # input fasta file with the wildtype and alternative alleles at a SNP site, including flanking sequences
    classification_tree = "%s/%s" % (parameter_file_dir, "classification_tree") # input classification tree

    # read in the input classification tree:
    (parents_dict, children_dict) = read_classification_tree(classification_tree)

    # make a blast database for the input assembly file, and run blast using the SNP alleles as queries:
    tmp_blast_output = make_blast_db_and_run_blast(input_assembly_file, snp_wildtype_and_alternative_alleles_plus_flanking_seqs_fasta_file)

    # read in the blast output file, to see which alleles were found in the assembly:
    alleles_found_set = read_blast_output(tmp_blast_output)

    # read in the file with files of SNPs to use for classifying:
    read_input_snps_to_use_for_classifying(input_snps_to_use_for_classifying_file, input_assembly_file, alleles_found_set, parents_dict, children_dict, output_file, give_verbose_output)

    # delete the temporary blast output file:
    os.unlink(tmp_blast_output)

    sys.stderr.write('FINISHED\n')

#====================================================================#

if __name__=="__main__":
    main()

#====================================================================#


