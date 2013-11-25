##
#
# FILE: TEST.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains test code for all components in
# phylogenetic profiling project.
#

from imports import *

import ParseFASTA
import Matching
import NormalLookup
import BLAST
import GUI
import Profiling
import Phylobuilder


## DEFINITIONS

TEST_FILE_1 = "../datasets/Saccharomyces_cerevisiae.EF4.73.pep.all.fa"
TEST_FILE_1_SPECIES = "Saccharomyces Cerevisiae"


## FUNCTIONS

# Test GUI functionality
def test_gui():

    print("rig")

    return


# Test loading and parsing a FASTA file
def test_load_parse():

    print("######################################")
    print("Testing FASTA file LOAD and PARSE functionality...")
    print("######################################")
    print()

    print("Filename: ", TEST_FILE_1)
    print("Species: ", TEST_FILE_1_SPECIES)
    print()
    
    file_reader = FASTAFile(TEST_FILE_1, TEST_FILE_1_SPECIES)

    num_genes = file_reader.get_num_unique_genes()
    print("Number of unique genes in ", file_reader.species, ": ", num_genes)
    print()

    gene_IDs = file_reader.get_unique_gene_IDs()
    print("Here are ten gene IDs from ", file_reader.species, ":")
    for index in range(10):
        print(gene_IDs[index])
    print()

    print("Here are the variants of the first three gene IDs:")
    for index in range(3):
        variants = file_reader.get_variants(gene_IDs[index])
        print("Gene ID: ", gene_IDs[index])
        for sub_index in range(len(variants)):
            print("Variant ID", sub_index, ": ", variants[sub_index])
    print()

    print()
    print("######################################")
    print("DONE testing FASTA file LOAD and PARSE functionality.")
    print("######################################")

    return


# Test matching objects
def test_matching():

    print("rig")

    return


# Test phylobuilder main
def test_builder():

    print("rig")

    return


# Run all tests
def test_project():

    test_load_parse()

    return


if __name__ == '__main__':
    test_project()
