##
#
# FILE: TEST.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains test code for all components in
# phylogenetic profiling project.
#

import time
import logging

from imports import *

from ParseFASTA import *
from Matching import *
from GUI import *
from Profiling import *
from PhyloBuilder import *
from BLAST import *


## DEFINITIONS

TEST_FILE_1 = "../datasets/Saccharomyces_cerevisiae.EF4.73.pep.all.fa"
TEST_FILE_1_SPECIES = "Saccharomyces Cerevisiae"

TEST_FILE_2 = "../datasets/Ciona_intestinalis.KH.73.pep.all.fa"
TEST_FILE_2_SPECIES = "Ciona Intestinalis"

LOG_FILE_PREFIX = "TEST_"
LOG_FILE_SUFFIX = ".txt"

LOGGER_NAME = "phylo_logger"


## FUNCTIONS

# Test GUI functionality
def test_gui(logger):

    logger.debug("rig")

    return


# Test loading and parsing a FASTA file
def test_load_parse(logger, filename, species):

    logger.debug("BEGIN testing FASTA file LOAD and PARSE functionality...")

    logger.debug("Filename: %s", filename)
    logger.debug("Species: %s", species)
    
    file_reader = FASTAFile(filename, species)

    num_genes = file_reader.get_num_unique_genes()
    logger.debug("Number of unique genes in %s: %d", file_reader.species, num_genes)

    gene_IDs = file_reader.get_unique_gene_IDs()
    logger.debug("Here are ten gene IDs from %s:", file_reader.species)
    for index in range(10):
        logger.debug("%s", gene_IDs[index])

    logger.debug("Here are the variants of the first three gene IDs:")
    for index in range(3):
        variants = file_reader.get_variants(gene_IDs[index])
        logger.debug("Gene ID: %s:", gene_IDs[index])
        for sub_index in range(len(variants)):
            variant_ID = variants[sub_index].variant
            logger.debug("Variant ID %d: %s", sub_index, variant_ID)

    logger.debug("DONE testing FASTA file LOAD and PARSE functionality.")
    logger.debug("")

    return file_reader


# Test matching objects
def test_matching(logger, file_list):

    logger.debug("BEGIN testing matching object creation from FASTA record...")

    tree = PhyloTree(file_list)

    logger.debug("Species in Tree:")
    for index in range(tree.num_species):
        logger.debug("Species %d: %s", index, tree.species_names[index])

    logger.debug("Here are the first 5 genes in each species and their variants.")
    for species_index in range(tree.num_species):
        logger.debug("Species %d: %s", species_index, tree.species[species_index])
        for gene_index in range(tree.species[species_index].num_genes):
            logger.debug("  Gene %d: %s", gene_index, tree.species[species_index].genes[gene_index].gene_ID)
            for variant_index in range(tree.species[species_index].genes[gene_index].num_variants):
                logger.debug("    Variant %d: %s", variant_index, tree.species[species_index].genes[gene_index].variant_IDs[variant_index])

    logger.debug("Here is all the information on a sample variant.")
    first_variant = tree.species[0].genes[0].variants[0]
    logger.debug("Variant ID: %s", first_variant.variant_ID)
    logger.debug("Peptide Type: %s", first_variant.peptide_type)
    logger.debug("Chromosome: %s", first_variant.chromosome)
    logger.debug("Sequence: %s", first_variant.sequence)
    logger.debug("Parent Gene ID: %s", first_variant.gene_ID)

    logger.debug("DONE testing matching object creation from FASTA record.")

    return


# Test phylobuilder main
def test_builder(logger):

    logger.debug("rig")

    return


# Returns a logger object
def get_logger():

    time_date_string = time.strftime("%d-%m-%y_%H:%M:%S")
    log_filename = LOG_FILE_PREFIX + time_date_string + LOG_FILE_SUFFIX

    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)

    handler = logging.FileHandler(log_filename)
    logger.addHandler(handler)

    logger.debug("Created logger object.")
    logger.debug("")

    return logger

# Run all tests
def test_project():

    # Get the logger object for all testing
    logger = get_logger()

    # Test file parsing and record creation
    records_list = []

    file_records_1 = test_load_parse(logger, TEST_FILE_1, TEST_FILE_1_SPECIES)
    file_records_2 = test_load_parse(logger, TEST_FILE_2, TEST_FILE_2_SPECIES)

    records_list.append(file_records_1)
    records_list.append(file_records_2)

    test_matching(logger, records_list)

    return


if __name__ == '__main__':
    test_project()
