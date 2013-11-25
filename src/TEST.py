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

LOG_FILE_PREFIX = "TEST_"
LOG_FILE_SUFFIX = ".txt"

LOGGER_NAME = "phylo_logger"


## FUNCTIONS

# Test GUI functionality
def test_gui(logger):

    logger.debug("rig")

    return


# Test loading and parsing a FASTA file
def test_load_parse(logger):

    logger.debug("Testing FASTA file LOAD and PARSE functionality...")

    logger.debug("Filename: %s", TEST_FILE_1)
    logger.debug("Species: %s", TEST_FILE_1_SPECIES)
    
    file_reader = FASTAFile(TEST_FILE_1, TEST_FILE_1_SPECIES)

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

    return


# Test matching objects
def test_matching(logger):

    logger.debug("rig")

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

    return logger

# Run all tests
def test_project():

    logger = get_logger()

    test_load_parse(logger)

    return


if __name__ == '__main__':
    test_project()
