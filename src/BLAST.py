##
#
# FILE: BLAST.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains functionality for implementing
# gapped BLAST subsequence alignment using the blosum62
# matrix.
#

from imports import *
import Matching


# Default word length "w"
DEFAULT_WORD_LEN = 3
DEFAULT_SUB_MATRIX = blosum62

# Score threshold percentage - i.e. match must be in the X highest percent of scores
# evaluated thus far or it is discarded
DEFAULT_THRESHOLD = 5

# Class for holding general information about the application of a BLAST algorithm
class BLASTInfo (object):

    # Initialize - takes a list of FASTAFile objects, word length, substitution matrix, and threshold
    def __init__(self, matrix=DEFAULT_SUB_MATRIX, word_length=DEFAULT_WORD_LEN, 
        threshold=DEFAULT_THRESHOLD):

        self.matrix = matrix
        self.word_length = word_length
        self.threshold = threshold


# Class for implementing BLAST alignment between two Species objects
class BLASTSpeciesPair (object):

    # Initialize - takes two Species objects and an info class object
    def __init__(self, species_a, species_b, info=None):

        self.species_a = species_a
        self.species_b = species_b

        # Info struct
        if info != None and isinstance(info, BLASTInfo):
            self.info = info

# Class for containing BLAST alignment data for a given pair of Variants
class BLASTVariantPair (object):

    # Initialize - takes two Variant objects and an info class object and
    # scores their alignment against one another
    def __init__(self, variant_a, variant_b, info=None):

        self.variant_a = variant_a
        self.variant_b = variant_b

        # Info struct
        if info != None and isinstance(info, BLASTInfo):
            self.info = info

    # Ungapped BLAST alignment function - takes two sequences, a substitution
    # matrix (default blosum62), and a percentage of best matches to calculate score with
    def ungapped_BLAST(self, seq_a, seq_b, matrix=DEFAULT_SUB_MATRIX, 
        word_length=DEFAULT_WORD_LEN, threshold):


