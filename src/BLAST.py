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
from math import *
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

        # Info about BLAST
        self.matrix = matrix
        self.word_length = word_length
        self.threshold = threshold
        self.all_words = self.get_all_possible_words(matrix, word_length)

        # To dynamically adjust score threshold
        self.total_match_score = 0
        self.total_match_squared_score = 0
        self.total_number_matches = 0
        self.threshold_score = 0

    # Here we get a list of all possible word_length length words in matrix and store it
    def get_all_words(self, matrix, word_length):
        all_words = []
        count = 0
        return self.get_all_words_helper(all_words, matrix, count, word_length)

    # Recursive helper method for getting all words
    # Takes a word list, sub matrix, char count so far, word length
    def get_all_words_helper(self, word_list, matrix, count, word_length):
        # Base case
        if count == word_length:
            return word_list

        # Recursive case
        words_so_far = len(word_list)
        for i in range(words_so_far):
            # Add all new possible permutations
            for j in range(len(matrix)):
                new_word = word_list[i] + matrix[j]
                word_list.append(new_word)
            # Remove old base word
            word_list.pop(i)

        # Call again with new word list and longer length
        return self.get_all_words_helper(word_list, matrix, count + 1, word_length)

    # Re-calculate the approximated threshold score after each new matching
    def recalculate_threshold_score(self):

        # Get Z-score
        mean = total_match_score / total_number_matches
        std_dev = math.pow(((pow(self.total_match_squared_score), 2) / self.total_number_matches)
            - pow(mean, 2))

        # Get approximate threshold score
        self.threshold_score = (((100 - self.threshold) / 100) * self.total_number_matches) + 0.5

    # Score two words of length "word_length" using the given substitution matrix
    def BLAST_score(self, variant_a, variant_b, matrix, word_length):
        score = 0
        for i in range(word_length):
            score += matrix[variant_a.sequence[i]][variant_b.sequence[i]]
        return score
        

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

    # Get the word_length length words in each of the subsequences
    def get_possible_words(self):
        possible_words = []
        for i in range(len(self.variant_a) - self.word_length + 1):
            if 

    # Ungapped BLAST alignment function - takes two sequences, a substitution
    # matrix (default blosum62), and a percentage of best matches to calculate score with
    def ungapped_BLAST(self, seq_a, seq_b, matrix=DEFAULT_SUB_MATRIX, 
        word_length=DEFAULT_WORD_LEN, threshold):


